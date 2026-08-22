#pragma once

#include <cstdint>
#include <cstring>

#include "align/HitReporter.h"
#include "align/ScoreTargetBatched.h"
#include "align/dispatch.h"
#include "align/first_row.h"
#include "align/int16_safety.h"
#include "align/optimization/BatchedQueryProfile.h"
#include "cli/cli.h"
#include "memory/ByteBuffer.hpp"
#include "memory/GrowableBuffer.hpp"

/**
 * Queries held back so that several can be swept against one target together.
 *
 * A query and its header are printed from here rather than from the read loop,
 * so that collecting queries before running them cannot move either. What comes
 * out is what the same queries produce one at a time, in the same order.
 */
class QueryBatch {
public:
    static constexpr unsigned kQueries = BatchedQueryProfile::kLanes;

    /**
     * The sequence and the name are copied: the read loop reuses its buffers for
     * the next record before this batch runs.
     */
    void add(const ByteBuffer& query_seq, const char* name, int query_count, std::uint32_t len)
    {
        Entry& e = m_entries[m_count];
        e.seq.clear();
        e.seq.append(query_seq.data(), query_seq.size());
        e.name.clear();
        e.name.append(name, std::strlen(name));
        e.name.terminate();
        e.query_count = query_count;
        e.len = len;
        m_count++;
    }

    bool full() const
    {
        return m_count == kQueries;
    }

    bool empty() const
    {
        return m_count == 0;
    }

    void run(const ByteBuffer& target_seq, Dsm& dsm, const char* tname, int target_count,
             const config_st& config)
    {
        m_exact_rows = config.doSubopt && config.vicinity > 0;
        const bool swept = sweep_impl(target_seq, dsm, config.min_score);

        for (auto k = 0u; k < m_count; k++) {
            const Entry& e = m_entries[k];
            if (config.printShort < 2) {
                printf("\n\nquery %d: %s (%u nts) vs. target %d: %s (%u nts)\n\n", e.query_count,
                       e.name.data(), e.len, target_count, tname,
                       static_cast<std::uint32_t>(target_seq.size()));
            }
            if (swept) {
                report_query(k, e, target_seq, dsm, tname, config);
            } else {
                run_alignment(e.seq, target_seq, dsm, e.name.data(), tname, config);
            }
        }
        clear();
    }

    void clear()
    {
        m_count = 0;
    }

    /* Sweeps the batch against one target, or answers false and leaves the
       queries to be aligned one at a time. Public so that what it produces can
       be checked query by query against the ordinary sweep. */
    bool sweep(const ByteBuffer& target_seq, Dsm& dsm, int threshold)
    {
        return sweep_impl(target_seq, dsm, threshold);
    }

    /* Valid after a sweep that answered true: one query's whole run, the best
       score ending at each target position and the query position it ended at. */
    const std::int16_t* scores(unsigned query) const
    {
        return m_hs.get() + static_cast<std::size_t>(query) * m_n;
    }

    const std::int16_t* positions(unsigned query) const
    {
        return m_hp.get() + static_cast<std::size_t>(query) * m_n;
    }

    int swept_length() const
    {
        return m_n;
    }

private:
    struct Entry {
        ByteBuffer seq;
        ByteBuffer name;
        int query_count = 0;
        std::uint32_t len = 0;
    };

    /* Below five queries the batch is behind the ordinary sweep -- it pays for
       sixteen lanes whatever it fills -- so a handful of ASOs, or the last few
       records of a file, go the ordinary way. */
    static constexpr unsigned kMinQueries = 5;


    /* What the kernel reads of the queries, gathered once. m is the longest of
       them, which is what the rows are sized on. */
    struct BatchInputs {
        const unsigned char* queries[kQueries];
        std::uint32_t lengths[kQueries];
        std::uint32_t m = 0;
    };

    /* Whether the batched kernel can take this batch at all, and what it reads
       of it. Everything here is decided before any of the run is set up, so a
       no costs nothing but the walk over the queries.

       The whole batch goes through the kernel or none of it does: one query the
       int16 bound does not hold for would have to be swept on its own anyway. */
    bool can_sweep(const ByteBuffer& target_seq, const Dsm& dsm, BatchInputs& inputs) const
    {
        if (!CPU_HAS_AVX2 || m_count < kMinQueries || target_seq.is_empty()) {
            return false;
        }

        for (auto k = 0u; k < m_count; k++) {
            const Entry& e = m_entries[k];
            if (e.seq.is_empty()) {
                return false;
            }
            inputs.queries[k] = e.seq.unsigned_data();
            inputs.lengths[k] = static_cast<std::uint32_t>(e.seq.size());
            if (!fits_int16(dsm, inputs.queries[k], inputs.lengths[k])) {
                return false;
            }
            inputs.m = MAX(inputs.m, inputs.lengths[k]);
        }

        /* The column a row max was reached at is carried in a lane of shorts. */
        return inputs.m <= 32000;
    }

    void allocate_sweep_buffers(std::size_t m, std::size_t n)
    {
        const auto query_stride = (m + 1) * kQueries;
        const auto target_stride = n * kQueries;

        m_m_rows.reserve(2 * query_stride);
        m_iy_rows.reserve(2 * query_stride);
        m_scan1.reserve(query_stride);
        m_hs16.reserve(target_stride);
        m_hp16.reserve(target_stride);
        /* The runs laid out by query are what a vicinity window reads; without
           one the reporting takes the sweep's own layout and these are never
           touched. At a long target they are the two largest buffers here. */
        if (m_exact_rows) {
            m_hs.reserve(target_stride);
            m_hp.reserve(target_stride);
        }
    }

    bool sweep_impl(const ByteBuffer& target_seq, Dsm& dsm, int threshold)
    {
#if RISEARCH1_HAS_AVX2
        BatchInputs inputs;
        if (!can_sweep(target_seq, dsm, inputs)) {
            return false;
        }

        const auto m = inputs.m;
        const auto n = static_cast<int>(target_seq.size());
        const unsigned char* const target = target_seq.unsigned_data();
        m_threshold = threshold;

        if (!m_profile.build(inputs.queries, inputs.lengths, m_count, dsm, m)) {
            return false;
        }

        allocate_sweep_buffers(m, n);

        std::int16_t first_score[kQueries];
        std::int16_t first_pos[kQueries];
        if (!init_first_row(inputs.queries, inputs.lengths, m, target[n - 1], dsm, first_score,
                            first_pos)) {
            return false;
        }
        for (auto lane = 0u; lane < kQueries; lane++) {
            m_hs16[lane] = first_score[lane];
            m_hp16[lane] = first_pos[lane];
        }

        run_sweep(target, n, first_score, first_pos);
        if (m_exact_rows) {
            deinterleave(n);
        } else {
            build_clear_bits(n, threshold);
        }
        m_n = n;
        return true;
#else
        (void)target_seq;
        (void)dsm;
        (void)threshold;
        return false;
#endif
    }


#if RISEARCH1_HAS_AVX2
    /* Target position 1 follows the recurrence that governs a first row, not the
       one the sweep states, so it is computed here in int32 and handed over: the
       rows as values, and its Ix as the terms that make the sweep's scan
       reproduce it. */
    bool init_first_row(const unsigned char* const* queries, const std::uint32_t* lengths,
                        std::uint32_t m, unsigned char t_last, Dsm& dsm, std::int16_t* first_score,
                        std::int16_t* first_pos)
    {
        const auto stride = static_cast<std::size_t>(m + 1) * kQueries;
        std::int16_t* const m_row = m_m_rows.get() + stride; /* the row read first */
        std::int16_t* const iy_row = m_iy_rows.get() + stride;
        std::int16_t* const scan_row = m_scan1.get();
        const std::int16_t* const ix_prefix = m_profile.ix_prefix();

        /* Column 0 is never read; give it a value so nothing is undefined. */
        for (auto lane = 0u; lane < kQueries; lane++) {
            m_m_rows[lane] = NEG_INF_SHORT;
            m_iy_rows[lane] = 0;
        }

        m_row1.reserve(3 * static_cast<std::size_t>(m + 1));
        int* const M = m_row1.get() + 0 * (m + 1);
        int* const Ix = m_row1.get() + 1 * (m + 1);
        int* const Iy = m_row1.get() + 2 * (m + 1);

        for (auto lane = 0u; lane < kQueries; lane++) {
            if (lane >= m_count) {
                for (auto i = 0u; i <= m; i++) {
                    m_row[i * kQueries + lane] = NEG_INF_SHORT;
                    iy_row[i * kQueries + lane] = NEG_INF_SHORT;
                    scan_row[i * kQueries + lane] = NEG_INF_SHORT;
                }
                first_score[lane] = NEG_INF_SHORT;
                first_pos[lane] = 1;
                continue;
            }

            const unsigned char* const q = queries[lane];
            const auto len = lengths[lane];

            const RunningVectorMax running_row_max =
                score_first_row<int>(M, Ix, Iy, q, len, t_last, dsm);

            for (auto i = 0u; i <= len; i++) {
                m_row[i * kQueries + lane] =
                    static_cast<std::int16_t>(MAX(M[i], static_cast<int>(NEG_INF_SHORT)));
                iy_row[i * kQueries + lane] =
                    static_cast<std::int16_t>(MAX(Iy[i], static_cast<int>(NEG_INF_SHORT)));
            }
            /* Column 1's step is never taken: the scan enters column 2 with the
               magic value Ix[1][1] is. */
            scan_row[0 * kQueries + lane] = NEG_INF_SHORT;
            scan_row[1 * kQueries + lane] = NEG_INF_SHORT;
            for (auto i = 2u; i <= len; i++) {
                const int term = Ix[i] - ix_prefix[i * kQueries + lane] - M[i - 1];
                if (term < -BatchedQueryProfile::kTermLimit ||
                    term > BatchedQueryProfile::kTermLimit) {
                    return false;
                }
                scan_row[i * kQueries + lane] = static_cast<std::int16_t>(term);
            }
            /* Past this query nothing is reachable, and every term the sweep
               reads there is the magic value too. */
            for (auto i = len + 1; i <= m; i++) {
                m_row[i * kQueries + lane] = NEG_INF_SHORT;
                iy_row[i * kQueries + lane] = NEG_INF_SHORT;
                scan_row[i * kQueries + lane] = NEG_INF_SHORT;
            }

            first_score[lane] = static_cast<std::int16_t>(running_row_max.score);
            first_pos[lane] = static_cast<std::int16_t>(running_row_max.pos_i);
        }
        return true;
    }

    /* One bit per query per target position: does that query's row clear the
       threshold? Held per query so a run is read at a sixteenth of what the
       scores are, with the OR over each group of thirty-two rows beside it so a
       group none of them clears is one test. */
    __attribute__((target("avx2"))) void build_clear_bits(int n, int threshold)
    {
        const auto words = (static_cast<std::size_t>(n) + 31) / 32;
        m_clears.reserve(words * kQueries);

        std::memset(m_clears.get(), 0, words * kQueries * sizeof(std::uint32_t));


        const __m256i thr = v_int_to_avx2<std::int16_t>(static_cast<std::int16_t>(
            threshold > SHRT_MAX ? SHRT_MAX : (threshold < SHRT_MIN ? SHRT_MIN : threshold)));

        for (auto j = 0; j < n; j++) {
            const __m256i v = v_vec_load(m_hs16.get() + static_cast<std::size_t>(j) * kQueries);
            /* One bit per query, out of a sixteen lane compare. */
            auto bits = v_lane_bits16(v_greater_than<std::int16_t>(v, thr));
            while (bits) {
                const auto q = static_cast<unsigned>(__builtin_ctz(bits));
                bits &= bits - 1;
                const auto base = static_cast<std::size_t>(q) * words;
                m_clears[base + j / 32] |= 1u << (j % 32);
            }
        }
        m_clear_words = words;
    }

    /* The sweep writes a target position's sixteen scores together; the
       reporting reads one query's whole run, so the two want opposite layouts. */
    __attribute__((target("avx2"))) void deinterleave(int n)
    {
        std::int16_t* hs_query[kQueries];
        std::int16_t* hp_query[kQueries];
        for (auto k = 0u; k < m_count; k++) {
            hs_query[k] = m_hs.get() + static_cast<std::size_t>(k) * n;
            hp_query[k] = m_hp.get() + static_cast<std::size_t>(k) * n;
        }
        const std::int16_t* hs16 = m_hs16.get();
        const std::int16_t* hp16 = m_hp16.get();

        /* j stays in scope so the tail can resume where the blocks stopped: the
           target positions past the last whole sixteen go one at a time. */
        auto j = 0;
        for (; j + static_cast<int>(kQueries) <= n; j += kQueries) {
            v_transpose16x16(hs16, hs_query, m_count);
            v_transpose16x16(hp16, hp_query, m_count);
            hs16 += kQueries * kQueries;
            hp16 += kQueries * kQueries;
            for (auto k = 0u; k < m_count; k++) {
                hs_query[k] += kQueries;
                hp_query[k] += kQueries;
            }
        }
        for (; j < n; j++) {
            for (auto k = 0u; k < m_count; k++) {
                *hs_query[k]++ = hs16[k];
                *hp_query[k]++ = hp16[k];
            }
            hs16 += kQueries;
            hp16 += kQueries;
        }
    }

    __attribute__((target("avx2"))) void run_sweep(const unsigned char* target, int n,
                                                   const std::int16_t* first_score,
                                                   const std::int16_t* first_pos)
    {
        BatchedRunningMax best;
        best.score = v_vec_load(first_score);
        best.pos_i = v_vec_load(first_pos);
        best.pos_j_lo = v_int_to_avx2<std::int32_t>(1);
        best.pos_j_hi = v_int_to_avx2<std::int32_t>(1);

        const auto stride = static_cast<std::size_t>(m_profile.m() + 1) * kQueries;
        std::int16_t* const M[2] = {m_m_rows.get(), m_m_rows.get() + stride};
        std::int16_t* const Iy[2] = {m_iy_rows.get(), m_iy_rows.get() + stride};

        /* A vicinity window takes the best of a reported row's neighbours, so
           where one is asked for, every row's score is read and none can be left
           standing as a bound. */
        if (m_exact_rows) {
            score_target_batched<false>(target, m_profile, M, Iy, m_scan1.get(), m_hs16.get(),
                                        m_hp16.get(), static_cast<std::size_t>(n), m_threshold,
                                        best);
        } else {
            score_target_batched<true>(target, m_profile, M, Iy, m_scan1.get(), m_hs16.get(),
                                       m_hp16.get(), static_cast<std::size_t>(n), m_threshold,
                                       best);
        }

        v_vec_store(m_best_score, best.score);
        v_vec_store(m_best_i, best.pos_i);
        v_vec_store(m_best_j, best.pos_j_lo);
        v_vec_store(m_best_j + 8, best.pos_j_hi);
    }
#endif

    /* One query's hits, reported exactly as the unbatched path reports them: the
       sweep's runs are widened into the int arrays HitReporter reads, and the
       traceback that follows is the int32 one. */
    void report_query(unsigned lane, const Entry& e, const ByteBuffer& target_seq, Dsm& dsm,
                      const char* tname, const config_st& config)
    {
        const auto m = static_cast<std::uint32_t>(e.seq.size());
        const auto* query = e.seq.unsigned_data();
        const auto* target = target_seq.unsigned_data();

        const QueryProfile<std::int32_t> profile(query, m, dsm, has_positive_gap(dsm));

        RunningMax running_max{};
        running_max.set(m_best_score[lane], m_best_i[lane], m_best_j[lane]);

        HitReporter reporter(query, target, m_n, profile, config, e.name.data(), tname);
        if (m_exact_rows) {
            const auto offset = static_cast<std::size_t>(lane) * m_n;
            reporter.report_sweep(m_hs.get() + offset, m_hp.get() + offset, config.min_score,
                                  running_max);
        } else {
            reporter.report_sweep_sparse(m_clears.get() + lane * m_clear_words, m_hs16.get() + lane,
                                         m_hp16.get() + lane, kQueries, running_max);
        }
    }

    Entry m_entries[kQueries];
    unsigned m_count = 0;

    BatchedQueryProfile m_profile;
    GrowableBuffer<std::int16_t> m_m_rows, m_iy_rows, m_scan1;
    GrowableBuffer<std::int16_t> m_hs16, m_hp16;
    GrowableBuffer<std::int16_t> m_hs, m_hp;
    GrowableBuffer<int> m_row1;
    GrowableBuffer<std::uint32_t> m_clears;
    std::size_t m_clear_words = 0;
    std::int16_t m_best_score[kQueries]{};
    std::int16_t m_best_i[kQueries]{};
    int m_best_j[kQueries]{};
    int m_n = 0;
    int m_threshold = 0;
    bool m_exact_rows = true;
};
