#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>

#include <cstdint>
#include <cstdio>
#include <cstring>

#include "InteractionAlignment.h"
#include "MatrixStore.h"
#include "RunningMax.h"
#include "cli/cli.h"
#include "energy.hpp"
#include "memory/ByteBuffer.hpp"
#include "memory/MallocRAII.hpp"
#include "operations.h"
#include "optimization/QueryProfile.h"
#include "traceback.h"


/* Turns a winning DP cell into printed output.
 *
 * RIs_linSpace finds where the best alignments end but keeps only two rows, so
 * it cannot say which nucleotides paired. This clips a window of at most tblen
 * around a hit, re-aligns that window with a full matrix to recover the pairing,
 * converts the score to an energy and prints it.
 *
 * It owns the scratch that re-alignment needs, so those buffers are allocated
 * once per query/target pair and reused for every hit -- there may be hundreds
 * of thousands of them.
 */
class HitReporter {
public:
    HitReporter(const unsigned char* query, const unsigned char* target, std::uint32_t n,
                const QueryProfile<std::int32_t>& profile, const config_st& config,
                const char* qname, const char* tname)
        : m_query(query), m_target(target), m_n(n), m_profile(profile), m_config(config),
          m_qname(qname), m_tname(tname), m_reference(reference_from_matrix(config.mat_name)),
          m_matrices(config.tblen), m_best(config.tblen + 1),
          m_hit(static_cast<int>(1.5 * config.tblen)),
          /* Seven decimal fields, the energy and the separators, plus the two
             names: everything in a line whose length is known up front. */
          m_line_fixed(128 + std::strlen(qname) + std::strlen(tname))
    {
        m_line.reserve(m_line_fixed);
    }

    /* pos_i, pos_j: the DP cell the hit ends at. score: what the sweep recorded
       for it, which is what gets printed -- see the note in print(). */
    void report(std::uint32_t pos_i, std::uint32_t pos_j, int score, bool is_suboptimal)
    {
        const auto w = extract(pos_i, pos_j);

        RIs(m_query + w.qbeg - 1, w.target, w.qlen, w.tlen, &m_hit, m_config, m_matrices.M(),
            m_matrices.Ix(), m_matrices.Iy(), m_profile, w.qbeg - 1, m_best.get());

        const auto energy =
            static_cast<double>(m_hit.max + m_config.extension_penalty * m_hit.nucleotide_count() -
                                m_reference) /
            (-100.0);

        if (energy <= m_config.max_energy) {
            print(w, pos_j, score, energy, is_suboptimal);
        }
    }

    /* The runs the sweep leaves behind -- the best score ending at each target
       position and the query position it ended at -- turned into reported hits.
       Which of them get printed is an output mode question. */
    /* The runs are read by value, so a sweep that keeps them narrower than the
       reporting does needs no widening pass of its own. */
    template<typename Score>
    void report_sweep(const Score* run_scores, const Score* run_positions, int threshold,
                      const RunningMax& running_max)
    {
        report_top_hit_if_needed(running_max);

        if (m_config.doSubopt) {
            auto j = static_cast<int>(m_n); /* the runs are 0-based over rows 1..n */
            while (j--) {
                if (run_scores[j] <= threshold) {
                    continue;
                }
                const auto window_bottom = j - window_size(j);
                const auto best = best_row_in_window(run_scores, j);

                report(run_positions[best], best + 1, run_scores[best], true);

                /* Resume below the whole window, not just below the row reported. */
                j = window_bottom;
            }
        }

        finalize_report();
    }

    /* The same walk, over a run given as one bit per target position rather than
       as scores.
     *
     * Almost no row of a run clears the threshold, and the batched sweep leaves
     * its scores sixteen queries to a row, so reading a whole run to find the few
     * that matter reads sixteen times what it needs. One bit per target position
     * is a sixteenth of that, and a word of them none of which is set skips
     * thirty-two rows at once.
     *
     * The threshold is already in the bits, so no score is compared again.
     *
     * Only for a vicinity of zero: a window takes the best of a reported row's
     * neighbours, which reads the score of a row that cleared nothing.
     */
    template<typename Score>
    void report_sweep_sparse(const std::uint32_t* clears, const Score* scores,
                             const Score* positions, std::size_t stride,
                             const RunningMax& running_max)
    {
        report_top_hit_if_needed(running_max);

        if (m_config.doSubopt) {
            const auto rows = static_cast<std::uint32_t>(m_n);
            const auto words = (rows + 31) / 32;
            for (auto w = words; w-- > 0;) {
                auto bits = clears[w];
                while (bits) {
                    /* Highest row of the word first, as the dense walk descends. */
                    const auto b = 31u - static_cast<unsigned>(__builtin_clz(bits));
                    bits &= ~(1u << b);
                    const auto j = w * 32 + b;
                    if (j >= rows) {
                        continue;
                    }
                    const auto at = static_cast<std::size_t>(j) * stride;
                    report(positions[at], static_cast<int>(j) + 1, scores[at], true);
                }
            }
        }

        finalize_report();
    }

    int hitcount() const
    {
        return m_hitcount;
    }

private:
    /* How far below row j a vicinity window reaches, which is `vicinity` rows
       unless the run ends first. */
    int window_size(int j) const
    {
        return m_config.vicinity < j ? m_config.vicinity : j;
    }

    /* The best row of the window that ends at row j. Ties go to the lowest row:
       the walk starts at the bottom of the window and only a strictly better
       row displaces what it holds. */
    template<typename Score>
    int best_row_in_window(const Score* run_scores, int j) const
    {
        auto best = j;
        for (auto row = j - window_size(j); row < j; row++) {
            if (run_scores[row] > run_scores[best]) {
                best = row;
            }
        }
        return best;
    }

    void report_top_hit_if_needed(const RunningMax& running_max)
    {
        if (!(m_config.doSubopt && (m_config.filter_e || m_config.printShort > 1))) {
            report(running_max.pos_i, running_max.pos_j, running_max.score, false);
        }
    }
    void finalize_report()
    {
        /* One line per query and target rather than per hit. */
        if (m_config.printShort == 3) {
            write_line(FMT_COMPILE("{}\t{}\t{}\n"), m_qname, m_tname, m_hitcount);
        }
    }

    struct Window {
        std::uint32_t qbeg; /* 1-based query position the window starts at */
        std::uint32_t qlen;
        std::uint32_t tlen;
        ReversedSequence target;
    };

    /* Clip to at most tblen back from the hit. Both slices are read where they
       already sit; the target's runs backwards, which its view accounts for. */
    Window extract(std::uint32_t pos_i, std::uint32_t pos_j)
    {
        const auto qbeg = pos_i > m_config.tblen - 1 ? pos_i - (m_config.tblen - 1) : 1;
        const auto tbeg = pos_j > m_config.tblen - 1 ? pos_j - (m_config.tblen - 1) : 1;
        return Window{qbeg, pos_i - qbeg + 1, pos_j - tbeg + 1,
                      ReversedSequence{m_target + m_n - tbeg}};
    }

    /* Formats one line into the reporter's buffer and writes it.
     *
     * A run prints one of these per hit and there can be hundreds of thousands
     * of them, so the format strings are compile-time and the buffer is sized
     * once rather than per line. `extra` is the length of any variable-length
     * field in the line -- the interaction string -- since the two names are the
     * only other unbounded part and they are known when the reporter is built.
     *
     * The line is written as it is finished, so anything else that prints stays
     * in order with it.
     */
    /* A line whose length the reserve made at construction already covers. */
    template<typename Format, typename... Args>
    void write_line(const Format& format, const Args&... args)
    {
        write_line(std::size_t{0}, format, args...);
    }

    /* extra is what the one unbounded field -- the interaction string -- adds
       beyond that reserve. */
    template<typename Format, typename... Args>
    void write_line(std::size_t extra, const Format& format, const Args&... args)
    {
        m_line.reserve(m_line_fixed + extra);
        char* const end = fmt::format_to(m_line.data(), format, args...);
        std::fwrite(m_line.data(), 1, static_cast<std::size_t>(end - m_line.data()), stdout);
    }

    void print(const Window& w, int pos_j, int score, double energy, bool is_suboptimal)
    {
        const auto qb = w.qbeg + m_hit.qbeg - 1;
        const auto qe = w.qbeg + m_hit.qend - 1;
        const auto tb = m_n - pos_j + m_hit.tbeg;
        const auto te = m_n - pos_j + m_hit.tend;

        switch (m_config.printShort) {
        case 1:
            /* The two callers have always disagreed here: the best hit prints
               five fields as (tbeg, tend); a suboptimal prints six, as
               (tend, tbeg) plus the interaction string. Preserved so that
               extracting this function cannot move any output. */
            if (is_suboptimal) {
                const char* const ia = m_hit.ali_ia.get();
                write_line(std::strlen(ia), FMT_COMPILE("{}\t{}\t{}\t{}\t{:.2f}\t{}\n"), qb, qe, te,
                           tb, energy, ia);
            } else {
                write_line(FMT_COMPILE("{}\t{}\t{}\t{}\t{:.2f}\n"), qb, qe, tb, te, energy);
            }
            break;

        case 2:
            write_line(FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n"), m_qname, qb, qe,
                       m_tname, tb, te, score, energy);
            break;

        case 3:
            m_hitcount += 1;
            break;

        default:
            /* score comes from the linear-space sweep, energy from the window
               re-alignment above. They can disagree when the alignment is longer
               than tblen, which is the open question the two TODOs marked. */
            printf("Free energy [kcal/mol]: %.2f (%d)\n", energy, score);
            printf("%d - %d\n", qb, qe); /* alignment in seq1 from to */
            printf("%s\n%s\n%s\n", m_hit.ali_seq1.get(), m_hit.ali_ia.get(), m_hit.ali_seq2.get());
            printf("%d - %d (3' <-- 5')\n", te, tb);
            break;
        }
    }

    const unsigned char* m_query;
    const unsigned char* m_target;
    std::uint32_t m_n;
    const QueryProfile<std::int32_t>& m_profile;
    const config_st& m_config;
    const char* m_qname;
    const char* m_tname;
    float m_reference;

    /* Scratch, reused across every hit. */
    MatrixStore m_matrices;
    /* Best M + close per query column; transpose_best_cell reads it. */
    MallocRAII<std::int32_t> m_best;
    IA m_hit;

    std::size_t m_line_fixed;
    ByteBuffer m_line;

    int m_hitcount = 0;
};