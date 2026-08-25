/***********************************************************
  RIsearch v 1.2   --   RNA-RNA interaction search
  Copyright 2012 Anne Wenzel <wenzel@rth.dk> (RIsearch v.1.0 and v.1.1)
  Copyright 2021 Giulia I Corsi <giulia@rth.dk> (Extension of RIsearch v.1.1 in RIsearch v.1.2)

  This file is part of RIsearch.

  RIsearch is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.

***********************************************************/

#include <malloc.h>
#include <unistd.h>

#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "FastaRAII.h"
#include "FastaRecord.h"
#include "align/QueryBatch.h"
#include "align/dispatch.h"
#include "cli/cli.h"
#include "dsm.h"
#include "nucleotide.h"
#include "pair_header.h"

void tune_glibc_allocator()
{
    /* Set the minimum amount of memory that gets given back to the OS, so a run
       page-faults less when it allocates and runs faster. This does not raise
       memory usage to the set amount -- the process just does not drop back
       below what it has already spiked to.

       Both lines are needed: setting only the trim threshold freezes the size
       above which an allocation is served by mmap instead of the heap, and the
       buffers a pair reuses then come from fresh pages every time, which faults
       more than leaving the allocator alone. */
#if defined(__GLIBC__)
    mallopt(M_TRIM_THRESHOLD, 64 * 1024 * 1024);
    mallopt(M_MMAP_THRESHOLD, 64 * 1024 * 1024);
#endif
}

/* One sequence ready to align. id is which record of its file the sequence was,
   counting the ones seq2ix rejected, so a query keeps the number it is reported
   under. A sequence given on the command line has no record, which is what id 0
   says, and it is named from_cli. */
struct SequenceItem {
    std::string name;
    ByteBuffer indices;
    int id{0};
};

std::vector<SequenceItem> load_sequences_from_cli(const char* cli_str, const char* role)
{
    SequenceItem item;
    item.name = "from_cli";
    if (!seq2ix(static_cast<std::uint32_t>(std::strlen(cli_str)), cli_str, item.indices,
                "from command line", role)) {
        exit(-1); /* non-alpha char in input -- no other sequence to move on to */
    }

    if (item.indices.is_empty()) {
        std::fprintf(stderr, "No %s seq in the one given on the command line\n", role);
        exit(-1);
    }

    std::vector<SequenceItem> sequences;
    sequences.push_back(std::move(item));
    return sequences;
}

std::vector<SequenceItem> load_sequences_from_file(const char* file_path, const char* role)
{
    FastaRAII fasta_sequences(file_path);
    if (!fasta_sequences.handle()) {
        /* seq2ix names the role in lower case in its own messages; this one has
           always begun with a capital. */
        std::fprintf(stderr, "%c%s file %s is not readable\n",
                     std::toupper(static_cast<unsigned char>(role[0])), role + 1, file_path);
        exit(-1);
    }

    std::vector<SequenceItem> sequences;
    FastaRecord seq_record;
    int count = 0;
    while (seq_record.read(fasta_sequences.handle())) {
        count++;
        ByteBuffer indices;
        if (seq2ix(seq_record.get_size(), seq_record.get_sequence(), indices, seq_record.get_name(),
                   role)) {
            sequences.push_back({seq_record.get_name(), std::move(indices), count});
        }
    }

    if (sequences.empty()) {
        std::fprintf(stderr, "No %s seq in %s\n", role, file_path);
        exit(-1);
    }
    return sequences;
}

std::vector<SequenceItem> load_sequences(const char* file_path, const char* cli_str,
                                         const char* role)
{
    if (file_path) {
        return load_sequences_from_file(file_path, role);
    }

    if (cli_str) {
        return load_sequences_from_cli(cli_str, role);
    }

    /* getArgs refuses a run with no query or no target, so this is unreachable
       -- alternative run seq against itself!? */
    std::fprintf(stderr, "No %s seq given!", role);
    exit(-1);
}

void print_header(const SequenceItem& query, const SequenceItem& target)
{
    print_pair_header(query.name.c_str(), query.id,
                      static_cast<std::uint32_t>(query.indices.size()), target.name.c_str(),
                      target.id, static_cast<std::uint32_t>(target.indices.size()));
}

void process_target(const SequenceItem& target, const std::vector<SequenceItem>& queries, Dsm& dsm,
                    const config_st& config, QueryProfileCache& profiles)
{
    /* Queries are swept together where the sweep is what runs. Only a file
       holds enough queries to fill a batch's lanes, and only a file numbers the
       records that pairing them off reads. */
    const bool from_files = target.id != 0 && queries.front().id != 0;

    // force start not supported for batching
    const bool batching = from_files && !uses_force_start(config);

    QueryBatch batch;

    for (const auto& query : queries) {
        /* Pairing a query with the target of the same number is only meaningful
           where both are numbered. */
        if (from_files && !config.all_vs_all && target.id != query.id) {
            continue;
        }
        const auto len_seq1 = static_cast<std::uint32_t>(query.indices.size());

        if (batching) {
            batch.add(query.indices, query.name.c_str(), query.id, len_seq1);
            if (batch.full()) {
                batch.run(target.indices, dsm, target.name.c_str(), target.id, config,
                          profiles);
            }
            continue;
        }

        if (config.printShort < 2) {
            print_header(query, target);
        }

        run_alignment(query.indices, target.indices, dsm, query.name.c_str(), target.name.c_str(),
                      config);
    }

    if (batching && !batch.empty()) {
        batch.run(target.indices, dsm, target.name.c_str(), target.id, config, profiles);
    }
}

int main(int argc, char* argv[])
{
    tune_glibc_allocator();

    /* values filled in by getArgs from the command line */
    static config_st config;

    Dsm dsm;

    getArgs(argc, argv, config);

    getMat(config.mat_name, &dsm[0][0][0][0], config.extension_penalty,
           config.transpose_matrix_flag);

    const auto queries = load_sequences(config.seq1_file_name, config.seq1_cli, "query");
    const auto targets = load_sequences(config.seq2_file_name, config.seq2_cli, "target");

    /* Built here rather than per target: a query's profile does not depend on
       the target it is swept against. */
    QueryProfileCache profiles(dsm);

    for (const auto& target : targets) {
        process_target(target, queries, dsm, config, profiles);
    }

    return 0;
}
