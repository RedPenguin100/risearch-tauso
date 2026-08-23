#pragma once

#include <cstdint>
#include <cstdio>

/* The banner printed before a query/target pair's hits.
 *
 * Four spellings, one for each way the pair was given: a sequence read from a
 * file prints its record number and its name, one given on the command line
 * prints neither. A record number of 0 is what says a sequence came from the
 * command line.
 *
 * Both the one-at-a-time path and the batched one print through this, so
 * holding queries back to fill a batch cannot move the banner.
 */
inline void print_pair_header(const char* qname, int qid, std::uint32_t qlen, const char* tname,
                              int tid, std::uint32_t tlen)
{
    if (qid == 0 && tid == 0) {
        std::printf("\n\nquery from_cli (%u nts) vs. target from_cli (%u nts)\n\n", qlen, tlen);
    } else if (qid == 0) {
        std::printf("\n\nquery from_cli (%u nts) vs. target %s (%u nts)\n\n", qlen, tname, tlen);
    } else if (tid == 0) {
        std::printf("\n\nquery %s (%u nts) vs. target from_cli (%u nts)\n\n", qname, qlen, tlen);
    } else {
        std::printf("\n\nquery %d: %s (%u nts) vs. target %d: %s (%u nts)\n\n", qid, qname, qlen,
                    tid, tname, tlen);
    }
}
