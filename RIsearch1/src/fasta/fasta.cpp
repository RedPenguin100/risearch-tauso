/* A FASTA reader with one line of lookahead.
 *
 * FASTA has no end-of-record marker: a record is finished only once the next
 * record's header line has been read. That line is kept in the FASTAFILE, so
 * the next call starts from it instead of having to push it back. Keeping it
 * per file rather than in a static is what lets two files be read at once.
 *
 *     ffp = OpenFASTA(seqfile);
 *     while (ReadFASTA(ffp, seq, name)) { ... }
 *     CloseFASTA(ffp);
 *
 * OpenFASTA answers null if the file cannot be read or is empty. ReadFASTA
 * answers false once the file holds no further record.
 *
 * A line may be at most FASTA_MAXLINE bytes. A longer sequence line is read as
 * two, which changes nothing, but a longer header has its name truncated.
 *
 * The name is the first whitespace-separated token after the '>', found with
 * strtok, which is not reentrant -- one thread may read at a time.
 *
 * Derived from the FASTA reader in SRE's Bio5495/BME537 course material.
 */

#include "fasta.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "fasta/ResidueTable.h"
#include "memory/alloc.hpp"


namespace {

/* Besides residues, a sequence line may hold the line ending, the spaces and
 * digits of a coordinate column, and the '-' and '.' of an alignment gap. A
 * byte outside that set is corruption, and dropping it quietly would return a
 * shorter sequence that still looks valid.
 */
bool tolerated_in_sequence(unsigned char c)
{
    return c == '\n' || c == '\r' || c == ' ' || c == '\t' || c == '-' || c == '.' ||
           (c >= '0' && c <= '9');
}

/* Out of line so that fprintf does not sit inside the per-character loop, where
 * it costs about a tenth of the throughput for a branch a readable file never
 * takes.
 */
[[noreturn]] __attribute__((noinline, cold)) void reject_corrupt_byte(unsigned char c,
                                                                      const char* name)
{
    fprintf(stderr, "Corrupt byte 0x%02x in the sequence of '%s'.\n", c, name);
    exit(1);
}

} // namespace


FASTAFILE* OpenFASTA(const char* seqfile)
{
    FASTAFILE* ffp = reinterpret_cast<FASTAFILE*>(malloc(sizeof(FASTAFILE)));
    if (ffp == nullptr) {
        out_of_memory("a FASTA reader", sizeof(FASTAFILE));
    }

    if (strcmp(seqfile, "-")) {        /*returns 0/FALSE if they are same! */
        ffp->fp = fopen(seqfile, "r"); /* checked below */
    } else {
        ffp->fp = stdin;
    }
    if (ffp->fp == nullptr) {
        free(ffp);
        return nullptr;
    }
    /* Prime the lookahead with the first line, which every record starts from. */
    if (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp) == nullptr) {
        free(ffp);
        return nullptr;
    }
    return ffp;
}

bool ReadFASTA(FASTAFILE* ffp, ByteBuffer& ret_seq, ByteBuffer& ret_name)
{
    /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
     */
    if (ffp->buffer[0] != '>')
        return false;

    /* Parse out the name: the first non-whitespace token after the >
     */
    const char* s = strtok(ffp->buffer + 1, " \t\r\n");
    ret_name.clear();
    ret_name.append(s, strlen(s));
    ret_name.terminate();

    /* Everything else 'til the next descline is the sequence. clear() keeps the
     * capacity the previous record grew, so reading a file of similar records
     * settles on one buffer instead of allocating per record.
     */
    ret_seq.clear();
    for (;;) {
        if (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp) == NULL) {
            /* End of file. fgets leaves the buffer alone when it fails, so the
             * descline this record came from would still be sitting in it and
             * the next call would hand back this same record again, for ever.
             * Empty the lookahead so that call reports the file as finished.
             */
            ffp->buffer[0] = '\0';
            break;
        }
        if (ffp->buffer[0] == '>')
            break; /* a-ha, we've reached the next descline */

        /* A sequence line is normally residues followed by a line ending, so
         * take the leading run in one append and only go character by character
         * over whatever follows it.
         */
        const unsigned char* line = reinterpret_cast<const unsigned char*>(ffp->buffer);
        std::size_t run = 0;
        while (kResidue.is[line[run]])
            run++;
        ret_seq.append(ffp->buffer, run);

        for (std::size_t i = run; line[i] != '\0'; i++) {
            const unsigned char c = line[i];
            if (!kResidue.is[c]) {
                if (!tolerated_in_sequence(c))
                    reject_corrupt_byte(c, ret_name.c_str());
                continue; /* accept any alphabetic character */
            }
            ret_seq.push_back(static_cast<char>(c));
        }
    }
    ret_seq.terminate();
    return true;
}

void CloseFASTA(FASTAFILE* ffp)
{
    fclose(ffp->fp);
    free(ffp);
}
