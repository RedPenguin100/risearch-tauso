#pragma once

#include "fasta.h"

/* Owns a FASTAFILE for the length of a scope. handle() is null if the file
   could not be opened, which is what a caller checks before reading. */
class FastaRAII {
public:
    explicit FastaRAII(const char* sequence_file);
    ~FastaRAII();

    FastaRAII(const FastaRAII&) = delete;
    FastaRAII(FastaRAII&&) = delete;
    FastaRAII& operator=(const FastaRAII&) = delete;
    FastaRAII& operator=(FastaRAII&&) = delete;

    FASTAFILE* handle() const;

private:
    FASTAFILE* m_fasta_file_handle;
};
