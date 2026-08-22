#pragma once

#include <cstddef>
#include <cstdio>
#include <cstdlib>

/* Every buffer in the program is sized from the input, so an allocation that
 * fails means the run cannot go on. There is nothing to fall back to and no
 * partial answer worth printing, so this reports and exits rather than
 * returning a null the caller would have to check at each of its uses.
 */
[[noreturn]] inline void out_of_memory(const char* what, std::size_t bytes)
{
    std::fprintf(stderr, "Out of memory: could not allocate %zu bytes for %s.\n", bytes, what);
    std::exit(1);
}
