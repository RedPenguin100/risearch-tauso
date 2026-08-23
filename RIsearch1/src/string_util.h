#pragma once

/* Reverses str_beg[0 .. last] in place. last is the index of the last
   character, not the length. */
inline void reverse_inplace(char* str_beg, int last)
{
    char* str_end = &str_beg[last];

    while (str_end > str_beg) {
        const char tmp = *str_beg;
        *str_beg++ = *str_end;
        *str_end-- = tmp;
    }
}
