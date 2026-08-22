#pragma once

#include <cstddef>
#include <cstdlib>

#include "memory/alloc.hpp"

/**
 * An array that grows to whatever a run needs and is reused by the next one.
 *
 * Sized rather than resized: growing throws away what the buffer held, because
 * every caller fills it before reading it. Holding the capacity here is what
 * separates this from MallocRAII, whose callers had to carry a matching
 * capacity of their own beside every buffer.
 *
 * ALIGNED, BECAUSE THE BATCHED SWEEP READS WHOLE REGISTERS AT REGISTER-SIZED
 * STRIDES. Its rows, runs and term tables are all indexed in units of sixteen
 * lanes, so every access is a 32-byte piece 32 bytes along from the last, and
 * the whole run is aligned exactly when the buffer is. malloc promises only 16:
 * half the small buffers and every large one then start 16 bytes into a step,
 * and each access is split in two. Measured at a 30% difference in the sweep.
 *
 * A cache line rather than the 32 bytes needed, which costs at most one line
 * per buffer and leaves nothing to revisit if a wider register is ever used.
 *
 * The single-query sweep gains nothing from this and does not use it: it reads
 * its diagonal predecessor at m_last + i - 1, one element off whatever the last
 * access was, so no choice of base makes those accesses aligned.
 */
template<typename T>
class GrowableBuffer {
public:
    static constexpr std::size_t kAlign = 64;

    GrowableBuffer() = default;

    ~GrowableBuffer()
    {
        std::free(m_data);
    }

    GrowableBuffer(const GrowableBuffer&) = delete;
    GrowableBuffer& operator=(const GrowableBuffer&) = delete;

    void reserve(std::size_t wanted)
    {
        if (wanted <= m_capacity) {
            return;
        }
        /* aligned_alloc wants a size that is a whole number of alignments. */
        const auto bytes = (wanted * sizeof(T) + kAlign - 1) / kAlign * kAlign;
        std::free(m_data);
        m_data = static_cast<T*>(std::aligned_alloc(kAlign, bytes));
        if (m_data == nullptr) {
            out_of_memory("a sweep buffer", bytes);
        }
        m_capacity = bytes / sizeof(T);
    }

    T* get() const
    {
        return m_data;
    }

    T& operator[](std::size_t i) const
    {
        return m_data[i];
    }

private:
    T* m_data = nullptr;
    std::size_t m_capacity = 0;
};
