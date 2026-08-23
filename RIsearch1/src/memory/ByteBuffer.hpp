#pragma once

#include <cstdlib>
#include <cstring>
#include <utility> // for std::exchange

#include "memory/alloc.hpp"

/**
 * Lightweight wrapper for char / byte pointers.
 *
 *
 * std::string is not good, we use this instead, faster and smaller binary size
 *
 * std::string needs to include locale and a lot of legacy stuff that make it heavy
 * std::string needs to deal with bunch of unrealistic edge cases which make it slow
 */
class ByteBuffer {
public:
    ByteBuffer() = default;

    ~ByteBuffer()
    {
        std::free(m_data);
    }

    ByteBuffer(const ByteBuffer&) = delete;
    ByteBuffer& operator=(const ByteBuffer&) = delete;


    // Move Constructor
    ByteBuffer(ByteBuffer&& other) noexcept
        : m_data(std::exchange(other.m_data, nullptr)), m_size(std::exchange(other.m_size, 0)),
          m_capacity(std::exchange(other.m_capacity, 0))
    {
    }

    // Move Assignment Operator
    ByteBuffer& operator=(ByteBuffer&& other) noexcept
    {
        if (this != &other) {
            std::free(m_data);
            m_data = std::exchange(other.m_data, nullptr);
            m_size = std::exchange(other.m_size, 0);
            m_capacity = std::exchange(other.m_capacity, 0);
        }
        return *this;
    }


    void clear()
    {
        m_size = 0;
    }

    [[nodiscard]] std::size_t size() const
    {
        return m_size;
    }
    [[nodiscard]] bool is_empty() const
    {
        return m_size == 0;
    }

    /* Valid only after terminate(). Empty rather than null before anything has
     * been written, so callers never have to check.
     */
    [[nodiscard]] const char* c_str() const
    {
        return m_data ? m_data : "";
    }

    /* Unchecked, like the array it replaces. */
    [[nodiscard]] char& operator[](std::size_t i)
    {
        return m_data[i];
    }
    [[nodiscard]] const char& operator[](std::size_t i) const
    {
        return m_data[i];
    }

    /* Null while the buffer is empty. Use this rather than &buffer[0] to take a
     * base pointer: indexing an empty buffer binds a reference to null.
     */
    [[nodiscard]] char* data()
    {
        return m_data;
    }
    [[nodiscard]] const char* data() const
    {
        return m_data;
    }


    [[nodiscard]] unsigned char* unsigned_data()
    {
        return reinterpret_cast<unsigned char*>(m_data);
    }
    [[nodiscard]] const unsigned char* unsigned_data() const
    {
        return reinterpret_cast<const unsigned char*>(m_data);
    }

    void append(const char* bytes, std::size_t count)
    {
        /* A blank line appends nothing, and memcpy may not be handed a null
         * destination even for zero bytes.
         */
        if (count == 0)
            return;
        reserve(m_size + count);
        std::memcpy(m_data + m_size, bytes, count);
        m_size += count;
    }

    void push_back(char byte)
    {
        reserve(m_size + 1);
        m_data[m_size++] = byte;
    }

    /* Writes the terminator without counting it, so c_str() can be handed to
     * anything expecting a C string while size() stays the byte count.
     */
    void terminate()
    {
        reserve(m_size + 1);
        m_data[m_size] = '\0';
    }


    void reserve(std::size_t wanted)
    {
        if (wanted <= m_capacity)
            return;
        std::size_t capacity = m_capacity ? m_capacity : 128;
        while (capacity < wanted)
            capacity *= 2;
        char* const grown = static_cast<char*>(std::realloc(m_data, capacity));
        if (grown == nullptr) {
            out_of_memory("a sequence buffer", capacity);
        }
        m_data = grown;
        m_capacity = capacity;
    }

private:
    char* m_data = nullptr;
    std::size_t m_size = 0;
    std::size_t m_capacity = 0;
};
