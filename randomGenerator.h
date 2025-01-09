#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace K {

class RandomGenerator {
public:
    virtual ~RandomGenerator();
    inline uint32_t uniformU32() {
        if (m_index == m_size) fillBuffer();
        return m_buffer[m_index++];
    }
    inline float uniformF32() {
        constexpr float kNorm = 1.f/4294967296.f;
        return kNorm * uniformU32();
    }
    static std::unique_ptr<RandomGenerator> defaultGenerator(uint32_t sequence = 0, uint32_t bufferSize = 256);
protected:
    RandomGenerator(uint32_t bufferSize);
    virtual void fillBuffer() = 0;

    std::vector<uint32_t> m_buffer;
    uint32_t              m_size;
    uint32_t              m_index;
};

}
