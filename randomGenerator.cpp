#include "randomGenerator.h"

#include <random>

namespace K {

RandomGenerator::~RandomGenerator() = default;

RandomGenerator::RandomGenerator(uint32_t bufferSize) : m_size(bufferSize), m_index(bufferSize) {
    m_buffer.resize(m_size);
}

// See xoshiro256+ in https://en.wikipedia.org/wiki/Xorshift
class Xoshiro final : public RandomGenerator {
    uint64_t s[4];
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
public:
    Xoshiro(uint32_t seed, uint32_t bufferSize) : RandomGenerator(bufferSize) {
        std::mt19937_64 rndm(seed);
        for (int i = 0; i < 4; ++i) s[i] = rndm();
    }

protected:
    virtual void fillBuffer() override {
        for (uint32_t i = 0; i < m_size; ++i) {
            m_buffer[i] = s[0] + s[3];
            uint64_t t = s[1] << 17;
            s[2] ^= s[0];
            s[3] ^= s[1];
            s[1] ^= s[2];
            s[0] ^= s[3];
            s[2] ^= t;
            s[3] = rotl(s[3], 45);
        }
        m_index = 0;
    }
};

std::unique_ptr<RandomGenerator> RandomGenerator::defaultGenerator(uint32_t sequence, uint32_t bufferSize) {
    return std::make_unique<Xoshiro>(5489u+sequence, bufferSize);
}

}

#ifdef TEST

#include <chrono>
#include <cstdio>

int main() {
    auto rng = K::RandomGenerator::defaultGenerator();
    printf("First 10 random numbers:\n");
    for (int i = 0; i < 10; ++i) printf("    %f\n", rng->uniformF32());

    auto t1 = std::chrono::steady_clock::now();
    const int n = 1000000;
    float sum = 0;
    for (int j = 0; j < n; ++j) {
        sum += rng->uniformF32();
    }
    auto t2 = std::chrono::steady_clock::now();
    printf("Average: %f\n", sum/n);
    printf("Time: %g us for %d samples\n", 1e-3*std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count(), n);
    return 0;
}

#endif
