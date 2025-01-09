#include "randomGenerator.h"

#include <cstdlib>
#include <chrono>
#include <cstdio>
#include <cmath>

template <typename Generator>
void checkPerformance(const char* msg, Generator& g, int nPoint) {
    auto t1 = std::chrono::steady_clock::now();
    double average_r2 = 0;
    for (int j = 0; j < nPoint; ++j) {
        auto [x, y] = g.generate();
        float r2 = x*x + y*y;
        average_r2 += r2;
    }
    auto t2 = std::chrono::steady_clock::now();
    printf("\n===================== %s(%d points)\n", msg, nPoint);
    printf("<r2> = %g\n", average_r2/nPoint);
    printf("Time: %g ms\n", 1e-6*std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
}

struct Rejection {
    Rejection(K::RandomGenerator& rng) : rng(rng) {}
    inline std::pair<float, float> generate() {
        float x, y;
        do {
            x = 2*rng.uniformF32();
            y = 2*rng.uniformF32();
        } while (x*x + y*y > 1);
        return std::make_pair(x, y);
    }
    K::RandomGenerator& rng;
};

inline std::pair<float, float> randomAzimuth(K::RandomGenerator& rng) {
    float x,x2,y,y2,r;
    do {
        x = 2*rng.uniformF32() - 1; x2 = x*x;
        y = rng.uniformF32(); y2 = y*y;
        r = x2 + y2;
    } while (r > 1);
    float ri = 1/r;
    return std::make_pair((x2 - y2)*ri, 2*x*y*ri);
}

struct Direct1 {
    Direct1(K::RandomGenerator& rng) : rng(rng) {}
    inline std::pair<float, float> generate() {
        float r = sqrt(rng.uniformF32());
        auto [cphi, sphi] = randomAzimuth(rng);
        return std::make_pair(r*cphi, r*sphi);
    }
    K::RandomGenerator& rng;
};

struct Direct2 {
    Direct2(K::RandomGenerator& rng) : rng(rng) {}
    inline std::pair<float, float> generate() {
        float r = std::max(rng.uniformF32(), rng.uniformF32());
        auto [cphi, sphi] = randomAzimuth(rng);
        return std::make_pair(r*cphi, r*sphi);
    }
    K::RandomGenerator& rng;
};

int main(int argc, char **argv) {
    int nPoint = argc > 1 ? atoi(argv[1]) : 1000000;

    auto rng = K::RandomGenerator::defaultGenerator();

    Rejection rej(*rng);
    checkPerformance("Rejection", rej, nPoint);

    Direct1 direct1(*rng);
    checkPerformance("Direct1", direct1, nPoint);

    Direct2 direct2(*rng);
    checkPerformance("Direct2", direct2, nPoint);

    return 0;
}
