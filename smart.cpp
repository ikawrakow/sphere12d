#include "randomGenerator.h"

#include <cstdlib>
#include <chrono>
#include <cstdio>
#include <algorithm>
#include <cmath>

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

inline int sample_point(K::RandomGenerator* rng, float * pos) {
    float z[6];
    for (int j = 0; j < 6; ++j) z[j] = rng->uniformF32();
    std::sort(z, z+6);
    float last = 0;
    for (int j = 0; j < 6; ++j) {
        float r = sqrt(z[j] - last);
        last = z[j];
        auto [cphi, sphi] = randomAzimuth(*rng);
        pos[2*j + 0] = r*cphi;
        pos[2*j + 1] = r*sphi;
    }
    return 1;
}

int main(int argc, char **argv) {
    int n_sample = argc > 1 ? atoi(argv[1]) : 100000;

    auto rng = K::RandomGenerator::defaultGenerator();
    float pos[12];
    uint64_t ntry = 0;
    auto t1 = std::chrono::steady_clock::now();
    double average_r2 = 0;
    for (int j = 0; j < n_sample; ++j) {
        ntry += sample_point(rng.get(), pos);
        float r2 = 0;
        for (int j = 0; j < 12; ++j) r2 += pos[j]*pos[j];
        average_r2 += r2;
    }
    auto t2 = std::chrono::steady_clock::now();
    printf("It took %zu attempts to sample %d points in a 12d sphere\n", ntry, n_sample);
    printf("Sampling efficiency: %g\n", 1.*n_sample/ntry);
    printf("<r^2> = %g\n", average_r2/n_sample);
    printf("Run time: %g ms\n", 1e-3*std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    return 0;
}
