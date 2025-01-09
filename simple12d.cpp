#include "randomGenerator.h"

#include <cstdlib>
#include <chrono>
#include <cstdio>

inline int sample_point(K::RandomGenerator* rng, float * pos, float& R2) {
    int ntry = 0;
    float r2;
    do {
        ++ntry;
        r2 = 0;
        for (int j = 0; j < 12; ++j) {
            pos[j] = 2*rng->uniformF32() - 1;
            r2 += pos[j]*pos[j];
        }
    } while (r2 > 1);
    R2 = r2;
    return ntry;
}

int main(int argc, char **argv) {
    int n_sample = argc > 1 ? atoi(argv[1]) : 100000;

    auto rng = K::RandomGenerator::defaultGenerator();
    float pos[12];
    uint64_t ntry = 0;
    auto t1 = std::chrono::steady_clock::now();
    double average_r2 = 0;
    for (int j = 0; j < n_sample; ++j) {
        float r2;
        ntry += sample_point(rng.get(), pos, r2);
        average_r2 += r2;
    }
    auto t2 = std::chrono::steady_clock::now();
    printf("It took %zu attempts to sample %d points in a 12d sphere\n", ntry, n_sample);
    printf("Sampling efficiency: %g\n", 1.*n_sample/ntry);
    printf("<r^2> = %g\n", average_r2/n_sample);
    printf("Run time: %g ms\n", 1e-3*std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
    return 0;
}
