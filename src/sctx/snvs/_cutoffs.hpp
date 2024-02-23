#pragma once
#include <stdint.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iomanip>
#include "_distimpl.hpp"

namespace scapi {

/*
inline void calcCutoffs(size_t N, size_t C, const uint32_t * DC, const uint32_t * cutoffs, uint32_t * out){
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < i; j++){
            // j is smaller than i
            size_t idx = condensed_idx(i, j, N);
            uint32_t c = DC[idx];
            for(size_t x = 0; x < C; x++){
                if(c >= cutoffs[x] && c <= cutoffs[x + 1]){
                    out[i * C + x]++;
                }
            }
        }
        for(size_t j = i + 1; j < N; j++){
            // j is bigger than i
            size_t idx = condensed_idx(i, j, N);
            uint32_t c = DC[idx];
            for(size_t x = 0; x < C; x++){
                if(c >= cutoffs[x] && c <= cutoffs[x + 1]){
                    out[i * C + x]++;
                }
            }
        }
    }
}
*/

inline void calc2dCutoffs(size_t N, size_t C1, size_t C2, const uint32_t * DC, const uint32_t * count_co, const double * cell_co, double * out){
    std::vector<uint32_t> counts;
    for(size_t i = 0; i < N; i++){
        counts.clear();
        counts.resize(C1);
        for(size_t j = 0; j < N; j++){
            if(i == j) continue;
            size_t idx = condensed_idx(i, j, N);
            uint32_t c = DC[idx];
            for(size_t x = 0; (x < C1 && c >= count_co[x]); x++){
                counts[x]++;
            }
        }

        for(size_t x = 0; x < C1; x++){
            double p = 1.0 * counts[x] / N;
            for(size_t y = 0; (y < C2 && p >= cell_co[y]); y++){
                out[x * C2 + y]++;
            }
        }
    }
    for(size_t x = 0; x < C1; x++){
        for(size_t y = 0; y < C2; y++){
            out[x * C2 + y] /= (N - 1);
        }
    }
}

inline void applyCutoffs(size_t N, double * DA, uint32_t * DC, uint32_t count_co, double cell_co, uint32_t * cidx, uint32_t * ret){
    size_t x = 0;
    std::vector<bool> check(N);
    for(size_t i = 0; i < N; i++){
        size_t cc = 0;
        for(size_t j = 0; j < N; j++){
            if(i == j) continue;
            size_t idx = condensed_idx(i, j, N);
            uint32_t c = DC[idx];
            if(c >= count_co){
                cc++;
            }
        }
        double perc = 1.0 * cc / (N - 1);
        if(perc >= cell_co){
            cidx[x++] = i;
            check[i] = true;
        }
    }
    ret[0] = x;
    size_t cc = 0;
    for(size_t i = 0; i < N; i++){
        for(size_t j = i + 1; j < N; j++){
            if(check[i] && check[j]){
                size_t idx = condensed_idx(i, j, N);
                DA[cc] = DA[idx];
                DC[cc] = DC[idx];
                if(DC[cc] < count_co){
                    DA[cc] = 1.0;
                }
                cc++;
            }
        }
    }
    ret[1]= cc;
}

}
