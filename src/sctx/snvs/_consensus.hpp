#pragma once
#include <stdint.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <numeric>
#include <fstream>
#include <cmath>
#include <limits>
#include "_cafimpl.hpp"
#include "_matrix.hpp"

namespace scapi{

void build_genotypes(const spCS<uint32_t> * counts, int32_t * labels, size_t G, uint32_t * out, double * summaries){
    auto & cc = *counts;
    //std::cout << "N = " << cc.N << " M = " << cc.M << " R = " << cc.R << " G = " << G << "\n";
    PtrMatrix<uint32_t> umat(cc.R, G * 3, out);
    PtrMatrix<double> dmat(cc.R, G * 4, summaries);
    for(size_t i = 0; i < cc.R; i++){
        int32_t s = cc.indptr[i];
        int32_t e = cc.indptr[i + 1];
        for(int32_t j = s; j < e; j++){
            uint32_t ci = cc.indices[j];
            int32_t l = labels[ci];
            if(l < 0 || l > 1) continue;
            umat(i, l * 3) += cc.ref[j];
            umat(i, l * 3 + 1) += cc.alt[j];
            umat(i, l * 3 + 2) += (cc.alt[j] > 0) | (cc.ref[j] > 0);
        }
    }

    std::vector<std::vector<double>> AFs(G);
    std::vector<double> tots(G);
    for(size_t i = 0; i < cc.R; i++){
        for(size_t l = 0; l < G; l++){
            //dmat(i, l * 2) = umat(i, l * 3 + 2) > 0 ? (1.0 * umat(i, l * 3 + 1) / (umat(i, l * 3) + umat(i, l * 3 + 1))) : 0.0;
            AFs[l].clear();
            tots[l] = 0.0;
        }
        int32_t s = cc.indptr[i];
        int32_t e = cc.indptr[i + 1];
        for(int32_t j = s; j < e; j++){
            uint32_t ci = cc.indices[j];
            int32_t l = labels[ci];
            if(l < 0 || l > 1) continue;
            double AF = 1.0 * cc.alt[j] / (cc.alt[j] + cc.ref[j]);
            tots[l] += cc.alt[j] + cc.ref[j];
            AFs[l].push_back(AF);
        }
        for(size_t l = 0; l < G; l++){
            if(!AFs[l].empty()){
                double mean = std::accumulate(AFs[l].begin(), AFs[l].end(), 0.0) / AFs[l].size();
                double rmse = 0.0;
                double mae = 0.0;
                for(auto a : AFs[l]){
                    rmse += (a - mean) * (a - mean);
                    mae += std::abs(a - mean);
                }
                dmat(i, l * 4) = mean;
                dmat(i, l * 4 + 1) = std::sqrt(rmse / AFs[l].size());
                dmat(i, l * 4 + 2) = mae / AFs[l].size();
                dmat(i, l * 4 + 3) = tots[l] / AFs[l].size();
            }
        }
    }
}

void cell_div(const spCS<uint32_t> * counts, const uint32_t * ridx, size_t N, const double * summary_AF, double * div, uint32_t * cnts){
    auto & cc = *counts;
    std::cout << "N = " << cc.N << " M = " << cc.M << " R = " << cc.R << " ridx N = " << N << "\n";
    for(size_t ii = 0; ii < N; ii++){
        size_t i = ridx[ii];
        int32_t s = cc.indptr[i];
        int32_t e = cc.indptr[i + 1];
        double saf = summary_AF[i];
        for(int32_t j = s; j < e; j++){
            uint32_t ci = cc.indices[j];
            double AF = 1.0 * cc.alt[j] / (cc.alt[j] + cc.ref[j]);
            cnts[ci]++;
            div[ci] += std::abs(saf - AF);
        }
    }
}

};
