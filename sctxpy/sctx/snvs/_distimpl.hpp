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
#include "_cafimpl.hpp"

namespace scapi {

inline size_t condensed_idx(size_t i, size_t j, size_t N){
    if(i < j) std::swap(i, j);
    return N * j - j * (j + 1) / 2 + i - 1 - j;
}

struct AFSNVDist{
    AFSNVDist(size_t N, const spCS<uint32_t> * counts, double * da, uint32_t * dc)
        : counts(counts), DA(da), DC(dc), N(N)
    {
    
    }

    const spCS<uint32_t> * counts;
    double               * DA;
    uint32_t             * DC;
    size_t                 N;

    double dist_matrix(size_t c1, size_t c2) const {
        auto idx = condensed_idx(c1, c2, N);
        return DA[idx];
    }

    void dist(size_t c1, size_t c2) const {
        auto & cc = *counts;
        int32_t s1 = cc.indptr[c1];
        int32_t e1 = cc.indptr[c1 + 1];
        int32_t s2 = cc.indptr[c2];
        int32_t e2 = cc.indptr[c2 + 1];
        uint32_t count = 0;
        double dist = 0.0;

        while(s1 != e1 && s2 != e2){
            auto i1 = cc.indices[s1];
            auto i2 = cc.indices[s2];
            if(i1 == i2){
                double AF1 = 1.0 * cc.alt[s1]  / (cc.ref[s1] + cc.alt[s1]);
                double AF2 = 1.0 * cc.alt[s2]  / (cc.ref[s2] + cc.alt[s2]);
                count++;
                dist += std::abs(AF1 - AF2);
                s1++;
                s2++;
            }else if(i1 < i2){
                s1++;
            }else{
                s2++;
            }
        }

        if(count > 0){
            dist = dist / count;
        }else{
            dist = 1.0;
        }
        
        auto idx = condensed_idx(c1, c2, N);
        DC[idx] = count;
        DA[idx] = dist;
    }
};

template <typename DD>
class DistCalculator{
    public:
        DistCalculator(DD * data, unsigned int threads)
            : data_(data), threads_(threads)
        {

        }

        void run();
    private:
        struct DistWorker{
            DistWorker(unsigned int start, unsigned int end,  DD * data) 
            : data(*data), start(start), end(end)
            {

            }

            void operator()();

            void run(){
                thread = std::thread(std::ref(*this));
            }

            void join() {
                thread.join();
            }

            std::thread            thread;
            DD                   & data;
            size_t                 start = 0;
            size_t                 end = 0;
        };
        DD                    * data_;
        std::mutex              mutex_;
        unsigned int            threads_ = 1;
};

template <typename DD>
void DistCalculator<DD>::run() {
    std::vector<DistWorker*> workers;
    if(threads_ < 1) threads_ = 1;
    workers.resize(threads_);
    auto N = data_->N;
    unsigned int X = (N * (N - 1)) / 2;
    unsigned int step = X / threads_;
    unsigned int count = 0;
    std::vector<unsigned int> chunks = {0};
    //std::cout << "N = " << N << "\n";
    for(size_t i = 0; i < N; i++){
        unsigned int Y = N - i  - 1;
        if((Y + count) >= step){
            //std::cout << "i = " << i << " Y = " << Y << " count = " << count << "\n";
            chunks.push_back(i);
            count = 0;
        }
        count += Y;
    }
    chunks.back() = N;

    for(size_t i = 0; i < workers.size(); i++){
        workers[i] = new DistWorker(chunks[i], chunks[i + 1], data_);
        //std::cout << "Worker " << i << " range = " << chunks[i] << " - " << chunks[i + 1] << "\n";
    }

    //std::cout << "starting run\n";
    for(DistWorker * w : workers) w->run();
    //std::cout << "starting join\n";
    for(DistWorker * w : workers) w->join();
    //std::cout << "starting cleanup\n";
    for(DistWorker * w : workers) delete w;

    workers.clear();

}

template <typename DD>
void DistCalculator<DD>::DistWorker::operator()() {
    for(size_t c1 = start; c1 < end; c1++){
        for(size_t c2 = c1 + 1; c2 < data.N; c2++){
            data.dist(c1, c2);
        }
    }
}

template<typename DD>
inline void calcDist(DD * data, unsigned int threads){
    DistCalculator<DD> dc(data, threads);
    dc.run();
}

/*
inline void calcRanges(size_t N, const double * DA, const uint32_t * DC, double * mm){
    for(size_t i = 0; i < N; i++){
        double sum = 0;
        double dsum = 0;
        for(size_t j = 0; j < i; j++){
            // j is smaller than i
            size_t idx = condensed_idx(i, j, N);
            uint32_t c = DC[idx];
            dsum += DA[idx];
            sum += c;
        }
        for(size_t j = i + 1; j < N; j++){
            // j is bigger than i
            size_t idx = condensed_idx(i, j, N);
            uint32_t c = DC[idx];
            dsum += DA[idx];
            sum += c;
        }
        mm[i * 2] = dsum / (N - 1);
        mm[i * 2 + 1] = sum / (N - 1);
    }
}
*/

}
