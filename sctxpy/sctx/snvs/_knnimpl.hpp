#pragma once
#include <stdint.h>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>
#include "_cafimpl.hpp"
#include "_distimpl.hpp"

namespace scapi {

template <typename DD>
class FindKNN{
    public:
        FindKNN(const DD * data, uint32_t * knni, double * knnd, uint8_t * bad_cells, unsigned int K, unsigned int threads)
            : data_(data), knni_(knni), knnd_(knnd), bad_cells_(bad_cells), K_(K), threads_(threads) 
        {

        }

        void run();

    private:

        struct KNNWorker{
            KNNWorker(unsigned int start, unsigned int end, const DD * data, unsigned int K) 
            : data(data), start(start), end(end), K(K)
            {

            }

            void set_out(uint32_t * ki, double * kd, uint8_t * bc){
                knni = ki;
                knnd = kd;
                bad_cells = bc;
            }


            void operator()();

            void run(){
                thread = std::thread(std::ref(*this));
            }

            void join() {
                thread.join();
            }

            std::thread      thread;
            const DD       * data;
            uint32_t       * knni = nullptr;
            double         * knnd = nullptr;
            uint8_t        * bad_cells = nullptr;
            size_t           start = 0;
            size_t           end = 0;
            unsigned int     K;
        };


        std::mutex        mutex_;
        const DD        * data_;
        uint32_t        * knni_ = nullptr;
        double          * knnd_ = nullptr;
        uint8_t         * bad_cells_ = nullptr;
        unsigned int      K_ = 5;
        unsigned int      threads_ = 1;
};

/*
inline size_t condensed_idx(size_t i, size_t j, size_t N){
    assert(j < i);
    return N * j - j * (j + 1) / 2 + i - 1 - j;
}
*/

template <typename DD>
inline void FindKNN<DD>::KNNWorker::operator()(){
    std::vector<double> dists;
    std::vector<uint32_t> index;
    std::vector<uint32_t> ai;
    for(size_t c1 = start; c1 != end; c1++){
        dists.clear();
        index.clear();
        ai.clear();
        for(size_t c2 = 0; c2 < c1; c2++){
            double dist = data->dist_matrix(c1, c2);
            if(dist < 1.0){
                dists.push_back(dist);
                index.push_back(c2);
            }
        }
        for(size_t c2 = c1 + 1; c2 < data->N; c2++){
            double dist = data->dist_matrix(c1, c2);
            if(dist < 1.0){
                dists.push_back(dist);
                index.push_back(c2);
            }
        }
        ai.resize(index.size());
        std::iota(ai.begin(), ai.end(), 0);

        if(index.size() < K){
            bad_cells[c1] = 1;
            continue;
        }

        std::partial_sort(ai.begin(), ai.begin() + std::min(ai.size(), (size_t)K), ai.end(),
           [&dists](size_t i1, size_t i2) {return dists[i1] < dists[i2];});

        size_t cidx = c1 * K;
        for(size_t i = 0; i < std::min(index.size(), (size_t)K); i++){
            knni[cidx + i] = index[ai[i]];
            knnd[cidx + i] = dists[ai[i]];
        }
    }
}

template <typename DD>
inline void FindKNN<DD>::run(){
    std::vector<KNNWorker*> workers;
    if(threads_ < 1) threads_ = 1;
    workers.resize(threads_);
    size_t N = data_->N;
    unsigned int step = N / threads_;
    //std::cout << "step = " << step << "\n";
    unsigned int start = 0;
    for(size_t i = 0; i < workers.size(); i++){
        unsigned int end = start;
        if(i == (workers.size() - 1)){
            end = N;
        }else{
            end += step;
        }
        //std::cout << "Workers[" << i <<"] start = " << start << " end = " << end << "\n";
        workers[i] = new KNNWorker(start, end, data_, K_);
        workers[i]->set_out(knni_, knnd_, bad_cells_);
        start = end;
    }

    for(KNNWorker * w : workers) w->run();
    for(KNNWorker * w : workers) w->join();
    for(KNNWorker * w : workers) delete w;

    workers.clear();
}

template <typename DD>
inline void findKNN(DD * data, uint32_t * knni, double * knnd, uint8_t * bad_cells, unsigned int K, unsigned int threads){
    FindKNN<DD> dd(data, knni, knnd, bad_cells, K, threads);
    dd.run();
}

}
