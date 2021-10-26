#pragma once
#include <stdint.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>

namespace scapi{

template <typename D>
class PtrMatrix{
    public:
        PtrMatrix(size_t R, size_t C, D * mat) : R_(R), C_(C), mat_(mat){

        }

        const D & operator()(size_t i, size_t j) const {
            return mat_[i * C_ + j];
        }

        D & operator()(size_t i, size_t j) {
            return mat_[i * C_ + j];
        }

        size_t rows() const {
            return R_;
        }

        size_t cols() const {
            return R_;
        }

    private:
        size_t R_;
        size_t C_;
        D * mat_;

};

}
