#pragma once
#include <stdint.h>
#include <iostream>
#include <fstream>

namespace scapi{

template <typename T>
struct spCS{
    spCS() {

    }

    spCS(size_t N, size_t M, size_t R, T * ref, T * alt, int32_t * indptr, int32_t * indices) 
        : N(N), M(M), R(R), ref(ref), alt(alt), indptr(indptr), indices(indices)
    {

    }

    size_t    N = 0;
    size_t    M = 0;
    size_t    R = 0;
    T       * ref = NULL;
    T       * alt = NULL;
    int32_t * indptr = NULL;
    int32_t * indices = NULL;
};


template <typename X>
inline void writeMM(const spCS<X> & ref, const std::string & ref_out, const std::string & alt_out){
    std::ofstream rout(ref_out.c_str());
    std::ofstream aout(alt_out.c_str());

    std::cout << "Alt out = " << alt_out << " Ref out = " << ref_out << "\n";

    rout << "%%MatrixMarket matrix coordinate integer general\n%\n";
    rout << ref.N << "\t" << ref.M << "\t" << ref.indptr[ref.R] << "\n";
    aout << "%%MatrixMarket matrix coordinate integer general\n%\n";
    aout << ref.N << "\t" << ref.M << "\t" << ref.indptr[ref.R] << "\n";

    for(size_t i = 0; i < ref.R; i++){
        int32_t s1 = ref.indptr[i];
        int32_t e1 = ref.indptr[i + 1];
        while(s1 != e1){
            auto i1 = ref.indices[s1];
            rout << (i + 1) << " " << (i1 + 1) << " " << ref.ref[s1] << "\n";
            aout << (i + 1) << " " << (i1 + 1) << " " << ref.alt[s1] << "\n";
            s1++;
        }
    }


    /*
                fp.write('%%MatrixMarket matrix coordinate integer general\n')
                fp.write('%\n')
                fp.write(f'{mat.shape[0]} {mat.shape[1]} {mat.shape[0] * mat.shape[1]}\n')
    */

}

}
