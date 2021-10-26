#cython: c_string_type=str, c_string_encoding=ascii, language_level=3str
from libc.stdint cimport int32_t, uint64_t, uint32_t, uint8_t
from libcpp cimport bool
from libcpp.string cimport string

cdef extern from './_cafimpl.hpp' namespace 'scapi':
    cdef cppclass spCS[T]:
        spCS(Py_ssize_t N, Py_ssize_t M, Py_ssize_t R, T * ref, T * alt, int32_t * indptr, int32_t * indices) 

    void writeMM[X](const spCS[X] & data, const string & ref_out, const string & alt_out)

cdef extern from './_knnimpl.hpp' namespace 'scapi':
    void findKNN[T](T * counts, uint32_t * knni, double * knnd, uint8_t * bad_cells, unsigned int K, unsigned int threads)

cdef extern from './_distimpl.hpp' namespace 'scapi':
    cdef cppclass AFSNVDist:
        AFSNVDist(Py_ssize_t N, const spCS[uint32_t] * counts, double * DA, uint32_t * DC)
        const spCS[uint32_t] * counts;
        double               * DA;
        uint32_t             * DC;

    void calcDist[T](T * counts, unsigned int threads)
    #void calcRanges(Py_ssize_t N, const double * DA, const uint32_t * DC, double * mm)

cdef extern from './_consensus.hpp' namespace 'scapi':
    void build_genotypes(const spCS[uint32_t] * counts, int32_t * labels, Py_ssize_t G, uint32_t * out, double * summary)
    void cell_div(const spCS[uint32_t] * counts, const uint32_t * ridx, size_t N, const double * summary_AF, double * div, uint32_t * cnt)

cdef extern from './_cutoffs.hpp' namespace 'scapi':
    void calc2dCutoffs(size_t N, size_t C1, size_t C2, const uint32_t * DC, const uint32_t * count_co, const double * cell_co, double * out)
    void applyCutoffs(size_t N, double * DA, uint32_t * DC, uint32_t count_co, double cell_co, uint32_t * cidx, uint32_t * ret)
