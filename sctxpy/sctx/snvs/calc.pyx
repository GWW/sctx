#cython: c_string_type=str, c_string_encoding=ascii, language_level=3str

#Copyright (c) 2018-2020 Gavin W. Wilson

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


import numpy as NP
cimport numpy as NP
import cython
cimport _caf as ccaf
from scipy.sparse import csr_matrix
from libc.stdint cimport int32_t, uint64_t, uint32_t

NP.import_array()

cdef ccaf.spCS[NP.uint32_t] * _DCS(N, M, NP.ndarray[NP.uint32_t] ref, NP.ndarray[NP.uint32_t] alt, NP.ndarray[NP.int32_t] indptr, NP.ndarray[NP.int32_t] indices):
    cdef Py_ssize_t R = indptr.shape[0]
    cdef ccaf.spCS[NP.uint32_t] * ret = new ccaf.spCS[NP.uint32_t](N, M, R - 1, &ref[0], &alt[0], &indptr[0], &indices[0])
    #print("N = ",N)
    #print("M = ",M)
    #print("R = ",R)
    return ret

def write_MTX(ref, alt, rout, aout):
    ref.data = ref.data.astype('uint32')
    alt.data = alt.data.astype('uint32')
    cdef const ccaf.spCS[NP.uint32_t] * a =_DCS(ref.shape[0], ref.shape[1], ref.data, alt.data, ref.indptr, ref.indices)
    ccaf.writeMM(cython.operator.dereference(a), rout, aout)
    del a

def calculate_AF_KNN(K, NP.ndarray[NP.double_t] DA, NP.ndarray[NP.uint32_t] DC, threads=1, ascending=False):
    #Actual distance
    #print("Calculating KNN")
    cdef Py_ssize_t N = int(NP.ceil(NP.sqrt(DA.shape[0] * 2)))
    if N * (N - 1) != DA.shape[0] * 2:
        raise ValueError('Incompatible distance matrix it must be condensed')
    cdef NP.ndarray[NP.uint32_t, ndim=2] Ki = NP.full((N, K), N, dtype='uint32')
    cdef NP.ndarray[NP.double_t, ndim=2] Kd = NP.full((N, K), N, dtype='double')
    cdef NP.ndarray[NP.uint8_t, ndim=1]  bad = NP.full(N, 0, dtype='uint8')
    cdef ccaf.AFSNVDist * dist = new ccaf.AFSNVDist(N, NULL, &DA[0], &DC[0])
    ccaf.findKNN[ccaf.AFSNVDist](dist, &Ki[0, 0], &Kd[0, 0], &bad[0], K, threads)
    del dist
    return Ki, Kd

def calculate_cutoffs_2d(NP.ndarray[NP.uint32_t] DC, NP.ndarray[NP.uint32_t] count_cutoffs, NP.ndarray[NP.double_t] cell_cutoffs):
    N = int(NP.ceil(NP.sqrt(DC.shape[0] * 2)))
    if N * (N - 1) != DC.shape[0] * 2:
        raise ValueError('Incompatible distance matrix it must be condensed')
    C1 = len(count_cutoffs)
    C2 = len(cell_cutoffs)
    cdef NP.ndarray[NP.double_t, ndim=2] out = NP.full((C1, C2), 0, dtype='double')
    ccaf.calc2dCutoffs(N, C1, C2,  <const uint32_t*>&DC[0], <const uint32_t*>&count_cutoffs[0], <const double*>&cell_cutoffs[0], &out[0, 0])
    return out

def apply_cutoffs_to_matrices(NP.ndarray[NP.double_t] DA, NP.ndarray[NP.uint32_t] DC, count_cutoff, cell_cutoff):
    N = int(NP.ceil(NP.sqrt(DC.shape[0] * 2)))
    if N * (N - 1) != DC.shape[0] * 2:
        raise ValueError('Incompatible distance matrix it must be condensed')
    cdef NP.ndarray[NP.uint32_t, ndim=1] cidx = NP.zeros(N, dtype='uint32')
    cdef NP.ndarray[NP.uint32_t, ndim=1] ret = NP.zeros(2, dtype='uint32')
    
    #void applyCutoffs(size_t N, double * DA, uint32_t * DC, uint32_t count_co, double cell_co, uint32_t * cidx, uint32_t * c1, uint32_t * c2)
    ccaf.applyCutoffs(N, &DA[0], &DC[0], count_cutoff, cell_cutoff, &cidx[0], &ret[0])
    DA = DA[:ret[1]]
    DC = DC[:ret[1]]
    cidx = cidx[:ret[0]]
    return DA, DC, cidx

def calculate_AF_dist(refs, alts, threads=1):
    #print("Making CSC array and output arrays")
    ref = refs.tocsc()
    alt = alts.tocsc()
    ref.data = ref.data.astype('uint32')
    alt.data = alt.data.astype('uint32')
    X = ref.shape[1] * (ref.shape[1] - 1) // 2
    #print("M = ",ref.shape[1], "X = ",X)
    cdef const ccaf.spCS[NP.uint32_t] * a =_DCS(ref.shape[0], ref.shape[1], ref.data, alt.data, ref.indptr, ref.indices)
    cdef NP.ndarray[NP.double_t] DA = NP.full(X, -1.0, dtype='double')
    cdef NP.ndarray[NP.uint32_t] DC = NP.zeros(X, dtype='uint32')
    cdef Py_ssize_t N = refs.shape[1]
    cdef ccaf.AFSNVDist * dist = new ccaf.AFSNVDist(N, a, &DA[0], &DC[0])
    ccaf.calcDist[ccaf.AFSNVDist](dist, threads)
    del a
    del dist
    return DA, DC

def calculate_genotypes(ref, alt, NP.ndarray[NP.int32_t] labels):
    LS = set(labels)
    cdef Py_ssize_t G = sum(1 for l in LS if l > -1 and l < 2)
    #print('G = ', G, 'LS = ',LS)
    ref.data = ref.data.astype('uint32')
    alt.data = alt.data.astype('uint32')
    cdef const ccaf.spCS[NP.uint32_t] * a =_DCS(ref.shape[0], ref.shape[1], ref.data, alt.data, ref.indptr, ref.indices)
    cdef NP.ndarray[NP.uint32_t, ndim=2] out = NP.zeros((ref.shape[0], G * 3), dtype='uint32')
    cdef NP.ndarray[NP.double_t, ndim=2] out2 = NP.zeros((ref.shape[0], G * 4), dtype='double')
    ccaf.build_genotypes(a, &labels[0], G, &out[0, 0], &out2[0, 0])
    del a

    data = {}
    cols = []
    dt = []
    for i in range(0, G):
        data[f'ref_{i}'] = out[:,i * 3 + 0].copy()
        data[f'alt_{i}'] = out[:,i * 3 + 1].copy()
        data[f'cells_{i}'] = out[:,i * 3 + 2].copy()
        #data[f'total_AF_{i}'] = out2[:,i * 2 + 0].copy()
        data[f'mean_AF_{i}'] = out2[:,i * 4].copy()
        data[f'mean_AF_RMSE_{i}'] = out2[:,i * 4 + 1].copy()
        data[f'mean_AF_MAE_{i}'] = out2[:,i * 4 + 2].copy()
        data[f'mean_count_{i}'] = out2[:,i * 4 + 3].copy()
    return data

def calculate_cell_div(ref, alt, ridx, NP.ndarray[NP.double_t] AFs):
    cdef const ccaf.spCS[NP.uint32_t] * a =_DCS(ref.shape[0], ref.shape[1], ref.data, alt.data, ref.indptr, ref.indices)
    cdef NP.ndarray[NP.uint32_t, ndim=1] aridx = ridx.astype('uint32')
    cdef NP.ndarray[NP.double_t, ndim=1] div = NP.zeros(ref.shape[1], dtype='double')
    cdef NP.ndarray[NP.uint32_t, ndim=1] counts = NP.zeros(ref.shape[1], dtype='uint32')
    ccaf.cell_div(a, <const uint32_t*>(&aridx[0]), aridx.shape[0], <const double*>(&AFs[0]), &div[0], &counts[0])
    del a
    return div, counts
