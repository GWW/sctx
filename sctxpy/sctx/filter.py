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
from .plt import format_subplots
import pylab as plt


def snv_qa_plot(ref, alt, out=None, close_fig=False):
    cell_counts = []
    alt_counts = []
    BB_counts = []
    AFs = []
    for s, e in zip(ref.indptr, ref.indptr[1:]):
        d1 = ref.data[s:e]
        d2 = alt.data[s:e]
        rsum = d1.sum()
        asum = d2.sum()
        AFs.append(100.0 * asum / (rsum + asum))
        alt_counts.append(NP.count_nonzero(d2))
        cell_counts.append(len(d2))
        BB_counts.append(((d2 > 0) & (d1 == 0)).sum() / len(d1)) 
    cc = NP.array(cell_counts)
    ac = NP.array(alt_counts)
    bc = 100.0 * NP.array(BB_counts)
    ac[ac > 100] = 100
    cc[cc > 100] = 100
    fig, axs = format_subplots(2, 2, bsize=4)
    axs = axs.flatten()
    fig.subplots_adjust(hspace=0.4)
    axs[0].hist(cc, bins=NP.linspace(0, 100, 101), weights=NP.ones(len(cc)) / len(cc), cumulative=True)
    axs[1].hist(ac, bins=NP.linspace(0, 100, 101), weights=NP.ones(len(ac)) / len(ac), cumulative=True)
    for ax in axs:
        ax.set_ylabel('SNVs', fontsize=8)
    axs[0].set_xlabel('Cells with expression', fontsize=8)
    axs[1].set_xlabel('Cells with alternative allele', fontsize=8)
    axs[2].set_xlabel('Cells with only alternative allele (%)', fontsize=8)
    axs[3].set_xlabel('Overall allele fraction', fontsize=8)
    axs[3].hist(AFs, bins=NP.linspace(0, 100, 101), weights=NP.ones(len(AFs)) / len(AFs), cumulative=True)
    for ax in axs:
        ax.set_xticks(NP.arange(0, 110, 10))
        ax.set_xticks(NP.arange(0, 101, 1), minor=True)
    axs[2].hist(bc, bins=NP.linspace(0, 101, 501), weights=NP.ones(len(bc)) / len(bc), cumulative=True)
    if out is not None:
        fig.savefig(out, bbox_inches='tight')
    if close_fig:
        plt.close(fig)

    return [NP.array(cell_counts), NP.array(alt_counts), bc, NP.array(AFs)]

def snv_filtering(ref, alt, snvs, qa_data, min_cells = None, min_alts = None, min_AF = None, max_homozygous = None, apply=False):
    idx = NP.full(len(snvs), True, dtype='bool')
    print(f'TOTAL SNVs before filtering {len(idx)}')
    if min_cells is not None:
        idx &= (qa_data[0] >= min_cells)
    if min_alts is not None:
        idx &= (qa_data[1] >= min_alts)
    if max_homozygous is not None:
        idx &= (qa_data[2] < max_homozygous)
    if min_AF is not None:
        #print((qa_data[3] > min_AF).sum())
        idx &= (qa_data[3] >= min_AF)
    sidx = NP.where(idx)[0]
    return sidx
