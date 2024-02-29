#Copyright(c) 2018-2020 Gavin W. Wilson

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
 
from .snvs import calc
import sctx.filter as filt
from .plt import format_subplots
import pylab as plt
import numpy as NP
import igraph, operator
from scipy.stats import median_abs_deviation
from .misc import knn_to_igraph
from collections import Counter
import pandas as pd
import leidenalg as lae
from scipy.stats import binom
from .refine import iter_snvs
from cyvcf2 import VCF


def iter_snvs(ref, alt):
    for i, (s, e) in enumerate(zip(ref.indptr, ref.indptr[1:])):
        yield i, ref.indices[s:e], ref.data[s:e], alt.data[s:e]

def bcutoff(data, iqrs=1.5, up=True):
    p = NP.percentile(data, [25, 75])
    iqr = p[1] - p[0]
    if up:
        co = data[data < (p[1] + iqr * 1.5)].max()
    else:
        co = data[data > (p[0] - iqr * 1.5)].min()
    return co

class TxDemulti(object):
    def __init__(self, ref, alt, snvs, barcodes, close_fig=False):
        self._ref = ref
        self._alt = alt
        self._snvs = snvs
        self._close_fig = close_fig
        self._cells = ref.shape[1]
        self._capplied = False
        self._sapplied = False
        self._divs = {}
        self._barcodes = barcodes
        self._labels = {}
        self._data = {}
        self._params = {}
        self._counts = {}

    @property
    def cells(self):
        return self._cells

    @property
    def snvs(self):
        return len(self._snvs)

    def calculate_initial_distances(self, threads=1):
        self._capplied = False
        self._ODA, self._ODC =  calc.calculate_AF_dist(self._ref, self._alt, threads=threads)

    def snv_qa_plot(self, out=None):
        if self._sapplied: 
            print('SNV cutoff already applied must make a new TxDemulti object if you want to adjust the cutoff')
            return
        self._snv_qa = filt.snv_qa_plot(self._ref, self._alt, out=out, close_fig=self._close_fig)

    def snv_filtering(self, min_cells = None, min_alts = None, min_AF = None, max_homozygous = None, apply=False):
        if self._sapplied: 
            print('SNV cutoff already applied must make a new TxDemulti object if you want to adjust the cutoff')
            return

        sidx = filt.snv_filtering(self._ref, self._alt, self._snvs, self._snv_qa, min_cells=min_cells, min_alts=min_alts, min_AF=min_AF, max_homozygous=max_homozygous)
        print(len(self._snvs), len(sidx), self._ref.shape)
        C = len(self._snvs) - len(sidx)
        P = 100.0 * C / len(self._snvs)
        self._params['snv_filtering'] = {'min_cells':min_cells, 'min_alts':min_alts, 'min_AF':min_AF, 'max_homozygous':max_homozygous}
        if apply:
            print(f'Cutoffs discarded {C:,}  / {len(self._snvs):,} SNVs [{P:.2f}%]')
            self._snvs = self._snvs.iloc[sidx]
            self._ref = self._ref[sidx].copy()
            self._alt = self._alt[sidx].copy()
            self._sapplied = True
        else:
            print(f'Cutoffs would discard {C:,}  / {len(self._snvs):,} SNVs [{P:.2f}%]')

    def set_genotypes(self, vcf, gt1, gt2, min_snvs=10, eps=1e-5):
        if self._sapplied: 
            print('SNV data already initialized must make a new TxDemulti object if you want to adjust the SNV data')
            return
        rows = []
        fo = VCF(vcf, samples=[gt1, gt2])
        s1 = fo.samples.index(gt1)
        s2 = fo.samples.index(gt2)
        gmap = {0:0.0, 1:0.5, 3:1.0, 2:None}
        for r in fo:
            if r.FILTER or r.is_indel or len(r.ALT) > 1:
                continue
            pos = r.POS - 1
            ref = r.REF 
            alt = r.ALT[0]
            chrom = r.CHROM

            g1 = r.gt_types[s1]
            g2 = r.gt_types[s2]

            if g1 == g2 or (g1 == 2 or g2 == 2):
                continue

            rows.append((chrom, pos, alt, gmap[g1], gmap[g2]))
        ngt = pd.DataFrame(rows, columns=['chrom', 'pos', 'alt', 'gt1', 'gt2'])
        keep = self._snvs.merge(ngt, on=['chrom', 'pos', 'alt'], how='inner')

        sidx = keep['snv_idx']
        self._snvs = keep
        self._ref = self._ref[sidx].copy()
        self._alt = self._alt[sidx].copy()
        self._sapplied = True

        m1 = keep['gt1'].values.copy()
        m2 = keep['gt2'].values.copy()
        means = [m1, m2]

        rr = self._ref.tocsc()
        aa = self._alt.tocsc()
        NL = 2
        scores = NP.zeros((rr.shape[1], NL), dtype='double')
        counts = NP.zeros((rr.shape[1], NL), dtype='int')
        eps1 = NP.log(eps)
        eps2 = NP.log(1 - eps)

        for i, snvs, r, a in iter_snvs(rr, aa):
            k = NP.isfinite(means[0][snvs]) & NP.isfinite(means[1][snvs])
            ks = k.sum()
            for l in range(NL):
                m = means[l][snvs]
                pr = binom.logpmf(a[k], (a + r)[k], m[k])
                NP.clip(pr, eps1, eps2, out=pr)
                scores[i, l] = (-pr).mean() if len(pr) >= min_snvs else NP.nan
                counts[i, l] = len(pr)

        fchk = NP.all(NP.isfinite(scores), axis=1)
        passed = NP.where(fchk)[0]
        labels = NP.full(len(scores), 3, dtype='int')
        labels[passed] = NP.argmin(scores[passed], axis=1)
        cc = Counter(labels)
        self._labels['initial_clustering'] = labels.copy()
        if cc[0] > 0 and cc[1] > 0:
            self._NL = 2
        else:
            self._NL = 1

        print(f'Kept {keep.shape[0]:,} SNVs with differing GT values from known genotyping data', cc)

    def make_qa_plot(self, ce_cutoffs = None, co_cutoffs = None, out=None):
        if ce_cutoffs is None:
            ce_cutoffs = NP.arange(0.0, 1.01, 0.01, dtype='double')
        if co_cutoffs is None:
            co_cutoffs = NP.arange(0, 1010, 10, dtype='uint32')    
        print("calculate 2d cutoffs")    
        h2d = 1 - calc.calculate_cutoffs_2d(self._ODC, co_cutoffs, ce_cutoffs)
        fig, ax = format_subplots(1, 1, figsize=(10,10))

        im = ax.imshow(h2d.T, aspect='equal', origin='lower', cmap=plt.cm.Spectral_r, interpolation='none', vmin=0, vmax=1.0) #, norm=LogNorm())\
        ax.set_xticks(NP.arange(0, h2d.shape[0])[::5])
        ax.set_xticks(NP.arange(0, h2d.shape[0]), minor=True)    
        ax.set_xticklabels([f'{x}' for x in co_cutoffs[::5]])
        ax.set_yticks(NP.arange(0, h2d.shape[1])[::5])
        ax.set_yticks(NP.arange(0, h2d.shape[1]), minor=True)    
        ax.set_yticklabels([f'{x:.2f}' for x in 100.0 * ce_cutoffs[::5]])

        cax = fig.add_axes([1.05, 0.5, 0.05, 0.25])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label('Cells Lost (%)', fontsize=12)
        ax.set_xlabel('Co-expressed SNVs Cutoff', fontsize=12)
        ax.set_ylabel('Percent of Cells Passing Co-expression Cutoff (%)', fontsize=12)
        if out is not None:
            fig.savefig(out, bbox_inches='tight')
        if self._close_fig:
            plt.close(fig)

        self._h2d = h2d
        self._cutoff_bins = [co_cutoffs, ce_cutoffs]

    def apply_cell_cutoffs(self, snv_cutoff, cell_cutoff):
        self._capplied = True
        cell_cutoff /= 100
        N = self._ref.shape[1]
        DC, DA = self._ODC.copy(), self._ODA.copy()
        self._DA, self._DC, self._filter_idx = calc.apply_cutoffs_to_matrices(DA, DC, snv_cutoff, cell_cutoff)
        passed = len(self._filter_idx)
        perc = 100.0 * passed / N
        self._params['apply_cell_cutoffs'] = {'snv_cutoff':snv_cutoff, 'cell_cutoff':cell_cutoff}
        print(f'Passed {passed:,} out of {N:,} [{perc:.2f}%]')

    def cluster(self, k=10, threads=2, resolution=0.4, seed=42, min_cells=5, networkx=False):
        if not self._capplied:
            print('Must apply a cutoff before initial clustering')
            return

        DA, DC = self._DA, self._DC
        knn_index, knn_dist = calc.calculate_AF_KNN(k, DA, DC, threads)
        G, X = knn_to_igraph(knn_index, knn_dist, networkx)

        part = lae.find_partition(G, lae.RBConfigurationVertexPartition, seed=seed,
                                      resolution_parameter = resolution, weights='weight')
        labels = NP.array(part.membership)
        labels = labels.astype('int32')

        cc = Counter(labels)
        if len(cc) > 2:
            print("Found more than 2 labels keeping the largest two clusters and marking the other cells as QA'd")
            print("Label counts:", end=" ")
            for k, v in sorted(cc.items()):
                print(f'({k} = {v:,})', end=' ')
            print("\n")
            clusters = sorted(cc.items(), key=operator.itemgetter(1), reverse=True)
            discard = [k[0] for k in clusters[2:]]
            labels[NP.isin(labels, discard)] = 3

        self._knn = (knn_index, knn_dist)

        #self._labels = labels
        olabels = NP.full(self._ref.shape[1], 3, dtype='int32')
        olabels[self._filter_idx] = labels
        olabels[olabels == -1] = 3

        print("Initial Label counts:", end=" ")
        CX = 0
        for k, v in sorted(Counter(olabels).items()):
            if 0 <= k < 3 and v >= min_cells:
                CX += 1
            print(f'({k} = {v:,})', end=' ')
        print("\n")

        if CX == 1:
            print('Only one genoptype found try rescue_genotype to find a second one')

        self._labels['initial_clustering'] = olabels.copy()
        self._params['cluster'] = {'k':k, 'resolution':resolution, 'seed':seed, 'min_cells':min_cells}
        self._G = G
        self._X = X
        self._NL = CX

    def rescue_genotype(self, out=None, eps=1e-5, upper_mads=10, lower_mads=3):
        if self._NL > 1:
            print('Two genotypes already identified')
            return
        params = {'eps':eps, 'upper_mads':upper_mads, 'lower_mads':lower_mads}
        self._genotypes_for_labels()
        ref = self._ref.tocsc() #[:,self._labels['initial_clustering'] == 0]
        alt = self._alt.tocsc() #[:,self._labels['initial_clustering'] == 0]
        scores = NP.zeros(ref.shape[1], dtype='double')
        counts = NP.zeros(ref.shape[1], dtype='int')
        means = self._ogts[f'mean_AF_0']
        eps1 = NP.log(eps)
        eps2 = NP.log(1 - eps)
        for i, snvs, r, a in iter_snvs(ref, alt):
            m = means[snvs]
            k = NP.isfinite(m)
            p = binom.logpmf(a[k], (a + r)[k], m[k])
            NP.clip(p, eps1, eps2, out=p)
            scores[i] = -p.mean()
            counts[i] = len(p)

        fig, axs = format_subplots(1, 1, bsize=4)
        med = NP.median(scores)
        mad = median_abs_deviation(scores)
        upper = med + upper_mads * mad
        lower = med + lower_mads * mad
        params['upper_cutoff'] = upper
        params['lower_cutoff'] = lower
        axs.axhline(upper, ls='--', color='r')
        axs.axhline(lower, ls='--', color='b')


        if out is not None:
            fig.savefig(out, bbox_inches='tight')
        if self._close_fig:
            plt.close(fig)

        labels = NP.full(len(scores), 3, dtype='int') #self._labels['initial_clustering'].copy()
        labels[(NP.isfinite(scores) & (scores >= upper))] = 1
        labels[(NP.isfinite(scores) & (scores <= lower))] = 0
        #self._labels = labels
        #self._rscore = scores
        self._data['rescue_score'] = scores
        self._labels['initial_clustering'] = labels.copy()
        self._params['rescue_genotype'] = params
        axs.set_ylabel('Genotype Divergence')
        axs.set_xlabel('Cell SNV counts')
        lnames = ['Primary Genotype', 'Secondary Genotype', 'Uncertain']
        for i, l in zip([0, 1, 3], lnames):
            lidx = labels == i
            axs.scatter(counts[lidx], scores[lidx], s=5, lw=0, edgecolor='k', color=plt.cm.Dark2.colors[i], rasterized=True, label=l)
        axs.set_xscale('log')
        axs.legend(loc='upper left', bbox_to_anchor=[1.05, 1])
        #self._rescue_labels = labels.copy()

        print("Rescued Label counts:", end=" ")
        CX = 0
        for k, v in sorted(Counter(labels).items()):
            if 0 <= k < 3:
                CX += 1
            print(f'({k} = {v:,})', end=' ')
        print("\n")
        if CX > 1:
            print('Rescued a second potential genotype')
            self._NL = CX
        else:
            print('Could not find a second potential genotype')


    def _genotypes_for_labels(self, key='initial_clustering'):
        self._ogts = calc.calculate_genotypes(self._ref, self._alt, self._labels['initial_clustering'].astype('int32'))

