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


from .snvs import calc
import sctx.filter as filt
from collections import Counter
import statsmodels.api as sm
from .plt import format_subplots, rand_jitter
from scipy.sparse import csr_matrix
import pylab as plt
from matplotlib.lines import Line2D
import numpy as NP
from matplotlib.colors import ListedColormap, Normalize
from scipy.stats import binom
from sklearn.neighbors import NearestNeighbors

#cmap = ListedColormap(NP.concatenate([[plt.cm.Dark2.colors[-1]], plt.cm.Dark2.colors[:-1]]))
cmap = plt.cm.Dark2
norm = Normalize(vmin=0, vmax=cmap.N - 1)

def iter_snvs(ref, alt):
    for i, (s, e) in enumerate(zip(ref.indptr, ref.indptr[1:])):
        yield i, ref.indices[s:e], ref.data[s:e], alt.data[s:e]

def bwhisks(vals):
        Q1, med, Q3 = NP.percentile(vals, [25, 50, 75])
        IQR = Q3 - Q1
        loval = Q1 - 1.5 * IQR
        hival = Q3 + 1.5 * IQR
        loval = NP.min(vals[vals >= loval])
        hival = NP.max(vals[vals <= hival])
        return loval, hival

class RefineGT(object):
    def __init__(self, ref, alt, snvs, barcodes, close_fig=True):
        self._ref = ref
        self._alt = alt
        self._snvs = snvs
        self._cells = ref.shape[1]
        self._close_fig = close_fig
        self._barcodes = barcodes
        self._labels = {}
        self._cmap = cmap
        self._norm = norm

    def set_params(self, **kwargs):
        for k, v in kwargs:
            if k not in self._params:
                print(f'Ignored invalid parameter {k} -> {v}')
            else:
                self._params[k] = v

    @property
    def cells(self):
        return self._cells

    @property
    def snvs(self):
        return len(self._snvs)

    def _calculate_scores(self, labels):
        ref = self._ref
        alt = self._alt
        ogts = calc.calculate_genotypes(ref, alt, labels.astype('int32'))
        #print(Counter(labels))
        x, y = ogts['cells_0'], ogts['cells_1']
        p = self._params
        idx = NP.where((x >= p['snv_min_cells_per_genotype']) & (y >= p['snv_min_cells_per_genotype']) & ((x + y) >= p['snv_total_cells']))[0]
        keep = idx[(NP.abs(ogts['mean_AF_0'][idx] - ogts['mean_AF_1'][idx]) >= p['snv_AF_diff'])]

        return self._scores(ogts, ref, alt, keep, labels=labels)

    def _scores(self, gts, ref, alt, keep, labels=None, NL = 2):
        p = self._params
        rr = ref[keep].tocsc()
        aa = alt[keep].tocsc()

        scores = NP.zeros((rr.shape[1], NL), dtype='double')
        counts = NP.zeros((rr.shape[1], NL), dtype='int')
        means = [gts[f'mean_AF_{l}'][keep] for l in range(NL)]
        eps = p['eps']
        eps1 = NP.log(eps)
        eps2 = NP.log(1 - eps)
        min_snvs = p['cell_min_snvs']

        for i, snvs, r, a in iter_snvs(rr, aa):
            k = NP.isfinite(means[0][snvs]) & NP.isfinite(means[1][snvs])
            ks = k.sum()
            for l in range(NL):
                m = means[l][snvs]
                pr = binom.logpmf(a[k], (a + r)[k], m[k])
                NP.clip(pr, eps1, eps2, out=pr)

                #scores[i, l] = (-pr).mean() if len(pr) >= min_snvs else NP.nan
                scores[i, l] = (-pr).mean() if len(pr) >= min_snvs else NP.nan
                counts[i, l] = len(pr)


        fchk = NP.all(NP.isfinite(scores), axis=1)
        finite = NP.where(fchk)[0]
        pf = scores[finite].max(axis=1) >= self._params['cell_min_score']
        passed = finite[pf]

        if labels is not None:
            npassed = finite[~pf]
            nlabels = labels.copy()
            nlabels[~fchk] = 3
            f = passed[labels[passed] == 3]
            labels[f] = NP.argmin(scores[f], axis=1)
            f = npassed[labels[npassed] != 3]
            labels[f] = 3
            return nlabels, scores, passed, keep, gts
        else:
            return scores, passed, keep, gts



    def refine_genotypes(self, rlabels,
            snv_min_cells_per_genotype=2, snv_total_cells=5, snv_AF_diff = 0.3, eps=1e-5,
            cell_min_snvs=10, cell_min_score=0.25, max_iterations=5, singlet_percentile=95,
            doublet_rate=0.10, doublet_simulate_multiplier=4, doublet_minimum=100, doublet_min_pval = 0.2,
            seed=42, out=None):

        #cell_initial_score_diff = 0.1, 
        self._params = {
            'snv_min_cells_per_genotype':snv_min_cells_per_genotype,
            'snv_total_cells':snv_total_cells,
            'snv_AF_diff':snv_AF_diff,
            'eps':eps,
            'cell_min_snvs':cell_min_snvs,
            'cell_min_score':cell_min_score,
            'max_iterations':max_iterations,
            'singlet_percentile':singlet_percentile,
            'doublet_simulate_multiplier':doublet_simulate_multiplier,
            'doublet_minimum':doublet_minimum,
            'doublet_rate':doublet_rate,
            'doublet_min_pval':doublet_min_pval,
            'seed':42
        }


        NP.random.seed(seed=seed)
        score_steps = []
        iteration = 1
        raw_ref = self._ref
        raw_alt = self._alt
        labels = rlabels.copy()
        score_steps = [None]
        dscores = None

        #print(p)
        print('initial: ', Counter(labels))

        scats = []

        _, scores, passed, passed_snvs, gts = self._calculate_scores(labels)
        labels = NP.full(len(scores), 3, dtype='int')
        labels[passed] = NP.argmin(scores[passed], axis=1)

        cx = Counter(labels)
        if cx[1] < 5:
            print('Found fewer than 5 cells with a second genotype giving up')
            return False

        cutoffs = []
        for l in range(2):
            lidx = NP.where(labels == l)[0]
            sl = scores[lidx,l]
            co = NP.percentile(sl, self._params['singlet_percentile'])
            cutoffs.append(co)
        p = (scores[passed, 0] >= cutoffs[0]) & (scores[passed, 1] >= cutoffs[1])
        labels[passed[p]] = 2


        label_steps = [labels]
        #fig, axs = format_subplots(max_iterations, 3, bsize=3)
        while iteration <= max_iterations:
            labels, scores, passed, passed_snvs, gts = self._calculate_scores(label_steps[-1])
            cc = Counter(labels)
            if cc[1] <= 20:
                nlabels = NP.full(len(scores), 3, dtype='int')
                nlabels[passed] = NP.argmin(scores[passed], axis=1)
                #ax = axs[iteration - 1]
                ss = scores[passed]
                #ax[0].scatter(ss[:,0], ss[:,1], c=labels[passed], s=4, lw=0, cmap=cmap, norm=norm, marker='o')
                #ax[2].scatter(ss[:,0], ss[:,1], c=nlabels[passed], s=4, lw=0, cmap=cmap, norm=norm, marker='o')
            else:
                nlabels, dscores = self._find_doublets([scores, passed, passed_snvs, gts], labels, None) #axs[iteration - 1])

            if passed.sum() < 10:
                print('Bailing out the second genotype was likely not found during the rescue')
                return False


            #axs[iteration - 1, 1].scatter(dcscores[dpassed,0], dcscores[dpassed, 1], color=plt.cm.Dark2.colors[3], s=3, cmap=cmap, lw= 0, marker='o', norm=norm) 

            #print(Counter(nlabels))

            ccounter = Counter(zip(label_steps[-1], nlabels))
            label_steps.append(nlabels)
            score_steps.append(scores)
            changed = {k:v for k, v in ccounter.items() if k[1] >= 0 and k[0] != k[1]}
            cs = sum(changed.values())
            print(f'Iteration {iteration}: number of changed labels: ',cs, Counter(nlabels))
            ks = changed.get((3, 0), 0) + changed.get((3, 1), 0) + changed.get((3, 2), 0)
            iteration += 1
            if (cs == 0) or (cs <= 3 and len(changed) == 1 and ks <= 3):
                break

        #for a in axs[iteration - 1:]:
        #    for ax in a:
        #        ax.remove()

        self._label_steps = label_steps
        self._final_labels = label_steps[-1]
        self._final_scores = score_steps[-1]

        fig, axs = format_subplots(2, 1, bsize=4)
        ax = axs[1]
        ss = score_steps[-1]
        labels = label_steps[-1]
        lnames = ['Genotype 1', 'Genotype 2', 'Doublet']
        for l in (0, 1, 2):
            passed = labels == l
            ax.scatter(ss[passed,0], ss[passed,1], color=plt.cm.Dark2.colors[l], marker='o', label=lnames[l])
        ax.set_ylabel('GT2 score')
        ax.set_xlabel('GT1 score')
        ax.legend(loc='upper left', bbox_to_anchor=[1.05, 1])

        ax = axs[0]
        ax.set_ylabel('GT2 score')
        ax.set_xlabel('GT1 score')
        passed = labels < 3
        ax.scatter(ss[passed, 0], ss[passed, 1], color='0.5', marker='o', label='Real Cells')
        ax.scatter(dscores[:, 0], dscores[:, 1], color=plt.cm.Dark2.colors[4], marker='o', label='Simulated Doublets')
        ax.legend(loc='upper left', bbox_to_anchor=[1.05, 1])


        if out is not None:
            fig.savefig(out, bbox_inches='tight')
        if self._close_fig:
            plt.close(fig)


        return True

    def _generate_doublets(self, sdata, labels):
        def get_snvs(ref, alt, cidx):
            s, e = ref.indptr[cidx:cidx+2]
            return ref.indices[s:e], ref.data[s:e], alt.data[s:e]

        l1 = NP.where(labels == 0)[0]
        l2 = NP.where(labels == 1)[0]

        p = self._params
        D = p['doublet_rate']
        D = int(min(len(l1), len(l2)) * D * p['doublet_simulate_multiplier'])
        D = max(D, p['doublet_minimum'])

        c1 = NP.random.choice(l1, D)
        c2 = NP.random.choice(l2, D)

        ref = self._ref.tocsc()
        alt = self._alt.tocsc()

        dscores = NP.zeros((D, 2), dtype='double')
        dcounts = NP.zeros((D, 2), dtype='int')

        indices = []
        indptr = []
        rdata = []
        adata = []

        snv_tots = NP.zeros(alt.shape[0], dtype='int')
        for i, (i1, i2) in enumerate(zip(c1, c2)):
            s1, r1, a1 = get_snvs(ref, alt, i1)
            s2, r2, a2 = get_snvs(ref, alt, i2)   

            snvs = NP.union1d(s1, s2)
            r = NP.zeros(len(snvs), a1.dtype)    
            a = NP.zeros(len(snvs), a1.dtype)

            _, s1i, ssi = NP.intersect1d(s1, snvs, assume_unique=True, return_indices=True)
            r[ssi] += r1
            a[ssi] += a1    

            _, s2i, ssi = NP.intersect1d(s2, snvs, assume_unique=True, return_indices=True)
            r[ssi] += r2
            a[ssi] += a2

            indptr.append(len(indices))
            indices.extend(snvs)
            rdata.extend(r)
            adata.extend(a)

        indptr.append(len(indices))
        dref = csr_matrix((rdata, indices, indptr), shape=(D, ref.shape[0])).T.tocsr()
        dalt = csr_matrix((adata, indices, indptr), shape=(D, ref.shape[0])).T.tocsr()
        scores, passed, passed_snvs, gts = sdata
        dscores, dpassed, dkeep, dgts = self._scores(gts, dref, dalt, passed_snvs)

        return dscores[dpassed]

    def _find_doublets(self, sdata, labels, ax):
        scores, passed, passed_snvs, gts = sdata
        #labs = labels[passed]
        labs = NP.full(len(scores), 3, dtype='int')
        labs[passed] = NP.argmin(scores[passed], axis=1)

        ds = self._generate_doublets(sdata, labs)
        dlabs = NP.argmin(ds, axis=1)
        from sklearn.neighbors import KNeighborsClassifier

        nlabels = labels.copy()
        dprob = labels.copy()

        for l in range(2):
            lidx = NP.where(labs == l)[0]
            didx = NP.where(dlabs == l)[0]
            sl = scores[lidx]

            if l == 0:
                cutoffs = ds[didx,0].min(), ds[didx,1].max()
                pp = (sl[:,0] >= cutoffs[0]) & (sl[:,1] <= cutoffs[1])
                l1 = lidx[pp]
                l2 = lidx[~pp]
            else:
                cutoffs = ds[didx,0].max(), ds[didx,1].min()
                pp = (sl[:,0] <= cutoffs[0]) & (sl[:,1] >= cutoffs[1])
                l1 = lidx[pp]
                l2 = lidx[~pp]

            nlabels[l2] = l
            if len(l1) > 0:
                training = NP.concatenate([scores[labels == l], ds[didx]])
                tlabs = NP.concatenate([NP.repeat(l, (labels == l).sum()), NP.repeat(2, len(didx))])
                neigh = KNeighborsClassifier(n_neighbors=5)
                neigh.fit(training, tlabs)
                nlabels[l1] = neigh.predict(scores[l1])
        ss = scores[passed]
        if ax is not None:
            ax[0].scatter(ss[:,0], ss[:,1], c=labels[passed], s=4, lw=0, cmap=cmap, norm=norm, marker='o')
            ax[1].scatter(ss[:,0], ss[:,1], color='0.5', s=4, lw=0, cmap=cmap, norm=norm, marker='o')
            ax[1].scatter(ds[:,0], ds[:,1], color=plt.cm.Dark2.colors[4], s=4, lw=0, cmap=cmap, norm=norm, marker='o')
            ax[2].scatter(scores[passed,0], scores[passed,1], c=nlabels[passed], s=4, lw=0, cmap=cmap, norm=norm, marker='o')
        return nlabels, ds
