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
from collections import Counter
from .plt import format_subplots, rand_jitter
from scipy.sparse import csr_matrix
import pylab as plt
import pandas as pd
from matplotlib.lines import Line2D
import numpy as NP
from sklearn.metrics import roc_curve, auc
import matplotlib
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.collections import PatchCollection
from scipy.stats import binom
from cyvcf2 import VCF
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def iter_snvs(ref, alt):
    for i, (s, e) in enumerate(zip(ref.indptr, ref.indptr[1:])):
        yield i, ref.indices[s:e], s, e, ref.data[s:e], alt.data[s:e]

class FindSNVs(object):
    def __init__(self, ref, alt, snvs, barcodes, labels, close_fig=True):
        self._ref = ref
        self._alt = alt
        self._snvs = snvs
        self._cells = ref.shape[1]
        self._close_fig = close_fig
        self._barcodes = barcodes
        self._labels = labels

    def score_snvs(self, snv_min_cells_per_genotype=2, snv_total_cells=5, snv_AF_diff = 0.3, top_snvs=5):
        self._params = {
            'snv_min_cells_per_genotype':snv_min_cells_per_genotype,
            'snv_total_cells':snv_total_cells,
            'snv_AF_diff':snv_AF_diff,
            #'eps':eps,
            'top_snvs':top_snvs,
            #'snv_score_diff':score_diff,
            #'snv_score_perc':score_perc
        }
        labels = self._labels
        ogts = calc.calculate_genotypes(self._ref, self._alt, labels.astype('int32'))
        x, y = ogts['cells_0'], ogts['cells_1']
        p = self._params
        idx = NP.where((x >= p['snv_min_cells_per_genotype']) & (y >= p['snv_min_cells_per_genotype']) & ((x + y) >= p['snv_total_cells']))[0]
        keep = idx[(NP.abs(ogts['mean_AF_0'][idx] - ogts['mean_AF_1'][idx]) >= p['snv_AF_diff'])]
        snv_scores, predictive, snv_index = self._scores(ogts, self._ref, self._alt, keep, labels=labels)
        self._keep = keep
        self._snv_scores = snv_scores
        self._predictive = predictive
        self._snv_index = snv_index

    def find_predictive(self, N = 5, max_snvs = 50, use_genotypes=False):
        totals = self._snv_scores['total_correct'].values
        gts = None
        if use_genotypes:
            gts = self._snv_scores['exome_gt1'].values
        ai = NP.argsort(totals)
        C = self._predictive.shape[1]
        cell_counts = NP.zeros(C, dtype='int')
        labs = self._labels
        glabs = NP.where(labs < 2)[0]
        pp = self._predictive
        snvs = []
        for i in ai[::-1]:
            if len(snvs) > max_snvs:
                break
            if use_genotypes and gts[i] == '':
                continue
            s, e = pp.indptr[i], pp.indptr[i + 1]
            ci = pp.indices[s:e]
            lidx = ci[labs[ci] < 2]
            updated = (cell_counts[lidx] < N).sum()
            if updated > 0:
                cell_counts[ci] += 1
                finished = (cell_counts[glabs] >= N).sum()
                snvs.append(i)
                if finished >= len(glabs):
                    break

        self._predictive_snvs = self._snv_scores.iloc[snvs]
        self._predictive_counts = cell_counts
        print(f'Needed {len(snvs)} SNVs to cover each singlet {N} times')

    def dotplot(self, label_map, out = None, close=False, use_genotypes=False):
        snvs = self._predictive_snvs
        idx = snvs['snv_idx'].values
        ref, alt = self._ref[idx], self._alt[idx]

        snames = NP.array([f'{r}: {p + 1}' for r, p in zip(snvs.chrom, snvs.pos)])

        good_index = [i for i, (b, l) in enumerate(zip(self._barcodes, self._labels)) if b in label_map and l < 3]
        print(f'Using {len(good_index)} out of {len(self._barcodes)} barcodes')
        labels_raw = [label_map[b] for b in self._barcodes[good_index]]
        label_names = NP.array(list(set(labels_raw)) + ['Doublets'])
        label_imap = {l:i for i,l in enumerate(label_names)}
        labels = NP.fromiter((label_imap[l] for l in labels_raw), dtype='int')
        ref, alt = ref[:,good_index], alt[:,good_index]

        NL = len(label_imap)
        gts = self._labels[good_index]

        labels[gts > 1] = NL - 1
        groups = sorted(Counter(zip(gts, labels)).items(), key=lambda k: (k[0][0], -k[1]), reverse=True)
        print(groups)
        gcounts = NP.fromiter((g[1] for g in groups), dtype='int')
        groups = NP.array([g[0] for g in groups])

        gmap = {tuple(d):i for i, d in enumerate(groups)}

        cafs = NP.zeros((len(groups), len(idx)), dtype='double')
        ctots = NP.zeros((len(groups), len(idx)), dtype='int')

        cell_rows = NP.fromiter((gmap[(g, l)] for g, l in zip(gts, labels)), dtype='int')

        for j, barcodes, start, end, r, a in iter_snvs(ref, alt):
            #print(cell_rows[barcodes])
            NP.add.at(ctots[:,j], cell_rows[barcodes], 1)
            AF = a / (r + a)
            NP.add.at(cafs[:,j], cell_rows[barcodes], AF)

        cperc = ctots / gcounts[:,NP.newaxis]
        R = (cperc * 0.90)  / 2
        A = NP.zeros(R.shape, dtype='double')
        A[ctots > 0] = cafs[ctots > 0] / ctots[ctots > 0]

        d1 = dist.pdist(A.T, metric='euclidean')
        rlink = sch.linkage(d1, method='ward') ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        cidx = sch.leaves_list(rlink) 

        A = A[:,cidx]
        R = R[:,cidx]

        fig, ax = format_subplots(1, 1, figsize=(8, 5), dpi=1000)
        ax.set_aspect('equal')

        # SNV's along the X axis 
        x, y = NP.meshgrid(NP.arange(R.shape[1]), NP.arange(R.shape[0]))
        circles = NP.array([plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)])
        norm = matplotlib.colors.Normalize(vmin=0.,vmax=1.)
        col = PatchCollection(circles, array=A.flatten(), cmap=plt.cm.Spectral_r, norm=norm, linewidth=0.5, edgecolor='k')
        ax.add_collection(col)
        ax.set_ylim(-0.5, R.shape[0] - 0.5)
        ax.set_xlim(-0.5, R.shape[1] - 0.5)
        ax.set_yticks(NP.arange(R.shape[0]))
        ax.set_xticks(NP.arange(R.shape[1]))
        ax.set_xticklabels(snames[cidx], fontsize=10, rotation=90)
        labs = [f'{l} N = {c}' for l, c in zip(label_names[groups[:,1]], gcounts)]
        ax.set_yticklabels(labs)
        ax.axhline(NP.where(groups[:,0] == 1)[0][0] - 0.5, color='k')
        ax.axhline(NP.where(groups[:,0] == 0)[0][0] - 0.5, color='k')
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)

        #for g, c in zip(groups, gcounts):
        #    print(g[0], label_names[g[1]], c)

        g1 = NP.where(groups[:,0] == 0)[0]
        g2 = NP.where(groups[:,0] == 1)[0]
        print(g1, g2)
        trans = ax.get_yaxis_transform() 
        ax.text(1.01, NP.median(g1), 'GT1', va='center', transform=trans, rotation=270, fontsize=15)
        ax.text(1.01, NP.median(g2), 'GT2', va='center', transform=trans, rotation=270, fontsize=15)

        pos = ax.get_position()

        W = pos.xmax - pos.xmin
        H = pos.ymax - pos.ymin
        WR = W / A.shape[1]
        HR = H / A.shape[0]
        cb_lgd = fig.add_axes([pos.xmax * 1.15, pos.ymax * 0.6, W / A.shape[1], H / A.shape[0] * 4 ])
        lcircles = [plt.Circle((0, i), radius=((r * 0.95) / 2), fill='k') for i, r in enumerate([0.25, 0.50, 0.75, 1])]
        lcol = PatchCollection(lcircles, array=NP.array([1] * len(lcircles)), cmap="binary", norm=norm, linewidth=0.5)
        cb_lgd.add_collection(lcol)
        cb_lgd.set_ylim(-0.5, 3.5)
        cb_lgd.set_xlim(-0.5, 0.5)
        cb_lgd.set_aspect('equal')
        cb_lgd.set_xticks([])
        cb_lgd.set_yticks([0, 1, 2, 3])
        cb_lgd.set_yticklabels([25, 50, 75, 100])
        cb_lgd.set_ylabel('Cells Expressing SNV (%)')

        pos = cb_lgd.get_position()
        cb_ax = fig.add_axes([pos.xmax * 1.1, pos.ymin, 0.05, 0.25 ])
        fig.colorbar(col, cax=cb_ax)
        cb_ax.yaxis.tick_left()
        cb_ax.set_ylabel('Mean Allele Fraction')
        cb_ax.yaxis.set_label_position("left")

        if use_genotypes:
            axd = make_axes_locatable(ax)
            gt_ax = axd.append_axes("top", size="20%", pad="1%")
            #gt_ax = fig.add_axes([pos.xmin, pos.ymax, W, 0.25 ])

            gt_labs = NP.zeros((2, len(snvs)), dtype='double')
            gt_labs[0] = [float(s) for s in snvs['exome_gt1'].values]
            gt_labs[1] = [float(s) for s in snvs['exome_gt2'].values]

            gt_ax.imshow(gt_labs[:,cidx], vmin=0, vmax=1, cmap=plt.cm.Spectral_r, aspect='auto')
            gt_ax.set_xticks([])
            gt_ax.set_yticks([0, 1])
            gt_ax.set_yticklabels(['GT1 Exome GT', 'GT2 Exome GT'], fontsize=8)


        if out is not None:
            fig.savefig(out, bbox_inches='tight')

        if close:
            plt.close(fig)

    def add_genotypes(self, vcf, gt1, gt2):
        fo = VCF(vcf, samples=[gt1, gt2])
        s1 = fo.samples.index(gt1)
        s2 = fo.samples.index(gt2)
        gmap = {0:0.0, 1:0.5, 3:1.0, 2:None}
        exome = {}
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
            exome[chrom, pos] = (ref, alt, gmap[g1], gmap[g2])
        dgt1 = []
        dgt2 = []
        for i, r in self._snv_scores.iterrows():
            t = (r.chrom, r.pos)
            d = exome.get(t)
            if d is not None and d[0] == r['ref'] and d[1] == r['alt']:
                dgt1.append(d[2])
                dgt2.append(d[3])
            else:
                dgt1.append('')
                dgt2.append('')
        self._snv_scores['exome_gt1'] = dgt1
        self._snv_scores['exome_gt2'] = dgt2

    def _scores(self, gts, ref, alt, keep, labels, NL = 2):
        p = self._params
        #rr = ref[keep].tocsc()
        #aa = alt[keep].tocsc()
        rr = ref[keep]
        aa = alt[keep]
        snvs = self._snvs.iloc[keep]
        N = p['top_snvs']
        means = [gts[f'mean_AF_{l}'][keep] for l in range(NL)]
        #eps = p['eps']
        #eps1 = NP.log(eps)
        #eps2 = NP.log(1 - eps)
        #smin = p['snv_score_diff']
        #pmin = p['snv_score_perc']
        #predictive = NP.full(len(rr.data), -1, dtype='int8')
        rows = []

        data = []
        indptr = [0]
        indices = []
        snv_index = []

        for i, barcodes, start, end, r, a in iter_snvs(rr, aa):
            if NP.isnan(means[0][i]) or NP.isnan(means[1][i]):
                continue
            snv = snvs.iloc[i]
            labs = labels[barcodes]
            lidx = NP.where(labs < 2)[0]
            ll = labels[barcodes[lidx]]

            AF = a[lidx] / (a + r)[lidx]
            m1 = means[0][i]
            m2 = means[1][i]

            plab = 1 if (m1 < m2) else 0
            fpr, tpr, thresholds = roc_curve(ll, AF, pos_label=plab) #, drop_intermediate=False)
            #a = NP.abs(auc(fpr, tpr) - 0.5) * 2
            a = auc(fpr, tpr)
            thresh = thresholds[NP.argmax(tpr - fpr)]

            g1c = (ll == plab) & (AF >= thresh)
            gt1c = g1c.sum()
            g2c = (ll != plab) & (AF < thresh)
            gt2c = g2c.sum()
            #predictive[lidx + start] = (g2c | g1c)

            kidx = NP.where(g2c | g1c)[0]
            indices.extend(barcodes[lidx[kidx]])
            #indices.extend(barcodes[lidx])
            #data.extend((g2c | g1c))
            data.extend(NP.ones(len(kidx), dtype='bool'))
            indptr.append(indptr[-1] + len(kidx))
            #indptr.append(indptr[-1] + len(lidx))
            snv_index.append(i)

            if plab == 1:
                gt1c, gt2c = gt2c, gt1c

            rows.append((keep[i], snvs.index[i][0], snvs.index[i][1], snv['ref'], snv['alt'], (ll == 0).sum(), (ll == 1).sum(), m1, m2, a, thresh, gt1c, gt2c, gt1c + gt2c))

        df = pd.DataFrame(rows, columns=('snv_idx', 'chrom', 'pos', 'ref', 'alt', 'gt1', 'gt2', 'mgt1', 'mgt2', 'auc', 'threshold', 'gt1_correct', 'gt2_correct', 'total_correct'))

        predictive = csr_matrix((data, indices, indptr), shape=(len(snv_index), ref.shape[1]))
        print(predictive.shape, len(predictive.data), len(predictive.indptr), len(snv_index))
        return df, predictive, NP.array(snv_index)
