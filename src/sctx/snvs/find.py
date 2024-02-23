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
import pandas as pd
from sklearn.metrics import roc_auc_score

def snv_calculate(cdata, snvs, alt, tot, min_expressed=10, gkey='genotype'):
    gts = cdata[gkey].values.copy()
    gt1 = NP.where(gts == 0)[0]
    gt2 = NP.where(gts == 1)[0]

    nz1 = NP.count_nonzero(tot[:,gt1], axis=1)
    nz2 = NP.count_nonzero(tot[:,gt2], axis=1)
    sidx = NP.where((nz1 >= min_expressed) & (nz2 >= min_expressed))[0]

    ridx = sidx[:,NP.newaxis]
    oaf1 = alt[ridx,gt1].sum(axis=1) / tot[ridx,gt1].sum(axis=1)
    oaf2 = alt[ridx,gt2].sum(axis=1) / tot[ridx,gt2].sum(axis=1)
    af1 = NP.nanmean(alt[ridx,gt1] / tot[ridx,gt1], axis=1)
    af2 = NP.nanmean(alt[ridx,gt2] / tot[ridx,gt2], axis=1)
   
    nz1 = nz1[sidx]
    nz2 = nz2[sidx]
    aucs = []
    gtp = NP.where((gts == 0) | (gts == 1))[0]
    gtu = gts[gtp]
    for i in sidx:
        n = NP.where(tot[i, gtp] > 0)[0]
        k = alt[i, gtp[n]] / tot[i, gtp[n]]
        aucs.append(roc_auc_score(gtu[n], k))
    t = NP.array(aucs) - 0.5

    dfo = snvs.loc[sidx, ['chrom', 'pos']].copy()
    dfo.reset_index(drop=True, inplace=True)
    dfo['sidx'] = sidx
    dfo['raw_auc'] = aucs
    dfo['auc'] = NP.abs(t) * 2
    dfo['auc_sign'] = NP.sign(t)       
    dfo['gt0_nz'] = nz1
    dfo['gt1_nz'] = nz2
    dfo['gt0_taf'] = oaf1
    dfo['gt1_taf'] = oaf2
    dfo['total_diff'] = oaf1 - oaf2
    dfo['gt0_maf'] = af1
    dfo['gt1_maf'] = af2
    dfo['mean_diff'] = af1 - af2
    return dfo

def find_set(cdata, sdf, alt, tot, min_covered = 1, min_diff=0.4, min_auc=0.75, gkey='genotype'):
    gts = cdata[gkey].values.copy()
    tidx = set()
    gtp = NP.where((gts == 0) | (gts == 1))[0]
    for g, s in zip((0, 1), (-1, 1)):
        covered = NP.zeros(len(gtp), dtype='int')
        keep = []
        used = set()
        idx = sdf[(sdf.auc >= min_auc) & (sdf.total_diff.abs() >= min_diff) & (sdf.auc_sign == s)]
        while True:
            mcov = 0
            mi = -1
            mz = None
            for i, ii in zip(idx.index, idx.sidx):
                if i in used:
                    continue
                nz = tot[ii, gtp] > 0
                cov = nz.sum()
                added = ((covered < min_covered) & nz).sum()                
                if added >= 5 and added > mcov:
                    mcov = added.sum()
                    mi = i
                    mz = nz
            if mi == -1:
                break
            
            used.add(mi)
            keep.append(mi)
            covered[mz] += 1
        N = len(tidx)
        tidx.update(keep)
        print(f'Genotype {g} Added = {len(tidx) - N} Total = {len(tidx)}')
    dfo = sdf.iloc[list(tidx)].copy()
    print(f'Covered Mean = {covered.mean():.2f} +/- {NP.std(covered):.2f} Zeros = {(covered == 0).sum()}')

    return dfo

def find_bad_cells(cdata, snvs, alt, tot, min_auc=0.75, min_diff=0.85, min_snvs=5, gkey='genotype', prefix=None):
    gts = cdata[gkey].values.copy()

    tset = snvs[(snvs.auc >= min_auc) & (snvs.total_diff.abs() >= min_diff)].copy()
    ttidx = tset.sidx.values.copy() 
    cidx = None
    taf = alt[ttidx] / tot[ttidx]

    cdevs = NP.zeros(taf.shape[1], dtype='double')
    cused = NP.zeros(taf.shape[1], dtype='int')
    odevs = NP.zeros(taf.shape[1], dtype='double')
    #print(taf.shape, len(gts))

    for gt in (0, 1, 2, 3):
        gidx = gts == gt
        if not NP.any(gidx):
            continue
        idx = NP.where(gidx)[0]
        if gt < 2:
            g = gt
            gn = 1 - gt
        else:
            g = 0
            gn = 1   
        g1 = tset[f'gt{g}_taf'].values
        #print(g1.min(), g1.max())
        g2 = tset[f'gt{gn}_taf'].values
        #print(min(idx), max(idx))
        for j in idx:
            d = taf[:,j]
            dev = -1
            odev = -1
            na = ~NP.isnan(d)
            used = na.sum()
            if used >= min_snvs:
                dev = NP.mean(NP.mean(NP.abs(g1[na] - d[na])))
                odev = NP.mean(NP.mean(NP.abs(g2[na] - d[na])))            

            if g < 2:
                cdevs[j] = dev
                odevs[j] = odev        
                cused[j] = used
            else:
                cdevs[j] = min(dev, odev)
                odevs[j] = max(dev, odev)
                cused[j] = used
    
    if prefix is None:
        prefix = ''
    else:
        prefix += '_'
    cdata[f'{prefix}match_deviance'] = cdevs
    cdata[f'{prefix}mismatch_deviance'] = odevs
    cdata[f'{prefix}snvs_used'] = cused
    return tset
