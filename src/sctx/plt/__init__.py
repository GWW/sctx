from .formatax import format_subplots, format_axis, add_ax
import matplotlib
import numpy as NP
import matplotlib
import pylab as plt
from .formatax import format_subplots as fs


def rand_jitter(N, mid, diff = 0.4, stdev = 0.120):
    rn = (mid - diff, mid + diff)
    ret = NP.zeros(N, dtype='float') + NP.random.randn(N) * stdev + mid
    ret[ret > rn[1]] = rn[1]
    ret[ret < rn[0]] = rn[0]
    return ret

def get_colors():
    cb = plt.cm.tab20(NP.linspace(0, 1, 20))
    c1 = cb[NP.arange(0, 20, 2)]
    c2 = cb[NP.arange(1, 20, 2)]
    c3 = plt.cm.Dark2(NP.linspace(0, 1, 8))
    vcolors = NP.concatenate((c3, c1, c2))
    return vcolors

def density_scatter(ax, x, y, s=10, cmap=plt.cm.plasma, logx=False, logy=False):
    lx, ly = x, y
    if logx:
        lx = NP.log10(x)
    if logy:
        ly = NP.log10(y)
    xy = NP.vstack([lx,ly])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    ret = ax.scatter(x, y, c=z, s=s, edgecolor='', cmap=cmap, marker='o')
    #locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=NP.arange(2, 10) * .1,
    #                                      numticks=100)
    if logx:
        ax.set_xscale('log')
        #ax.xaxis.set_major_locator(locmaj)
        #ax.xaxis.set_minor_locator(locmin)
        #ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if logy:
        ax.set_yscale('log')
        #ax.yaxis.set_major_locator(locmaj)
        #ax.yaxis.set_minor_locator(locmin)
        #ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
