import pylab as plt
import matplotlib

def format_axis(ax, hide_xticks = False, hide_yticks = False):
    ax.get_yaxis().tick_left()
    ax.get_xaxis().tick_bottom()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.get_yaxis().set_tick_params(direction='out')
    ax.get_xaxis().set_tick_params(direction='out')
    if hide_xticks:
        for t in ax.xaxis.get_ticklines(): 
            t.set_visible(False)
    if hide_yticks:
        for t in ax.yaxis.get_ticklines(): 
            t.set_visible(False)
    return ax

def format_subplots(r, c, bsize = None, dpi=150, **kwargs):
    if bsize is not None:
        kwargs['figsize'] = (c * bsize, r * bsize)
    fig, axs = plt.subplots(r, c, dpi=dpi, **kwargs)
    if kwargs.get('squeeze', True) and r == 1 and c == 1:
        format_axis(axs)
    else:
        for ax in axs.flat:
            format_axis(ax)
    return fig, axs

def add_ax(fig, lft, top, width, height, **kwargs):
    ax = fig.add_axes([lft, 1 - top, width, height], **kwargs)
    return format_axis(ax)

fsp = format_subplots
fax = format_axis
