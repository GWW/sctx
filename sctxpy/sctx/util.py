import numpy as NP
import re
import numpy as NP
from operator import itemgetter


def boxtail(x, top=True, iqrs=1.5):
    Q1, median, Q3 = NP.percentile(x, [25, 50, 75])
    IQR = Q3 - Q1
    if top:
        return Q3 + iqrs * IQR
    return Q1 - iqrs * IQR

def sorted_nicely_key(l, key):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda item: [ convert(c) for c in re.split('([0-9]+)', key(item)) ]
    return sorted(l, key = alphanum_key)

def sorted_nicely_tuple(l):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda item: [ convert(c) for c in re.split('([0-9]+)', it) for it in item]
    return sorted(l, key = alphanum_key)

def sorted_nicely(l):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda item: [ convert(c) for c in re.split('([0-9]+)', item) ]
    return sorted(l, key = alphanum_key)


