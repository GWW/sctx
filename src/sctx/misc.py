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

import igraph
import numpy as NP
import networkx as nx

def knn_to_igraph(knn, knnd, nx_graph = False):
    E = []
    W = []
    bad = 0
    for i, kk in enumerate(knn):
        for j, k in enumerate(kk):
            E.append((i, k))
            W.append(knnd[i, j])
    W = NP.array(W)
    g = igraph.Graph(n=knn.shape[0], edges=E, directed=True)
    g.es['weight'] = W
    x = None
    if nx_graph:
        x = nx.Graph()
        x.add_nodes_from(NP.arange(0, len(knn)))
        x.add_weighted_edges_from(((x[0], x[1], w) for x, w in zip(E, W)))

    return g, x

