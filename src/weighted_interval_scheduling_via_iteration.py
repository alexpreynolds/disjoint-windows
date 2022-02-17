#!/usr/bin/env python

'''
In this approach, we use a weighted-interval scheduling algorithm to 
find an optimal packing of non-overlapping/disjoint/"mutually compatible"
elements, which are sorted lexicographically, i.e. by BEDOPS 'sort-bed' or
similar.

This uses dynamic programming with optimization from genomic intervals 
ordered by their stop position. Based on how we are creating our input, 
we can safely assume the stop position is sorted. If windows are of odd 
sizes, we may need to add sorting for a general use case.

Iteration is used to walk forwards through an ideal score path, and then 
backwards from the last interval. The backwards pass draws a list of high-
scoring intervals, with scores and windows optimally distributed over the 
input space.

This result contains non-overlapping elements. If the length of this result is 
greater than specified k, then we further constrain the result set using a 
max-heap or priority queue, which pulls off the best-scoring k elements from
this superset.

Notes
-----

Deterministic, although min-heap may have sort stability issues for ties. See 
max_heap.py for notes.

To handle multiple chromosomes, we convert input from a UCSC-like chromosomal
coordinate system to an "absolute" coordinate system. This lets us order 
coordinates linearly, from left to right, scanning for optimal weight solutions
over the whole genome.

We swap out recursion for an iterative approach, to avoid hitting Python stack 
depth limits.

'''

import pyranges as pr
import pandas as pd
import numpy as np
import bisect
import collections
import sys
import heapq
from io import StringIO
import click

@click.command()
@click.option('--input', help='input filename')
@click.option('--k', type=int, default=-1, help='samples')
def main(input, k):
    '''
    Use merged intervals to build absolute coordinate space later on.
    '''
    d = pr.read_bed(input)
    m = d.merge()
    df = d.as_df()
    df = df.rename(columns={"Name" : "Score"})
    n = len(df.index)
    '''
    Build accumulated sizes dictionary for calculating absolute coordinates.
    '''
    acc_chr = m.as_df().loc[:, 'Chromosome'].copy().values
    acc_abs = np.concatenate([[0], m.as_df().loc[:, 'End'].copy().values[:-1]])
    j = len(acc_abs) - 1
    while j > 0:
        acc_abs[j] = np.sum(acc_abs[0:j+1])
        j -= 1
    acc = dict(zip(acc_chr, acc_abs))
    '''
    For every interval j, compute the rightmost disjoint interval i, where i < j.
    '''
    s = df.loc[:, 'Start'].copy().values
    e = df.loc[:, 'End'].copy().values
    '''
    We first need to translate the start (s) and end (e) genomic coordinates to
    absolute coordinate space. This lets us search for an optimal path over the 
    whole genome, without having to do complicated things to handle multiple 
    chromosomes.
    '''
    c = df.loc[:, 'Chromosome'].copy().values
    translate_genomic_coords_to_abs = np.vectorize(lambda x, y: y + acc[x])
    s = translate_genomic_coords_to_abs(c, s)
    e = translate_genomic_coords_to_abs(c, e)
    
    '''
    Important! 
    ----------
    Intervals in our df must be sorted by their end coordinates, which we can safely
    assume here, as the input windows are built by stepping across the genome. 
    
    If we are working with other windows that are differently sized or of odd sizes, 
    sorting values by the end value will be required for the DP trace to work correctly.
    '''

    '''
    Next, we look for the nearest disjoint interval to the left of the current interval.
    '''
    p = []
    for j in range(n):
        i = bisect.bisect_right(e, s[j]) - 1
        p.append(i)

    '''
    Set up initial opt(imum) table values.
    '''
    opt = collections.defaultdict(int)
    opt[-1] = 0
    opt[0] = 0
    '''
    Forwards pass that puts the best-scoring path into opt.
    '''
    for j in range(1, n):
        score = df.loc[df.index[j], 'Score']
        opt[j] = max(score + opt[p[j]], opt[j - 1])
    '''
    Backwards trace to retrieve path.
    '''
    # given opt and p, find solution intervals in O(n)
    q = []
    j = n - 1
    while j >= 0:
        score = df.loc[df.index[j], 'Score']
        if score + opt[p[j]] > opt[j - 1]:
            q.append(j)
            j = p[j] # jump to the nearest disjoint interval to the left
        else:
            j -= 1 # try the "next" interval, one to the left
    '''
    Sort qualifying values so that we can read the dataframe by index.
    '''
    q.sort()
    '''
    Write elements to standard output, constraining to k elements, if possible.
    In either case, we are guaranteed non-overlapping elements.
    '''
    o = StringIO()
    if k == -1 or len(q) < k:
        sys.stderr.write('{} elements were found\n'.format(len(q)))
        df.iloc[q].to_csv(o, sep='\t', index=False, header=False)
        sys.stdout.write('{}'.format(o.getvalue()))
    else:
        sys.stderr.write('{} elements were found, and we pull {} best from this set\n'.format(len(q), k))
        # use heap to get k best-scoring element subset of larger superset
        h = []
        df_subset = df.iloc[q]
        w_subset = df_subset.loc[:, 'Score'].values
        '''
        While we want to prioritize higher scores, the built-in Python 
        priority queue uses a min-heap structure, so we need to take the 
        inverse of the score. Popping then gives us the highest scoring 
        element.
        '''
        for i, v in enumerate(w_subset):
            heapq.heappush(h, (-v, i))
        q_subset = []
        while k > 0:
            (v, i) = heapq.heappop(h)
            q_subset.append(i)
            k -= 1
        q_subset.sort()
        df_subset.iloc[q_subset].to_csv(o, sep='\t', index=False, header=False)
        sys.stdout.write('{}'.format(o.getvalue()))

if __name__ == '__main__':
  main()