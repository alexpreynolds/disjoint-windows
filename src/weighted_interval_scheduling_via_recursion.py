#!/usr/bin/env python

'''
In this approach, we use a weighted-interval scheduling algorithm to 
find an optimal packing of non-overlapping/disjoint/"mutually compatible"
elements.

This uses dynamic programming with a O(n log n) sorting optimization on
genomic intervals ordered by their stop position. Recursion is used to 
"walk backwards" from an interval, looking for the best overall weight
distribution for a set of intervals.

Notes
-----

Recursive method for calculating solution causes Python stack limit to be
hit quickly on chromosome-sized input. Appears to work reasonably well on 
smaller inputs.

'''

import bisect
import pandas as pd
import collections
import sys
import threading
from io import StringIO
import click

@click.command()
@click.option('--input', help='input filename')
def main(input):
    df = pd.read_csv(input, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score'])
    n = len(df.index)
    '''
    For every interval j, compute the rightmost disjoint interval i, where i < j
    '''
    s = df.loc[:, 'Start'].values
    f = df.loc[:, 'End'].values
    p = []
    for j in range(n):
        i = bisect.bisect_right(f, s[j]) - 1
        p.append(i)
    '''
    Use DP to schedule weighted intervals. 
    '''
    opt = collections.defaultdict(int)
    opt[-1] = 0
    opt[0] = 0
    for j in range(1, n):
        score = df.loc[df.index[j], 'Score']
        opt[j] = max(score + opt[p[j]], opt[j - 1])
    # given opt and p, find solution intervals in O(n)
    q = []
    sys.setrecursionlimit(10**7)
    threading.stack_size(2**27)
    def compute_solution(j):
        if j >= 0:  # will halt on opt[-1]
            score = df.loc[df.index[j], 'Score']
            if score + opt[p[j]] > opt[j - 1]:
                q.append(j)
                compute_solution(p[j])
            else:
                compute_solution(j - 1)
    compute_solution(n - 1)
    q.sort()
    o = StringIO()
    df.iloc[q].to_csv(o, sep='\t', index=False, header=False)
    sys.stdout.write('{}'.format(o.getvalue()))

if __name__ == '__main__':
  main()