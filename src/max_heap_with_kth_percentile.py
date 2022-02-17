#!/usr/bin/env python

'''
In this approach, we use a priority queue (min-heap) of scores and 
indices.

We read score values into the heap, along with an index key. This
key refers back to an interval.

We start by popping the current highest scoring element off 
the heap.

We look for overlaps of this candidate element with all other 
elements in our "qualifying list", using a "rejection vector". 

The rejection vector contains a boolean flag for each interval. 
We can look upstream and downstream of some index into this vector 
by N intervals, where N is defined by `--window-span`.

In our use case, we have 4.2kb intervals that step by 200nt increments
across the genome.

This means that for a given bin, we look at the rejection flags for 
10 intervals upstream and 10 intervals downstream, looking for bins 
that were already marked as "qualifying".

In other words, any values that are "True" indicate that the 
candidate bin will have an overlap with some bin upstream 
or downstream, which was previously marked as "qualifying". 

In case no such overlap is found (such as the very first pop), we 
add the candidate element to the qualifying list, then pop another 
element and repeat the test.

In case an overlap is found somewhere, and if it is smaller than the
kth-percentile value passed into this script, we skip this candidate 
and pop another element off the heap.

At the end, we print out all `k` elements that qualified.

Notes
-----

 - Deterministic, but need to investigate how stable this is, i.e., 
   in case of ties, which element gets popped? This could change which 
   interval gets picked first, with knock-on effects for calling overlaps 
   downstream. Python docs have some notes on sort stability.

 - This could discard qualifying intervals at 4.2kb edges of chromosomes.
   Could be fixed by using a chromsizes file to define edge indices, which
   in turn modify the bounds given to the np.any() test.

 - Are high scores concentrated in some chromosomes, or distributed evenly
   over the genome? Concentrated scores might bias the set of elements 
   that get popped off the heap, tested, and qualified.

'''

import sys
import heapq
import pandas as pd
import numpy as np
from io import StringIO
import click

@click.command()
@click.option('--input', help='input filename')
@click.option('--k', type=int, help='samples')
@click.option('--window-span', type=int, default=10, help='number of windows used for overlap/rejection testing')
@click.option('--kth-percentile', type=float, default=0, help='value of kth percentile weight (e.g., 95%)')
@click.option('--is-final-pass', default=False, is_flag=True, help='is this a final pass on data')
def main(input, k, window_span, kth_percentile, is_final_pass):
    column_names=['Chromosome', 'Start', 'End', 'Score']
    df = pd.read_csv(input, sep='\t', header=None, names=column_names)
    n = len(df.index)
    r = np.zeros(n, dtype=np.bool_)
    w = df.loc[:, 'Score'].values
    h = []
    '''
    We want to prioritize higher scores, but the built-in Python 
    priority queue uses a min-heap structure, so we need to take the 
    inverse of the score. Popping then gives us the highest scoring 
    element.
    '''
    for i, v in enumerate(w):
        heapq.heappush(h, (-v, i))
    '''
    Build a qualifying list, until we have as many elements as we want,
    or until we have emptied the heap of testable, non-zero-score indices.
    '''
    q = []
    while k > 0:
        try:
            (v, i) = heapq.heappop(h)
            if v == 0:
                break
            start_i = (i - window_span) if (i - window_span) > 0 else 0
            stop_i = (i + window_span + 1) if (i + window_span + 1) <= n else n
            # sys.stderr.write('k {} | testing key {} | bounds [{}:{}] -> {}'.format(k, i, start_i, stop_i, np.any(r[start_i:stop_i])))
            '''
            If the test weight is gte to the kth-percentile, we keep it
            '''
            if not np.any(r[start_i : stop_i]) or (-v >= kth_percentile and not is_final_pass):
                r[i] = True
                q.append(i)
                k -= 1
        except IndexError:
            k = 0

    '''
    Sort indices and print.
    '''
    q.sort()
    o = StringIO()
    df.iloc[q].to_csv(o, sep='\t', index=False, header=False)
    sys.stdout.write('{}'.format(o.getvalue()))

if __name__ == '__main__':
  main()