#!/usr/bin/env python

'''
This uses the method in max_heap.py. However, we add jitter to the score 
column before pop-and-test. This performs returns fewer elements before
exhausting the heap.
'''

import sys
import heapq
import pandas as pd
import numpy as np
import click


@click.command()
@click.option('--input', help='input filename')
@click.option('--k', type=int, help='samples')
@click.option('--window-span', type=int, default=23, help='number of windows used for overlap/rejection testing')
def main(input, k, window_span):
    df = pd.read_csv(input, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score'])
    n = len(df.index)
    r = np.zeros(n, dtype=np.bool)
    w = df.loc[:, 'Score'].copy().values
    j = np.random.normal(np.mean(w), np.std(w), w.size)
    w += j
    h = []
    '''
    While we want to prioritize higher scores, the built-in Python 
    priority queue uses a min-heap structure, so we need to take the 
    inverse of the score. Popping then gives us the highest scoring 
    element.
    '''
    for i, v in enumerate(w):
        heapq.heappush(h, (-v, i))
    '''
    Build a qualifying list, until we have as many elements as we want,
    or until we have emptied the heap.
    '''
    q = []
    while k > 0:
        try:
            (v, i) = heapq.heappop(h)
            start_i = i - window_span if i - window_span > 0 else 0
            stop_i = i + window_span if i + window_span <= n else n
            # sys.stderr.write('k {} | testing key {} | bounds [{}:{}] -> {}'.format(k, i, start_i, stop_i, np.any(r[start_i:stop_i])))
            if not np.any(r[start_i:stop_i]):
                r[start_i:stop_i] = True
                q.append(i)
                k -= 1
        except IndexError:
            k = 0

    '''
    Sort indices and print.
    '''
    q.sort()
    sys.stdout.write(df.iloc[q].to_string())

if __name__ == '__main__':
  main()