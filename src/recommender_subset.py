#!/usr/bin/env python

'''
Recommender subset generator

Input:   We start with a four-column BED file like `removeOverlap.bed`, which 
         contains score data for each 1kb bin. We can optionally specify other
         parameters for changing `k` and window characteristics.

Output:  The output contains `k` or fewer intervals that are high-scoring and 
         do not overlap one another.
'''

import collections
import numpy as np
import pandas as pd
import pyranges as pr
import bisect
import click

def build_cumulative_sizes(ranges):
    '''
    Generates key-value pairs, where keys are chromosome names,
    and sizes are the cumulative size as one goes through a list
    of ordered chromosomes. 
    
    The order of chromosomes should reflect their ordering in pyranges-
    formatted intervals, as this table will be used to put intervals 
    onto one continuous line segment, ordered by their end coordinates.

    Input:

    ranges (PyRanges object) - 25kb intervals

    Output:

    cs (dict) - dictionary of chromosome names and cumulative sizes
    '''
    merged_ranges = ranges.merge()
    chroms = merged_ranges.as_df().loc[:, 'Chromosome'].copy().values
    starts = merged_ranges.as_df().loc[:, 'Start'].copy().values
    ends = merged_ranges.as_df().loc[:, 'End'].copy().values
    diffs = ends - starts
    diffs = np.concatenate([[0], diffs[:-1]])
    j = len(diffs) - 1
    while j > 0:
        diffs[j] = np.sum(diffs[0:j+1])
        j -= 1
    cs = dict(zip(chroms, diffs))
    return cs

def translate_chrom_to_abs_coords(ranges, cs):
    '''
    Translate start and end coordinates from chromosomal coordinate system
    to absolute coordinate system. This lets us search for an optimal path 
    over the genome, without having to do tricky things to handle multiple 
    chromosomes.

    Input:

    ranges (PyRanges object) - 25kb intervals to linearize
    cs (dict) - dictionary of chromosome names and cumulative sizes

    Output:

    starts (numpy array) - absolute start coordinates
    ends (numpy array) - absolute end coordinates
    '''
    chroms = ranges.as_df().loc[:, 'Chromosome'].copy().values
    starts = ranges.as_df().loc[:, 'Start'].copy().values
    ends = ranges.as_df().loc[:, 'End'].copy().values
    translate_genomic_coords_to_abs = np.vectorize(lambda x, y: y + cs[x])
    starts = translate_genomic_coords_to_abs(chroms, starts)
    ends = translate_genomic_coords_to_abs(chroms, ends)
    return [starts, ends]

def build_disjoint_interval_indices(starts, ends, n):
    '''
    Build indices of nearest leftmost disjoint intervals, given absolute
    start and end coordinates.

    Input:

    starts (numpy array) - absolute start coordinates
    ends (numpy array) - absolute end coordinates
    n (int) - number of intervals

    Output:

    di (list) - list of nearest leftmost disjoint interval indices 
    '''
    di = []
    for j in range(n):
        idx = bisect.bisect_right(ends, starts[j]) - 1
        di.append(idx)
    return di

def build_optimum_scoring_indices(ranges, disjoint_intervals, n):
    '''
    Build an optimum score path through our input intervals, and return
    a list of qualifying indices to use to filter ranges.

    Input:

    ranges (PyRanges object) - 25kb intervals to linearize
    disjoint_intervals (list) - list of nearest disjoint intervals
    n (int) - number of intervals

    Output:

    qi (list) - indices of high-scoring, non-overlapping intervals 
    '''
    df = ranges.as_df()
    '''
    Set up initial opt(imum) table values
    '''
    opt = collections.defaultdict(int)
    opt[-1] = 0
    opt[0] = 0
    '''
    Forwards pass that builds best-score path
    '''
    for j in range(1, n):
        score = df.loc[df.index[j], 'Score']
        opt[j] = max(score + opt[disjoint_intervals[j]], opt[j - 1])
    '''
    Backwards trace to retrieve path indices
    '''
    qi = []
    j = n - 1
    while j >= 0:
        score = df.loc[df.index[j], 'Score']
        '''
        If high-scoring, add qualifying index and jump to the 
        next disjoint interval to the left, otherwise just try the 
        next interval to the left
        '''
        if score + opt[disjoint_intervals[j]] > opt[j - 1]:
            qi.append(j)
            j = disjoint_intervals[j]
        else:
            j -= 1
    '''
    Sort qualifying values so that we can later read our original 
    dataframe by index
    '''
    return sorted(qi)

def no_overlaps(df, k, extend, window_size):
    '''
    Build an optimum score path through our input intervals, and return
    a list of qualifying indices to use to filter ranges.

    Input:

    df (Pandas dataframe) - dataframe with `Chromosome`, `Start`, `End`, and 
      `Score` columns, minimally
    k (int) - number of desired intervals meeting score and overlap criteria
    extend (int) - number of bases to extend 1k intervals up- and downstream
    window_size (int) - width of input interval (1k, for instance)

    Output:

    locs (Pandas dataframe) - best-scoring, non-overlapping intervals
    '''
    gr = pr.PyRanges(df)
    '''
    Expand and filter all pyranges object intervals to those 25kb in width
    '''
    gr = gr.slack(extend)
    span_size = 2 * extend + window_size
    gr = gr.apply(lambda df: df[(df['End'] - df['Start']) == span_size])
    '''
    Build cumulative size dictionary for calculating absolute coordinates
    and then linearize the 25kb intervals to a fake chromosome
    '''
    cs = build_cumulative_sizes(gr)
    [abs_starts, abs_ends] = translate_chrom_to_abs_coords(gr, cs)
    is_sorted = lambda a: np.all(a[:-1] <= a[1:])
    assert(is_sorted(abs_ends)) # make sure end positions are ordered
    '''
    Get indices of nearest disjoint intervals along linearized pseudo-chromosome
    '''
    n = len(gr.as_df().index)
    di = build_disjoint_interval_indices(abs_starts, abs_ends, n)
    '''
    Trace optimum score path and get qualifying indices of elements in 25kb set
    '''
    qi = build_optimum_scoring_indices(gr, di, n)
    '''
    Make dataframe of qualifying intervals
    '''
    qf = gr.as_df().iloc[qi]
    '''
    Return dataframe of k, best-scoring intervals. We might get fewer intervals 
    than we'd like (although that's not the case here) so we take the minimum of 
    either k or the number of qualifying indices
    '''
    k = min(k, len(qi))
    locs = qf.sort_values(['Score'], ascending=False).iloc[:k].sort_values(['Chromosome', 'Start', 'End'])
    return locs

@click.command()
@click.option('--input', help='input filename for 1k intervals with scores')
@click.option('--k', type=int, default=100000, help='number of desired intervals')
@click.option('--extend', type=int, default=12000, help='number of bases used to extend input intervals')
@click.option('--window-size', type=int, default=1000, help='number of bases in each input interval')
def main(input, k, extend, window_size):
    '''
    Note
    ----
    I am building the dataframe `df` from the provided BED4 file 
    (`removeOverlaps.bed`). This may be a good place to construct a 
    dataframe with `Chromosome`, `Start`, `End`, and `Score` columns 
    from the original script's 'Chromosome', 'Start', 'Score', and 
    'IDX' columns.
    '''
    df = pd.read_csv(input, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score'])
    locs = no_overlaps(df, k, extend, window_size)
    print(locs)

if __name__ == '__main__':
    main()