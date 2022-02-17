#!/usr/bin/env python

import sys
import numpy as np

window_size = int(sys.argv[1])
bin_size = int(sys.argv[2])

EXEMPLAR_SIZE = 2 * window_size + bin_size
EXEMPLAR_HALF_SIZE = int(EXEMPLAR_SIZE / 2)
FIELD_DELIM = '|'

assert(EXEMPLAR_HALF_SIZE * 2 == EXEMPLAR_SIZE)
previous_line_range = None

def point_1d_centroid(data):
    x = data
    l = len(x)
    return sum(x) / l

def rand_jitter(arr):
    stdev = 0.00001 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev

for line in sys.stdin:
    (chrom, start, stop, scores_str) = line.rstrip().split('\t')
    start = int(start)
    stop = int(stop)
    if stop - start == EXEMPLAR_SIZE:
        sys.stdout.write('{}'.format(line))
        previous_line_range = (chrom, str(start), str(stop))
    elif stop - start > EXEMPLAR_SIZE:
        scores = [float(x) for x in scores_str.split(FIELD_DELIM)]
        # coerce to list if single-element
        if not isinstance(scores, list):
            scores = [float(scores)]
        scores = list(rand_jitter(scores))

        ordered_scores = sorted(scores)

        centroid = point_1d_centroid(ordered_scores)
        dist_scores = [abs(x - centroid) for x in ordered_scores]
        min_dist_scores = min(dist_scores)
        min_dist_scores_index = dist_scores.index(min_dist_scores)
        score_of_interest = ordered_scores[min_dist_scores_index]
        score_of_interest_index = scores.index(score_of_interest)
        new_start = start + EXEMPLAR_SIZE * score_of_interest_index
        new_stop = new_start + EXEMPLAR_SIZE

        if new_start == new_stop:
            new_stop += EXEMPLAR_SIZE

        if new_stop > stop:
            new_stop = stop
            new_start = new_stop - EXEMPLAR_SIZE

        new_line = '{}\n'.format('\t'.join([chrom, str(new_start), str(new_stop), str("{:.5f}".format(score_of_interest))]))
        new_line_range = (chrom, start, stop)
        if new_line_range != previous_line_range:
            sys.stdout.write(new_line)
        previous_line_range = (chrom, start, stop)

    else:
        print(line)
        raise ValueError("stop and start values are unexpected")