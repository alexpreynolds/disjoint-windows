#!/usr/bin/env python

import sys

FIELD_DELIM = ';'

for line in sys.stdin:
    (chrom, start, stop, scores_str) = line.rstrip().split('\t')
    scores = [float(x) for x in scores_str.split(FIELD_DELIM)]
    max_score = max(scores)
    max_score_index = scores.index(max_score)
    score = str(max_score)
    sys.stdout.write('{}\n'.format('\t'.join([chrom, start, stop, score])))