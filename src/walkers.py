#!/usr/bin/env python

'''
Walker's Alias method for weighted sampling with replacement, using
rejection sampling to retrieve desired number of disjoint elements

References:

 - https://en.wikipedia.org/wiki/Alias_method
 - http://www.nrbook.com/devroye/Devroye_files/chapter_three.pdf (pg. 107)

Motivation for using this approach:

 - O(n) preprocessing time
 - O(1) sampling time
 - rejection sampling preserves distribution

We have intervals that have a score associated with each of them.

We use this score as a weight and perform weighted sampling by way 
of Walker's Alias.

We use rejection sampling to progressively build a pool of non-
overlapping, high-scoring intervals, stopping when we have as many 
as desired (say, 100k).

As we build this pool, for each sample, we either accept or reject 
that sampled interval, depending on whether that sample overlaps 
any other samples already in our pool.

We have 3M 25k windows to consider. All windows overlap each other 
by 1kb and span the genome.

Each window has a score, so we have a 3M-long vector of weights
and another 3M-long vector of keys (line indices).

We keep a 3M-element vector of booleans, intialized to False values. 
Each element is treated as a line index, which refers back to an 
interval in the original interval file. 

As we sample, we mark the sampled key/row-index as True, but only if 
there are no other True elements in the flag vector within 23 indices
(in both directions, so 46 tests). 

Any such True values indicate an overlap with another interval 
upstream or downstream.

Because this is a numpy boolean vector, we should be able to slice 
this vector down to a view, and then quickly test for any True values 
within this view. To further speed this up, any True value that is 
found will trigger an early rejection.

If we don't find any True values nearby, we keep the interval. 

If we do, we throw that interval away and keep looking/sampling.

Once we have the desired number of samples, we're done. We write out
the sample to a dataframe or file for further inspection.

Issues
------

Non-deterministic, unless you use the same random seeds. You'll get a
different answer on each run. In principle, you should get good, high-
scoring intervals.

One complication is chromosomal borders, which we ignore for now. The 
ends of chromosomes probably have low weights, so even if we reject a 
sample that (falsely) is called an overlap with an element on the start 
or end of another chromosome, it was probably unlikely we were going 
to sample many (or any) of these in the first place.

As we are sampling with replacement, we might get unlucky and pick
the same interval twice. However, if the distribution of weights is 
reasonable, and given that we are starting with 3M elements, we
expect that this is unlikely to occur too often.

We reject elements which overlap. As overlaps get wider, the odd are
that we will overlap with other elements upstream or downstream and 
trigger a rejection. We may end up spinning our wheels, continually
sampling high-score elements which are simply too close to other
high-score elements.

'''

import pandas as pd
import numpy as np
import random
import click

def construct_walkers_alias_table(keys, weights):
  global weighted_probs
  global aliases
  num_weights = len(weights)
  sum_weights = np.sum(weights)
  weighted_probs_f = lambda w: w * (num_weights / sum_weights) # average weight of 1
  weighted_probs_f_vec = np.vectorize(weighted_probs_f)
  weighted_probs = weighted_probs_f_vec(weights)
  #print(weighted_probs)
  aliases = np.full(num_weights, -1)
  short_items = np.where(weighted_probs < 1)[0]
  long_items = np.where(weighted_probs > 1)[0]
  #print(short_items)
  #print(long_items)
  while short_items.size > 0 and long_items.size > 0:
    j, short_items = short_items[-1], short_items[:-1]
    k = long_items[-1]
    aliases[j] = k
    weighted_probs[k] -= (1 - weighted_probs[j])
    if weighted_probs[k] < 1:
      np.append(short_items, k)
      long_items = long_items[:-1]
  #print(weighted_probs)
  #print(aliases)

'''
The first uniform random variable selects a component in 
the equiprobable list, while the second decides which 
segment of the two-point distribution is selected. A lookup
table returns the sampled key.
'''
def sample_key(num_keys):
  u = random.uniform(0, 1)
  j = random.randint(0, num_keys - 1)
  idx = j if u <= weighted_probs[j] else aliases[j]
  return keys[idx]

@click.command()
@click.option('--input', help='input filename')
@click.option('--k', type=int, help='samples')
@click.option('--window-span', type=int, default=23, help='number of windows used for overlap/rejection testing')
@click.option('--max-attempts', type=int, default=100000, help='number of rejections before quitting early')
def main(input, k, window_span, max_attempts):
  global keys
  global num_keys
  df = pd.read_csv(input, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score'])
  # Prepare metadata
  num_keys = len(df.index)
  keys = list(range(0, num_keys))
  # Prepare rejection vector
  rejection_v = np.zeros(num_keys, dtype=np.bool)
  # Scrape score column into weights
  w = df.loc[:, 'Score'].values
  # Preprocess weights into alias table
  construct_walkers_alias_table(keys, w)
  # Generate samples
  num_samples = k
  qualifying_keys = []
  # Keep sampling until we get as many as we need
  while num_samples > 0:
    # Apply rejection criteria
    while True:
      if max_attempts == 0:
        #print(np.argwhere(rejection_v == False))
        print(len(qualifying_keys))
        raise ValueError("Ran into maximum sample attempts; reduce --k parameter")
      key = sample_key(num_keys)
      '''
      Add the sample to the qualifying list, if it is not 
      within 47 1k intervals (-23L/1C/+23R). 
      
      To do this, we first set up start and stop indices, 
      checking boundaries between [0, num_keys-1]. 
      
      Over this index range, we check values of |rejection_v|.
      
      If any values are True, we |continue| the loop. 
      
      If all are False, we set these same values all True, add
      add the key to |qualifying_keys|, decrement |num_samples|
      and |break| to look for the next qualifying sample.
      '''
      if key not in qualifying_keys:
        start_index = key - window_span if key - window_span > 0 else 0
        stop_index = key + window_span if key + window_span <= num_keys else num_keys
        # print('testing key {} bounds [{}:{}] -> {}'.format(key, start_index, stop_index, np.any(rejection_v[start_index:stop_index])))
        if np.any(rejection_v[start_index:stop_index]):
          max_attempts -= 1
          continue
        else:
          rejection_v[start_index:stop_index] = True
          qualifying_keys.append(key)
          num_samples -= 1
          break
  # Sort the keys into numerical order to retrieve dataframe rows by index
  qualifying_keys.sort()
  # Write out sample for inspection (could instead drop this to a Pandas dataframe)
  print(df.iloc[qualifying_keys])
  

if __name__ == '__main__':
  main()
