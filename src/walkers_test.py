#!/usr/bin/env python

'''
Walker's Alias method for weighted sampling with replacement
https://en.wikipedia.org/wiki/Alias_method

O(n) preprocessing time
O(1) sampling time
'''

import sys
import random
import numpy as np

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

def sample_key(num_keys):
  u = random.uniform(0, 1)
  j = random.randint(0, num_keys - 1)
  random_int = j if u <= weighted_probs[j] else aliases[j]
  # lookup
  return keys[random_int]

if __name__ == '__main__':
  global keys
  global num_keys
  num_keys = 5
  keys = list(range(0, num_keys)) 
  #print(keys)
  w = np.array([0, 2, 4, 2, 2])
  # preprocess weights to probability table
  construct_walkers_alias_table(keys, w)
  # generate samples
  num_samples = int(sys.argv[1])
  for i in range(num_samples):
    print(sample_key(num_keys))
