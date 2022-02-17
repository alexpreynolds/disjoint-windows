#!/usr/bin/env python

import sys

distance = int(sys.argv[1])
curr = None
prev = None

def verify_neighbors(a, b):
  if not a or not b:
    return
  if a["chrom"] == b["chrom"]:
    try:
      assert(b["start"] - a["stop"] >= distance)
    except AssertionError:
      print("a -> {}:{}-{}".format(a["chrom"], a["start"], a["stop"]))
      print("b -> {}:{}-{}".format(b["chrom"], b["start"], b["stop"]))
      sys.exit(-1)

for line in sys.stdin:
  elems = line.rstrip().split('\t')
  curr = {
    "chrom" : elems[0],
    "start" : int(elems[1]),
    "stop" : int(elems[2])
  }
  if prev:
    verify_neighbors(prev, curr)
  prev = {
    "chrom" : elems[0],
    "start" : int(elems[1]),
    "stop" : int(elems[2])
  }