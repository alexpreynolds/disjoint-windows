# disjoint-windows

We have the following goals:

1. Find a good, potentially optimal distribution of scored intervals over the genome, which should be disjoint, non-overlapping or, in other settings, called "mutually compatible".

2. We would like the selected method to prioritize high-scoring intervals.

3. We would like 100k such elements out of the full superset, which meet the first two criteria.

## Methods

Review targets in `src/makefile` for more concrete run instructions.

### `extend_windows`

We start with 1kb windows spanning the genome. Each window has a score.

The `extend_windows` target creates 25kb windows, which step across the genome in 1kb units. Each interval contains a score for that 25kb span, centered on the middle-1k window. There are ~3M such elements.

These windows are contained in `data/windows.fixed.25k.bed`. This file is used as input for other targets.

### `walkers`

In the `walkers` target, we perform weighted sampling with replacement using the [Walker's Alias](https://en.wikipedia.org/wiki/Alias_method) technique. We use the score of each input element as its weight, so samples are high-scoring. We reject samples which overlap the genomic space of elements already sampled.

We repeat sampling until we get 100k elements, or until we hit some iteration limit, as in this case it is unlikely we'll find any more qualifying elements.

### `priority_queue`

In the `priority_queue` target, we put elements into a priority queue, ordered by score. We pop the highest-scoring element off the queue, keep it if it does not overlap any other popped elements, and reject it if it does. 

We repeat this until the queue is empty, or until we get 100k elements. 

### `priority_queue_with_jitter` 

This target use the same method as `priority_queue` but adds noise to the weights, before constructing the heap. The idea is to get more separation between entries in the heap.

### `wis_via_iteration`

In the `wis_via_iteration` target, we use a dynamic programming approach to build a list of high-scoring elements which do not overlap one another.

The locally optimal path contains elements which do not overlap. If we have more than 100k such elements in this path, we use a priority queue to get the top-100k such elements.

## Results

### Baseline

Here is a summary of the baseline score distribution:

```
$ cut -f4 ../data/windows.fixed.25k.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05247 0.30294 0.36386 0.92198 0.73372 9.92041 
```

We expect that approaches we use should give better overall attributes, except for the maximum. We also expect that the approaches used return at least some elements that have this maximum score.

### `walkers`

Our Walker's Alias approach retrieves about 2.8% of input elements, before it retrieves the same high-scoring elements and must quit early. This amounts to ~86k of the 100k element subset we require, so we do not compare this method here.

### `priority_queue`

The max-heap approach returns 94972 of the requested 100k elements (95%). The heap gets exhausted faster than we can retrieve the elements we want. This may be due to high-scoring elements being close to one another, such that they are popped off the top of the heap and rejected too quickly. 

A `priority_queue_with_jitter` target adds noise to the weights, to try to get more separation before constructing the heap. This performed worse (92067 elements returned).

### `wis_via_iteration`

The weighted-interval scheduling approach of `wis_via_iteration` generates ~104k non-overlapping elements over our test input, and so we are able to return 100k high-scoring, non-overlapping intervals. We can use the `wis_via_iteration_summaries` target to compare the baseline signal distribution in the original windows against the distribution in our 104k- and 100k-window subsets.

We observe that the ~104k `wis_via_iteration` subset (no `k` specified) has a better signal profile, as we are picking a locally optimal set of high-score intervals:

```
$ cut -f4 ../results/output.wis_via_iteration.all.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.09265 0.38719 0.92402 2.17034 3.50988 9.92041 
```

As expected, the use of a priority queue to filter these 104k elements down to a 100k subset returns a better overall signal distribution:

```
$ cut -f4 ../results/output.wis_via_iteration.100k.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.2884  0.4100  1.0608  2.2542  3.6608  9.9204 
```