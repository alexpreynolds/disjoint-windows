# disjoint-windows

We have the following goals:

1. Find a good, potentially optimal distribution of scored intervals over the genome, which should be disjoint, non-overlapping or, in other settings, called "mutually compatible".

2. We would like the selected method to prioritize high-scoring intervals.

3. We would like 100k such elements out of the full superset, which meet the first two criteria.

## Methods

Review targets in `src/makefile` for more concrete run instructions.

### `max_windows`

We start with 1kb windows spanning the genome. Each window has a score.

The `max_windows` target creates 25kb windows, which step across the genome in 1kb units. Each interval contains the maximum score over the 25kb span. There are ~3M such elements.

These windows are contained in `data/windows.max.25k.bed`.

### `walkers`

In the `walkers` target, we perform weighted sampling with replacement using the [Walker's Alias](https://en.wikipedia.org/wiki/Alias_method) technique. We use the score of each input element as its weight, so samples are high-scoring. We reject samples which overlap the genomic space of elements already sampled.

We repeat sampling until we get 100k elements, or until we hit some iteration limit, as in this case it is unlikely we'll find any more qualifying elements.

### `priority_queue`

In the `priority_queue` target, we put elements into a priority queue, ordered by score. We pop the highest-scoring element off the queue, keep it if it does not overlap any other popped elements, and reject it if it does. 

We repeat this until the queue is empty, or until we get 100k elements.

### `wis_via_iteration`

In the `wis_via_iteration` target, we use a dynamic programming approach to build a list of high-scoring elements which do not overlap one another.

The locally optimal path contains elements which do not overlap. If we have more than 100k such elements in this path, we use a priority queue to get the top-100k such elements.

## Results

### Bad news

The `walkers` and `priority_queue` elements prioritize getting high-scoring elements over getting more elements across the genomic space. 

These "greedy" methods stop return elements around 50k elements, and so do not meet the 100k criterion. In the case of Walker's, it is likely that sampling with replacement will bias samples towards picking the same high-scoring elements, which will always get rejected after the first insertion. The heap method, on the other hand,does not sample randomly, but picks the best elements first. These high-value intervals tend to clump together. Those that are subsequently popped will be likely get rejected, being very close to best elements that have been inserted.

### Good news

The weighted-interval scheduling approach of `wis_via_iteration` generates ~112k non-overlapping elements over our test input, and so we are able to return 100k non-overlapping intervals.

We can look at the signal distribution in the original windows:

```
$ cut -f4 ../data/windows.max.25k.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1354  0.4100  1.1727  2.3122  3.7394  9.9204 
```

We would expect (and we observe) that the ~112k `wis_via_iteration` subset has a better signal profile, as we are picking high-score intervals:

```
$ cut -f4 ../results/output.all.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1673  0.5035  2.0201  2.8866  4.8378  9.9204 
```

As expected, the use of a priority queue to filter these down to 100k elements returns a better overall signal:

```
$ cut -f4 ../results/output.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.4046  0.8628  2.6131  3.2338  5.2530  9.9204 
```