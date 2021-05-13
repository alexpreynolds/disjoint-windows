SHELL=/bin/bash

DATA=${PWD}/../data
RESULTS=${PWD}/../results

INPUT_FN=${DATA}/windows.fixed.25k.bed

LEFT_WINDOW_WIDTH=12000
RIGHT_WINDOW_WIDTH=12000
EXTEND_WIDTH=12000
WINDOW_WIDTH=25000
INTERVAL_STEP=1000
SAMPLES=100000

all: extend_windows

extend_windows:
	${PWD}/extend_windows.sh ${LEFT_WINDOW_WIDTH} ${RIGHT_WINDOW_WIDTH} ${WINDOW_WIDTH} ${DATA}/removeOverlap.bed.starch ${INPUT_FN}

recommender:
	unstarch ${DATA}/removeOverlap.bed.starch > ${DATA}/removeOverlap.bed
	${PWD}/recommender_subset.py --input=${DATA}/removeOverlap.bed --k=${SAMPLES} --extend=${EXTEND_WIDTH} --window-size=${INTERVAL_STEP}

wis_via_iteration:
	${PWD}/weighted_interval_scheduling_via_iteration.py --input=${INPUT_FN} | sort-bed - > ${RESULTS}/output.wis_via_iteration.all.bed
	${PWD}/weighted_interval_scheduling_via_iteration.py --input=${INPUT_FN} --k=${SAMPLES} | sort-bed - > ${RESULTS}/output.wis_via_iteration.100k.bed

wis_via_iteration_summaries:
	cut -f4 ${INPUT_FN} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
	cut -f4 ${RESULTS}/output.wis_via_iteration.all.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
	cut -f4 ${RESULTS}/output.wis_via_iteration.100k.bed | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

priority_queue:
	${PWD}/max_heap.py --input=${INPUT_FN} --k=${SAMPLES} | sort-bed - > ${RESULTS}/output.max_heap.100k.bed

priority_queue_with_jitter:
	${PWD}/max_heap_with_jitter.py --input=${INPUT_FN} --k=${SAMPLES} | sort-bed - > ${RESULTS}/output.max_heap_with_jitter.100k.bed

walkers:
	${PWD}/walkers.py --input=${INPUT_FN} --k=${SAMPLES} | sort-bed - > ${RESULTS}/output.walkers.100k.bed

#
# test Walker's Alias method
#
# keys           [1, 2, 3, 4, 5] 
# weights        [0, 2, 4, 2, 2] -- note that this doesn't need to sum to 1
# sample size    100000
#
# resulting frequency table of per-key samples should reflect (normalized) input weights
#

TEST_SAMPLES=1000000

walkers_test:
	${PWD}/walkers_test.py ${TEST_SAMPLES} | awk -v N=${TEST_SAMPLES} '{a[$$1]++}END{for(n in a){ print n, a[n], a[n]/N }}' | sort -n
