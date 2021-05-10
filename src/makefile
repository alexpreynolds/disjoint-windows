SHELL=/bin/bash

DATA=${PWD}/../data
RESULTS=${PWD}/../results

WINDOW_WIDTH=25000
WINDOW_STEP=1000

all: max_windows

max_windows:
	${PWD}/max_windows.sh ${WINDOW_WIDTH} ${WINDOW_STEP} ${DATA}/removeOverlap.bed ${DATA}/windows.max.25k.bed

wis_via_iteration:
	${PWD}/weighted_interval_scheduling_via_iteration.py --input=${DATA}/windows.max.25k.bed --k=100000 | sort-bed - > ${RESULTS}/output.wis_via_iteration.bed

priority_queue:
	${PWD}/max_heap.py --input=${DATA}/windows.max.25k.bed --k=100

walkers:
	${PWD}/walkers.py --input=${DATA}/windows.max.25k.bed --k=100

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
