#!/bin/bash

left_window_width=${1}
right_window_width=${2}
window_width=${3}
in_starch_fn=${4}
out_bed_fn=${5}

module add bedops

unstarch ${in_starch_fn} \
    | bedops --range -${left_window_width}:${right_window_width} --everything - \
    | awk -v w=${window_width} '(($3-$2) == w)' \
    > ${out_bed_fn}