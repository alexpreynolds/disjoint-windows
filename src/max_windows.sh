#!/bin/bash

module add bedops

WIDTH=${1}
STEP=${2}
IN=${3}
OUT=${4}

#
# chop up the genomic space by specified width, then map max signal from modified input to this space
#

bedmap --echo --max --delim '\t' <( bedops --merge ${IN} | bedops --chop ${WIDTH} --stagger ${STEP} - ) <( awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, ".", $4 }' ${IN} ) > ${OUT}
