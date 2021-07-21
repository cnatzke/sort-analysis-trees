#!/bin/bash

# create histograms
function unpack_data(){
    for file in $DATA_DIR/analysis"$1"_*.root; do
        fullname=${file##*/}
        filename=${fullname%%.*}
        subrun=${filename##*_}
        hist_file=$HIST_DIR/run_$1_$subrun.root
        if [[ ! -f $hist_file ]]; then
             #echo "$SORT_CODE $CAL_FILE $file"
             $SORT_CODE $CAL_FILE $file
             mv histograms.root $hist_file
        elif [[ -f $hist_file ]]; then
            echo "Run $1 subrun $subrun histograms already exist, skipping ..."
        fi
    done
}

#--- Main Function ---------------------------------------
TOP_DIR=$(pwd)
CAL_FILE="$TOP_DIR/July2020_calibration.cal"
HIST_DIR="/data/histograms/subruns"
DATA_DIR="/data/analysis-trees"
SORT_CODE="/vagrant/sort-analysis-tree/my-build/SortAnalysisTrees"

if [ $# -eq 0 ]; then
    echo "Pass one arguement to sort a single run or two for a range"
elif [ $# -eq 1 ]; then
    unpack_data $1
elif [ $# -eq 2 ]; then
    for ((i=$1;i<=$2;i++)); do
        unpack_data $i
    done
fi
