#!/bin/bash

# create histograms
function unpack_data(){
    for file in $DATA_DIR/analysis"$1"_*.root; do
        fullname=${file##*/}
        filename=${fullname%%.*}
        subrun=${filename##*_}
        hist_file=$HIST_DIR/run_$1_$subrun.root
        if [[ ! -f $hist_file ]]; then
             # echo "$SORT_CODE $cal_file $file"
             cd $SORT_DIR
             $SORT_CODE $CAL_FILE $file
             mv histograms.root $hist_file
        elif [[ -f $hist_file ]]; then
            echo "Run $1 subrun $subrun histograms already exist, skipping ..."
        fi
    done
}

#--- Main Function ---------------------------------------
TOP_DIR=$(pwd)
#CAL_FILE="/data1/S9038/current-sort/data/cal-files/run-agnostic/July2020_calibration.cal"
#HIST_DIR="/data1/S9038/current-sort/data/histograms/compton-algorithm/subruns"
#DATA_DIR="/data1/S9038/current-sort/data/analysis-trees"
#SORT_CODE="/data1/S9038/current-sort/analysis/sort-analysis-trees/my-build/SortAnalysisTrees"
#SORT_DIR="/data1/S9038/current-sort/ongoing-sort"

CAL_FILE="/data1/S9038/current-sort/data/S2036/cal-files/CalibrationFile.cal"
HIST_DIR="/data1/S9038/current-sort/data/S2036/histograms"
DATA_DIR="/data1/S9038/current-sort/data/S2036/analysis-trees"
SORT_CODE="/data1/S9038/current-sort/analysis/sorting/bgo-suppression-factor/my-build/SortAnalysisTrees"
SORT_DIR="/data1/S9038/current-sort/ongoing-sort"

if [ $# -eq 0 ]; then
    echo "Pass one arguement to sort a single run or two for a range"
elif [ $# -eq 1 ]; then
    unpack_data $1
elif [ $# -eq 2 ]; then
    for ((i=$1;i<=$2;i++)); do
        unpack_data $i
    done
fi
