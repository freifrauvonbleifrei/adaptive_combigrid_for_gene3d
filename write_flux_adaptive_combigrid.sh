#!/bin/bash

. ./adaptation_parameters.sh

MODE="oldSet.csv"
if [ ! -z $1 ]; then
    MODE=$1
fi

STARTTIME=${ADAPTATION_PARAM_CROP_TIME}
if [ ! -z $2 ]; then
    STARTTIME=$2
fi

ENDTIME=${ADAPTATION_PARAM_MIN_RUNTIME}
if [ ! -z $3 ]; then
    ENDTIME=$3
fi

python3 write_flux_adaptive_combigrid.py $MODE $STARTTIME $ENDTIME
