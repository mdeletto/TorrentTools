#!/bin/bash


echo "hey"
#echo ${RESULTS_DIR}[17:]
python ${DIRNAME}/Uniformity.py ${ANALYSIS_DIR} ${RESULTS_DIR} ${TSP_ANALYSIS_NAME}
