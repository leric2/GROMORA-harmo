#!/bin/bash

Yesterday_date=$(date -d "yesterday" '+%Y-%m-%d')
cd calibration/
echo "calibration_GROMORA"| matlab -nodisplay -nodesktop -nosplash >> /storage/tub/instruments/gromos/level1/GROMORA/calibration_log/cal_log_${Yesterday_date}.txt
