#!/bin/bash

for i in `find $PWD -maxdepth 1 -type d -name "seed_*"`
do 
  cd $i
  ~/bbin/linear_optics seed_results_*
  ~/bbin/beta_beat ../../girders.bmad seed_results_*
done
