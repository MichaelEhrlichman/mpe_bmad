#!/bin/bash

for i in `find $PWD -type d -name "seed_*"`
do 
  cd $i
  ~/bbin/orm ../orm.in seed_results_*
  ~/bbin/prm ../prm.in orbit_correction.lat
done
