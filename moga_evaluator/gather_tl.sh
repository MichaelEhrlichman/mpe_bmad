#!/usr/bin/env bash

grep -v hours seed_*/00touschek/tl-6D.dat | awk '{print $2}' > tl-6D.dat
wc -l tl-6D.dat
~/bmad_dist_2019_0521/mpe_bmad/moga_evaluator/col_avg.gp tl-6D.dat 1

