#!/usr/bin/env gnuplot

set term eps noenhanced size 15,10.5 font ",12"
set output 'avg_con_auto.eps'

set multiplot layout 4,4

stats '<tail -n +2 constraint_report.avg' nooutput
ncol = STATS_columns
set key autotitle columnheader bottom right

set xrange[40:]

do for [i=1:ncol:2] {
plot 'constraint_report.avg' u i w d,\
'' u (column(i+1)>0.0?0:column(i+1)) notitle w d
}


