#!/usr/bin/env gnuplot

#set term png size 1900,960 font "Arial,10"
set term png size 1900,1060 font "Arial,10"
set output 'vars.png'

datafile='moga_results.out'

countCommand='head -3 '.datafile.' | tail -1 | wc -w'
ncols = system(countCommand)
stats datafile

nplot=100

ngen = int(STATS_blocks-1)
stepsize = ngen/nplot

# fG = ngen%stepsize
# fG = (fG==0)?stepsize:fG
# nGenToPlot = ceil(1.0*ngen/stepsize)

set style fill transparent solid 0.10 noborder
set style circle radius 0.12

set xrange[0:nplot+1]

set multiplot layout 4,5

set nokey

do for [col=2:ncols-5] {set title system("awk '/^#/ {next}; {print $".(col+1).";exit}' ".datafile); plot for [IDX=0:nplot-1] datafile i IDX*stepsize u (IDX):(column(col)) w circles lc rgb 'blue' }

set key center box font "Arial,6"
set notitle
unset border
unset xtics
unset ytics

plot [0:1] [0:1] for [IDX=0:nplot-1] NaN t sprintf("%2d:  %5d",IDX,IDX*stepsize) 
