#!/usr/bin/env gnuplot

set term eps size 15,10.5 font ",14"
set output 'avg_con.eps'

set multiplot layout 4,4

set xrange[10:]

plot 'constraint_report.avg' u 1 title 'mag str' w d,\
'' u 2 notitle w d

plot 'constraint_report.avg' u 3 title 'mat -' w d,\
'' u 4 notitle w d

plot 'constraint_report.avg' u 5 title 'mat +' w d,\
'' u 6 notitle w d

plot 'constraint_report.avg' u 7 title 'co -' w d,\
'' u 8 notitle w d

plot 'constraint_report.avg' u 9 title 'co +' w d,\
'' u 10 notitle w d

plot 'constraint_report.avg' u 11 title 'wk pt x' w d,\
'' u 12 notitle w d

plot 'constraint_report.avg' u 13 title 'wk pt y' w d,\
'' u 14 notitle w d

plot 'constraint_report.avg' u 15 title 'ID beta x' w d,\
'' u 16 notitle w d

plot 'constraint_report.avg' u 17 title 'ID beta y' w d,\
'' u 18 notitle w d

plot 'constraint_report.avg' u 19 title 'ID eta x' w d,\
'' u 20 notitle w d

plot 'constraint_report.avg' u 21 title 'global beta x' w d,\
'' u 22 notitle w d

plot 'constraint_report.avg' u 23 title 'coupling' w d,\
'' u 24 notitle w d

plot 'constraint_report.avg' u 25 title 'chrom x' w d,\
'' u 26 notitle w d

plot 'constraint_report.avg' u 27 title 'chrom y' w d,\
'' u 28 notitle w d
