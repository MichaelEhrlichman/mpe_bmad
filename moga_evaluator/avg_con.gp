#!/usr/bin/env gnuplot

set term eps size 15,10.5 font ",14"
set output 'avg_con.eps'

system("~/scripts/avg_by_index_con.py constraint_report.out > constraint_avg.out")

set multiplot layout 4,4

set xrange[500:]

plot 'constraint_avg.out' u 3 title 'mag str' w d,\
'' u 4 notitle w d

plot 'constraint_avg.out' u 5 title 'mat x' w d,\
'' u 6 notitle w d

plot 'constraint_avg.out' u 7 title 'mat y' w d,\
'' u 8 notitle w d

plot 'constraint_avg.out' u 9 title 'co -' w d,\
'' u 10 notitle w d

plot 'constraint_avg.out' u 11 title 'co +' w d,\
'' u 12 notitle w d

plot 'constraint_avg.out' u 13 title 'wk pt x' w d,\
'' u 14 notitle w d

plot 'constraint_avg.out' u 15 title 'wk pt y' w d,\
'' u 16 notitle w d

plot 'constraint_avg.out' u 17 title 'ID beta x' w d,\
'' u 18 notitle w d

plot 'constraint_avg.out' u 19 title 'ID beta y' w d,\
'' u 20 notitle w d

plot 'constraint_avg.out' u 21 title 'ID eta x' w d,\
'' u 22 notitle w d

plot 'constraint_avg.out' u 23 title 'global beta x' w d,\
'' u 24 notitle w d

plot 'constraint_avg.out' u 25 title 'coupling' w d,\
'' u 26 notitle w d

plot 'constraint_avg.out' u 27 title 'chrom x' w d,\
'' u 28 notitle w d

plot 'constraint_avg.out' u 29 title 'chrom y' w d,\
'' u 30 notitle w d
