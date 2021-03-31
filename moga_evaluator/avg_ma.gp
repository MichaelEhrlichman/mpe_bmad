#!/usr/bin/env gnuplot

set term eps size 15,10.5 font ",24"
set output 'avg_obj.eps'

system("~/scripts/avg_by_index.py objective_report.out > objective_avg.out")

set multiplot layout 2,2

set xrange[29:]

plot 'objective_avg.out' u 9 title 'emittance' w d,\
'' u 10 notitle w d,\
0.5 t 'v20r'
#-4.96605512E-04 t 'v20r'

plot 'objective_avg.out' u 3 title 'DA 0%' w d,\
'' u 4 notitle w d
#8.73809542E-01 t 'v20r'

plot 'objective_avg.out' u 7 title 'pos ma' w d,\
'' u 8 notitle w d,\
0.7662 t 'v20r'

plot 'objective_avg.out' u 5 title 'neg ma' w d,\
'' u 6 notitle w d,\
0.7375 t 'v20r'

