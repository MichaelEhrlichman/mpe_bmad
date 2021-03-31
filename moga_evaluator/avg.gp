set term eps size 15,10.5 font ",26"
set output 'avg_obj.eps'

system("~/scripts/avg_by_index.py objective_report.out > objective_avg.out")

set multiplot layout 2,2

plot 'objective_avg.out' u 9 title 'emittance' w d,\
'' u 10 notitle w d

plot 'objective_avg.out' u 3 title 'DA 0%' w d,\
'' u 4 notitle w d

plot 'objective_avg.out' u 7 title 'DA 3%' w d,\
'' u 8 notitle w d

plot 'objective_avg.out' u 5 title 'DA -3%' w d,\
'' u 6 notitle w d

