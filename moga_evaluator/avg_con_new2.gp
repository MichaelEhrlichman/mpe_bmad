#!/usr/bin/env gnuplot

set term eps size 15,10.5 font ",14"
set output 'avg_con.eps'

set multiplot layout 4,4

set xrange[20:]

plot 'constraint_report.avg' u 'max|k2|[avg]' title 'mag str' w d,\
'' u  'max|k2|[best]' notitle w d

plot 'constraint_report.avg' u 'mats_-[avg]' title 'mat -' w d,\
'' u 'mats_-[best]' notitle w d

plot 'constraint_report.avg' u 'mats_+[avg]' title 'mat +' w d,\
'' u 'mats_+[best]' notitle w d

plot 'constraint_report.avg' u 'co@-de[avg]' title 'co -' w d,\
'' u 'co@-de[best]' notitle w d

plot 'constraint_report.avg' u 'co@+de[avg]' title 'co +' w d,\
'' u 'co@+de[best]' notitle w d

plot 'constraint_report.avg' u 'wk_pt_x[avg]' title 'wk pt x' w d,\
'' u 'wk_pt_x[best]' notitle w d

plot 'constraint_report.avg' u 'wk_pt_y[avg]' title 'wk pt y' w d,\
'' u 'wk_pt_y[best]' notitle w d

plot 'constraint_report.avg' u 'beta_x[avg]' title 'ID beta x' w d,\
'' u 'beta_x[best]' notitle w d

plot 'constraint_report.avg' u 'beta_y[avg]' title 'ID beta y' w d,\
'' u 'beta_y[best]' notitle w d

plot 'constraint_report.avg' u 'eta_x[avg]' title 'ID eta x' w d,\
'' u 'eta_x[best]' notitle w d

plot 'constraint_report.avg' u 'glo_beta_x[avg]' title 'global beta x' w d,\
'' u 'glo_beta_x[best]' notitle w d

plot 'constraint_report.avg' u 'coupling[avg]' title 'coupling' w d,\
'' u 'coupling[best]' notitle w d

plot 'constraint_report.avg' u 'chrom_x[avg]' title 'chrom x' w d,\
'' u 'chrom_x[best]' notitle w d

plot 'constraint_report.avg' u 'chrom_y[avg]' title 'chrom y' w d,\
'' u 'chrom_y[best]' notitle w d
