set term eps enhanced font "Arial,12" size 7,3

set output "cdf.eps"

a=1./30.

set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1
#set multiplot layout 2,1 margins 0.10,0.95,.15,.95 spacing 0,0
set multiplot

set xrange[0:3.5]
unset xtics
unset ytics
unset key
set size 0.85,0.05
set origin 0.1,0.9
set y2tics ("ideal" 1)
plot 'candidates.dat' u 1:(1) t 'v20r' lw 3 ps 0.6 lc rgb 'black' dt 1,\
'' u 2:(1) t '380' lw 3 ps 0.6 lc rgb 'black' dt 2,\
'' u 3:(1) t '697301' lw 1 ps 0.6,\
'' u 4:(1) t '421' lw 1 ps 0.6,\
'' u 5:(1) t '541' lw 1 ps 0.6,\
'' u 6:(1) t '1085653' lw 1 ps 0.6,\
'' u 7:(1) t '123925' lw 1 ps 0.6,\
'' u 8:(1) t '171' lw 1 ps 0.6,\
'' u 9:(1) t '505097' lw 1 ps 0.6,\
'' u 10:(1) t '994971' lw 1 ps 0.6

set xlabel 'Touschek Lifetime (hr)'
set ylabel 'CDF (30 seeds)'
set key top left
unset y2tics
set xtics
set ytics
set size 0.85,0.75
set origin 0.1,0.15
plot 'errors_seed_0/tl-6D.dat' u 1:(1./26.) smooth cumul title 'v20r' w lp lw 3 ps 0.2 lc rgb 'black' dt 1,\
'errors_seed_special380/tl-6D.dat' u 1:(a) smooth cumul title '380' w lp lw 3 ps 0.2 lc rgb 'black' dt 2,\
'errors_seed_697301/tl-6D.dat' u 1:(a) smooth cumul title 'no har 697301' w lp lw 1 ps 0.2,\
'errors_seed_421/tl-6D.dat' u 1:(a) smooth cumul title 'no har 421' w lp lw 1 ps 0.2,\
'errors_seed_541/tl-6D.dat' u 1:(a) smooth cumul title 'no har 541' w lp lw 1 ps 0.2,\
'errors_seed_1085653/tl-6D.dat' u 1:(a) smooth cumul title 'no har 1085653' w lp lw 1 ps 0.2,\
'errors_seed_123925/tl-6D.dat' u 1:(a) smooth cumul title 'no har 123925' w lp lw 1 ps 0.2,\
'errors_seed_171/tl-6D.dat' u 1:(a) smooth cumul title 'no har 171' w lp lw 1 ps 0.2,\
'errors_seed_505097/tl-6D.dat' u 1:(a) smooth cumul title 'no har 505097' w lp lw 1 ps 0.2,\
'errors_seed_994971/tl-6D.dat' u 1:(a) smooth cumul title 'no har 994971' w lp lw 1 ps 0.2

