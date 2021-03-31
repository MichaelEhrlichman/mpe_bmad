set term postscript enhanced color solid size 10,5
set output 'co.eps'

set ylabel 'x (mm)'
set xlabel 'Location (m)

set title '4D + pz closed orbit at +-4%'

plot '00co/co_-0.04.dat' u 1:($2*1000) t 'pz = -0.04' w l,\
'00co/co_0.04.dat' u 1:($2*1000) t 'pz = 0.04' w l




