set term postscript enhanced color solid size 12.0,6.0
#set size 0.7,0.7
set size square
set output 'trace.eps'

set multiplot layout 1,2

set xlabel 'tr(x)'
set ylabel 'tr(y)'

set xrange[-2.5:2.5]
set yrange[-2.5:2.5]

set object 1 rect from -2,-2 to 2,2 

plot '00footprint/footprint.dat' u 4:5 notitle w l lt 2 lc rgb 'blue'

set auto x
set xlabel '% Energy Offset'
set ylabel 'trace'
set yrange [-2.5:2.5]
unset object 1

plot '00footprint/footprint.dat' u ($1*100):4 t 'x' w l lt 2 lc rgb 'blue',\
'' u ($1*100):5 t 'y' w l lt 2 lc rgb 'red',\
-2 notitle lt 0,\
2 notitle lt 0
