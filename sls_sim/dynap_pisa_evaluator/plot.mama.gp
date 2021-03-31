set term postscript enhanced color solid size 10,3.5
set output 'mama.eps'

set xlabel 'Location (m)'
set ylabel 'Momentum Aperture (%)'

set key box opaque #maxrows 2

set yrange[-15:20]
set xrange[0:96]

list=system('ls seed_*/00touschek6D/aperture.dat')

plot '../00touschek/linear_ma_hd.dat' u 2:($3*100) t 'Linear' w l lt 2 lc rgb 'black',\
'' u 2:($4*100) notitle w l lt 2 lc rgb 'black',\
for [file in list] file u 2:($3*100) notitle w l lt 2 lc rgb 'magenta',\
for [file in list] file u 2:($4*100) notitle w l lt 2 lc rgb 'magenta',\
NaN t 'MA & Cor.' w l lt 2 lc rgb 'magenta',\
'../00touschek6D/aperture.dat' u 2:($3*100) t 'Ideal' w l lt 2 lc rgb 'blue',\
'' u 2:($4*100) notitle w l lt 2 lc rgb 'blue',\
5 t '5%' w l lt 2 lc rgb 'green',\
-5 notitle w l lt 2 lc rgb 'green'

