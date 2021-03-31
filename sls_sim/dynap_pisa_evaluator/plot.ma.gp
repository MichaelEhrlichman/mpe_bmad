set term postscript enhanced color solid size 10,3.5
set output 'ma.eps'

set xlabel 'Location (m)'
set ylabel 'Momentum Aperture (%)'

set key box opaque maxrows 2

#set key right center
set key right top

#set yrange[-15:20]
#set xrange[0:96]
set xrange[0:24.2]

plot '00touschek/aperture.dat' u 2:($3*100) t '6D Momentum Aperture' w l lt 2 lc rgb 'blue',\
'' u 2:($4*100) notitle w l lt 2 lc rgb 'blue',\
'00touschek/linear_ma.dat' u 2:($3*100) t 'Linear Aperture' w l lt 2 lc rgb 'magenta',\
'' u 2:($4*100) notitle w l lt 2 lc rgb 'magenta'

