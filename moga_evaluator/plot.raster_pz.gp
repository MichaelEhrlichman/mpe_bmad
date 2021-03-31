set term postscript enhanced color size 6,6
set output 'raster_pz.eps'

set palette defined (0.5 0.5 0.5 0.5, 1 1 1 1)

set xlabel 'x (mm)'
set ylabel 'pz (%)'

set xrange[-6:6]
set yrange [-7:7]

set title 'x-pz aperture'
unset colorbox

plot '00xpz_raster/raster_xpz.dat' nonuniform matrix using ($2*1000):($1*100):3 notitle with image

#'00da_linear/da_pz.out' u ($2*1000):($1*100) notitle w l lt 1 lc rgb 'blue',\
#'00da_linear/da_pz.out' u ($3*1000):($1*100) notitle w l lt 1 lc rgb 'blue'

