set term postscript enhanced color size 6,6
set output 'fma_pz.eps'

set palette rgbformulae 22,13,-31

set xlabel 'x (mm)'
set ylabel 'pz (%)'

#set xrange[-10:10]
#set yrange [0:10]
set cbrange[-20:-5]

plot '00xpz_raster/raster_xpz_fma.dat' nonuniform matrix using ($2*1000):($1*100):($3>-90 ? $3 : NaN) notitle with image,\
'00da_linear/da_pz.out' u ($2*1000):($1*100) title 'physical aperture' w l lt 8 lw 3,\
'' u ($3*1000):($1*100) notitle w l lt 8 lw 3

