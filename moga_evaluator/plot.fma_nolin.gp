set term postscript enhanced color size 10,15
set output 'fma_nolin.eps'

set multiplot layout 3,1

set palette rgbformulae 22,13,-31

set xlabel 'x (mm)'
set ylabel 'y (mm)'

#set xrange[-10:10]
#set yrange [0:10]
set cbrange[-20:-5]
set autoscale fix

key1=system("grep '^1' 00da_raster/raster.key | awk '{print $2*100}'")
key2=system("grep '^2' 00da_raster/raster.key | awk '{print $2*100}'")
key3=system("grep '^3' 00da_raster/raster.key | awk '{print $2*100}'")

set title key1.'%'
plot '00da_raster/raster_1fma.dat' nonuniform matrix using ($2*1000):($1*1000):($3>-90 ? $3 : NaN) notitle with image

set title key2.'%'
plot '00da_raster/raster_2fma.dat' nonuniform matrix using ($2*1000):($1*1000):3 notitle with image

set title key3.'%'
plot '00da_raster/raster_3fma.dat' nonuniform matrix using ($2*1000):($1*1000):3 notitle with image

