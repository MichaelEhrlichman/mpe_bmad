set term postscript enhanced color size 10,15
set output 'raster.eps'

set multiplot layout 3,1

set palette defined (0.5 0.5 0.5 0.5, 1 1 1 1)

set xlabel 'x (mm)'
set ylabel 'y (mm)'

set xrange[-30:30]
set yrange [0:16]

key1=system("grep '^1' 00da_raster/raster.key | awk '{print $2*100}'")
key2=system("grep '^2' 00da_raster/raster.key | awk '{print $2*100}'")
key3=system("grep '^3' 00da_raster/raster.key | awk '{print $2*100}'")

set title key1.'%'
plot '00da_raster/raster_1.dat' nonuniform matrix using ($2*1000):($1*1000):3 notitle with image,\
'00da_linear/da.out' i 0 u ($2*1000):($3*1000) notitle w l lt 1 lc rgb 'blue'

set title key2.'%'
plot '00da_raster/raster_2.dat' nonuniform matrix using ($2*1000):($1*1000):3 notitle with image,\
'00da_linear/da.out' i 1 u ($2*1000):($3*1000) notitle w l lt 1 lc rgb 'blue'

set title key3.'%'
plot '00da_raster/raster_3.dat' nonuniform matrix using ($2*1000):($1*1000):3 notitle with image,\
'00da_linear/da.out' i 2 u ($2*1000):($3*1000) notitle w l lt 1 lc rgb 'blue'

# '00da_ndp/raster.grid' i 0 u ($1*1000):($2*1000) notitle w p lt 1 lc rgb 'red',\
# '00da_ndp/raster.grid' i 1 u ($1*1000):($2*1000) notitle w p lt 1 lc rgb 'red',\
# '00da_ndp/raster.grid' i 2 u ($1*1000):($2*1000) notitle w p lt 1 lc rgb 'red',\
