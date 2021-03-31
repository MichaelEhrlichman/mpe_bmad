set term postscript enhanced color solid size 10,5
#set size 0.7,0.7
set output 'adts_x.eps'

set xlabel 'x (mm)'
set ylabel 'Horizontal Tune'

plot '00adts/tracker_adts.out' i 0 u ($1*1000):($3>0?$3:1./0.) notitle w l lt 2 lc rgb 'blue'

set output 'adts_y.eps'

set xlabel 'y (mm)'
set ylabel 'Vertical Tune'

plot '00adts/tracker_adts.out' i 1 u ($2*1000):($4>0?$4:1./0.) notitle w l lt 2 lc rgb 'blue'

