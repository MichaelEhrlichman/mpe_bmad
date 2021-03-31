set term postscript enhanced color font 18

set output 'pretty.da.eps'

#set colorsequence default

set xlabel 'x (mm)'
set ylabel 'y (mm)'

set yrange[0:]

set key maxrows 4 opaque box

plot '00da/da.out' i 0 u ($2*1000):($3*1000) t ' 0% DA' w lp lt 2 lc rgb 'red',\
'' i 1 u ($2*1000):($3*1000) t '-3% DA' w lp lt 2 lc rgb 'blue',\
'' i 2 u ($2*1000):($3*1000) t '+3% DA' w lp lt 2 lc rgb 'green',\
'/afs/psi.ch/user/e/ehrlichman_m/non_bmad_code/circle_data/circle_data.dat' u ($1*10):($2*10) t 'chamber' w l lt -1,\
'00da_linear/da.out' i 0 u ($2*1000):($3*1000) t '0% LA' w l lt 1 lc rgb 'red',\
'' i 1 u ($2*1000):($3*1000) t '-3% LA' w l lt 1 lc rgb 'blue',\
'' i 2 u ($2*1000):($3*1000) t '+3% LA' w l lt 1 lc rgb 'green' 

