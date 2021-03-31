set term postscript enhanced color font 18

set output 'pretty.mada.eps'

#set colorsequence default

set xlabel 'x (mm)'
set ylabel 'y (mm)'

set yrange[0:]

set key maxrows 4 opaque box

list=system('ls seed_*/00da/da.out')

plot '/afs/psi.ch/user/e/ehrlichman_m/non_bmad_code/circle_data/circle_data.dat' u ($1*10):($2*10) t 'Chamber' w l lt -1,\
'../00da_linear/da.out' i 0 u ($2*1000):($3*1000) t 'On-energy LA' w l lt 1 lc rgb 'blue',\
'../00da/da.out' i 0 u ($2*1000):($3*1000) t ' On-energy DA' w lp lt 2 lc rgb 'red',\
for [file in list] file i 0 u ($2*1000):($3*1000) notitle w lp lt 2 lc rgb 'green',\
NaN title 'Misaligned DA' w lp lt 2 lc rgb 'green'

