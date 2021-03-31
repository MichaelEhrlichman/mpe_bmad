set term postscript enhanced color 
set output 'da.eps'

#set colorsequence default


plot '00da/da.out' i 0 u 2:3 t '0% DA' w lp lt 2 lc rgb 'red',\
'' i 1 u 2:3 t '-dE' w lp lt 2 lc rgb 'blue',\
'' i 2 u 2:3 t '+dE' w lp lt 2 lc rgb 'green',\
'00da_linear/da.out' i 0 u 2:3 t '0% LA' w l lt 1 lc rgb 'red',\
'' i 1 u 2:3 t '-dE' w l lt 1 lc rgb 'blue',\
'' i 2 u 2:3 t '+dE' w l lt 1 lc rgb 'green'
