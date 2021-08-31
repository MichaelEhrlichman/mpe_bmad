plot 'lik.dat' u 1:2 t 'likelihood' w l,\
'meas.dat' u 1:2 t 'measurements' w p lt 7 ps 2,\
'pr.dat' u 1:2 t 'prior distribution' w l,\
'po.dat' u 1:2 t 'posterior distribution' w l

pause -1

