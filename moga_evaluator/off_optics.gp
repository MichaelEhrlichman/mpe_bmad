#!/usr/bin/env gnuplot

set term eps size 10,5 font "Arial,22"

set xrange[0:16.375]

set xlabel 'Location (s)'

set key above

set output 'off_optics_Bx.eps'
col = 4
set ylabel '{/Symbol b}_x (m)'
plot 'off_twiss_0.01.out'     u 2:col title '1%' w l dt 1 lc rgb 'blue' lw 2,\
'off_twiss_-0.01.out'    u 2:col title '-1%' w l dt 2 lc rgb 'blue' lw 2,\
'off_twiss_0.02.out'     u 2:col title '2%' w l dt 1 lc rgb 'red' lw 2,\
'off_twiss_-0.02.out'    u 2:col title '-2%' w l dt 2 lc rgb 'red' lw 2,\
'off_twiss_0.03.out'     u 2:col title '3%' w l dt 1 lc rgb 'green' lw 2,\
'off_twiss_-0.03.out'    u 2:col title '-3%' w l dt 2 lc rgb 'green' lw 2,\
'off_twiss_0.0325.out'     u 2:col title '3.25%' w l dt 1 lc rgb 'brown' lw 2,\
'off_twiss_-0.0325.out'    u 2:col title '-3.25%' w l dt 2 lc rgb 'brown' lw 2,\
'off_twiss_0.0.out' u 2:col title 'on-energy' w l lt -1 lw 2

set output 'off_optics_nx.eps'
set ylabel '{/Symbol h}_x (cm)'
plot 'off_twiss_0.01.out'     u 2:($7*1e2) title '1%' w l dt 1 lc rgb 'blue' lw 2,\
'off_twiss_-0.01.out'    u 2:($7*1e2) title '-1%' w l dt 2 lc rgb 'blue' lw 2,\
'off_twiss_0.02.out'     u 2:($7*1e2) title '2%' w l dt 1 lc rgb 'red' lw 2,\
'off_twiss_-0.02.out'    u 2:($7*1e2) title '-2%' w l dt 2 lc rgb 'red' lw 2,\
'off_twiss_0.03.out'     u 2:($7*1e2) title '3%' w l dt 1 lc rgb 'green' lw 2,\
'off_twiss_-0.03.out'    u 2:($7*1e2) title '-3%' w l dt 2 lc rgb 'green' lw 2,\
'off_twiss_0.0325.out'     u 2:($7*1e2) title '3.25%' w l dt 1 lc rgb 'brown' lw 2,\
'off_twiss_-0.0325.out'    u 2:($7*1e2) title '-3.25%' w l dt 2 lc rgb 'brown' lw 2,\
'off_twiss_0.0.out' u 2:($7*1e2) title 'on-energy' w l lt -1 lw 2

set output 'off_optics_By.eps'
col = 12
set ylabel '{/Symbol b}_y (cm)'
plot 'off_twiss_0.01.out'     u 2:col title '1%' w l dt 1 lc rgb 'blue' lw 2,\
'off_twiss_-0.01.out'    u 2:col title '-1%' w l dt 2 lc rgb 'blue' lw 2,\
'off_twiss_0.02.out'     u 2:col title '2%' w l dt 1 lc rgb 'red' lw 2,\
'off_twiss_-0.02.out'    u 2:col title '-2%' w l dt 2 lc rgb 'red' lw 2,\
'off_twiss_0.03.out'     u 2:col title '3%' w l dt 1 lc rgb 'green' lw 2,\
'off_twiss_-0.03.out'    u 2:col title '-3%' w l dt 2 lc rgb 'green' lw 2,\
'off_twiss_0.0325.out'     u 2:col title '3.25%' w l dt 1 lc rgb 'brown' lw 2,\
'off_twiss_-0.0325.out'    u 2:col title '-3.25%' w l dt 2 lc rgb 'brown' lw 2,\
'off_twiss_0.0.out' u 2:col title 'on-energy' w l lt -1 lw 2

set output 'off_optics_xco_minus_eta.eps'
set ylabel '(x_{co}-{/Symbol h}_x*{/Symbol d}E) (mm)'
plot 'off_twiss_0.01.out'     u 2:(($20-$7*0.01)*1e3) title '1%' w l dt 1 lc rgb 'blue' lw 2,\
'off_twiss_-0.01.out'    u 2:(($20-$7*-0.01)*1e3) title '-1%' w l dt 2 lc rgb 'blue' lw 2,\
'off_twiss_0.02.out'     u 2:(($20-$7*0.02)*1e3) title '2%' w l dt 1 lc rgb 'red' lw 2,\
'off_twiss_-0.02.out'    u 2:(($20-$7*-0.02)*1e3) title '-2%' w l dt 2 lc rgb 'red' lw 2,\
'off_twiss_0.03.out'     u 2:(($20-$7*0.03)*1e3) title '3%' w l dt 1 lc rgb 'green' lw 2,\
'off_twiss_-0.03.out'    u 2:(($20-$7*-0.03)*1e3) title '-3%' w l dt 2 lc rgb 'green' lw 2,\
'off_twiss_0.0325.out'     u 2:(($20-$7*0.0325)*1e3) title '3.25%' w l dt 1 lc rgb 'brown' lw 2,\
'off_twiss_-0.0325.out'    u 2:(($20-$7*-0.0325)*1e3) title '-3.25%' w l dt 2 lc rgb 'brown' lw 2,\
'off_twiss_0.0.out' u 2:($20*1e3) title 'on-energy' w l lt -1 lw 2

set output 'off_optics_xco.eps'
set ylabel 'x_{co} (mm)'
plot 'off_twiss_0.01.out'     u 2:($20*1e3) title '1%' w l dt 1 lc rgb 'blue' lw 2,\
'off_twiss_-0.01.out'    u 2:($20*1e3) title '-1%' w l dt 2 lc rgb 'blue' lw 2,\
'off_twiss_0.02.out'     u 2:($20*1e3) title '2%' w l dt 1 lc rgb 'red' lw 2,\
'off_twiss_-0.02.out'    u 2:($20*1e3) title '-2%' w l dt 2 lc rgb 'red' lw 2,\
'off_twiss_0.03.out'     u 2:($20*1e3) title '3%' w l dt 1 lc rgb 'green' lw 2,\
'off_twiss_-0.03.out'    u 2:($20*1e3) title '-3%' w l dt 2 lc rgb 'green' lw 2,\
'off_twiss_0.0325.out'     u 2:($20*1e3) title '3.25%' w l dt 1 lc rgb 'brown' lw 2,\
'off_twiss_-0.0325.out'    u 2:($20*1e3) title '-3.25%' w l dt 2 lc rgb 'brown' lw 2,\
'off_twiss_0.0.out' u 2:($20*1e3) title 'on-energy' w l lt -1 lw 2
