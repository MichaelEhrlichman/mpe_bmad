set term postscript enhanced color solid font 28
set output 'ps_x.eps'

set title 'On-Energy'

set xlabel 'x (mm)'
set ylabel 'p_x'

plot '00phase_space/phase_space_0.0_x.out' u ($4*1000):($5*1000) notitle w p pt 7 ps 0.2 lc rgb 'black'

set output 'ps_y.eps'

set title 'On-Energy'

set xlabel 'y (mm)'
set ylabel 'p_y'

plot '00phase_space/phase_space_0.0_y.out' u ($6*1000):($7*1000) notitle w p pt 7 ps 0.2 lc rgb 'black'


 set output 'ps_0.03_x.eps'
 
 set title '+3% Energy Offset'
 
 set xlabel 'x (mm)'
 set ylabel 'p_x'
 
 plot '00phase_space/phase_space_0.03_x.out' u ($4*1000):($5*1000) notitle w p pt 7 ps 0.2 lc rgb 'black'
 
 set output 'ps_0.03_y.eps'
 
 set title '+3% Energy Offset'
 
 set xlabel 'y (mm)'
 set ylabel 'p_y'
 
 plot '00phase_space/phase_space_0.03_y.out' u ($6*1000):($7*1000) notitle w p pt 7 ps 0.2 lc rgb 'black'
 
 
 set output 'ps_-0.03_x.eps'
 
 set title '-3% Energy Offset'
 
 set xlabel 'x (mm)'
 set ylabel 'p_x'
 
 plot '00phase_space/phase_space_-0.03_x.out' u ($4*1000):($5*1000) notitle w p pt 7 ps 0.2 lc rgb 'black'
 
 set output 'ps_-0.03_y.eps'
 
 set title '-3% Energy Offset'
 
 set xlabel 'y (mm)'
 set ylabel 'p_y'
 
 plot '00phase_space/phase_space_-0.03_y.out' u ($6*1000):($7*1000) notitle w p pt 7 ps 0.2 lc rgb 'black'

