set term postscript enhanced color solid font 18

set output 'ibs_x.eps'

set xlabel 'Current (mA)'
set ylabel 'Horizontal Emittance (nm)'

plot 'emittance.dat' u ($1*1000):($2*1e9) notitle w lp

set output 'ibs_y.eps'

set xlabel 'Current (mA)'
set ylabel 'Vertical Emittance (pm)'

plot 'emittance.dat' u ($1*1000):($3*1e12) notitle w lp

set output 'ibs_e.eps'

set xlabel 'Current (mA)'
set ylabel 'Energy Spread (10^{-4})'

plot 'emittance.dat' u ($1*1000):($4*1e4) notitle w lp
