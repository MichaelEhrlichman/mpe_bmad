set term postscript enhanced color font 21

set output 'beta_beating_hulls.eps'

set xlabel 'Average Percent {/Symbol b}_x Beat'
set ylabel '1/2 Raw DA Area (mm^2)'

plot 'ideal/stats.dat' u 1:($2*1e6) t 'Ideal Lattice' w p pt 7 ps 2 lc rgb 'black',\
'rms_05/stats.dat' u 1:($3*1e6) t 'K_1 %RMSE = 0.05' w p lc rgb 'green' pt 7 ps 1,\
'rms_05/50_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'green' lw 3 lt 1,\
'rms_05/90_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'green' lw 3 lt 3,\
'rms_10/stats.dat' u 1:($3*1e6) t 'K_1 %RMSE = 0.10' w p lc rgb 'blue' pt 7 ps 1,\
'rms_10/50_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'blue' lw 3 lt 1,\
'rms_10/90_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'blue' lw 3 lt 3,\
'rms_15/stats.dat' u 1:($3*1e6) t 'K_1 %RMSE = 0.15' w p lc rgb 'red' pt 7 ps 1,\
'rms_15/50_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'red' lw 3 lt 1,\
'rms_15/90_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'red' lw 3 lt 3,\
'rms_20/stats.dat' u 1:($3*1e6) t 'K_1 %RMSE = 0.20' w p lc rgb 'magenta' pt 7 ps 1,\
'rms_20/50_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'magenta' lw 3 lt 1,\
'rms_20/90_pct_hull.dat' u 1:($2*1e6) notitle w l lc rgb 'magenta' lw 3 lt 3
