head -1 -q seed_*/avg_beta_beats.dat > x_beats.dat
tail -q -n 1 seed_*/avg_beta_beats.dat > y_beats.dat
grep metric seed_*/00da/da.out | awk '{print $(NF-1) "   " $NF}' > metric.dat
paste x_beats.dat y_beats.dat metric.dat > beats_metric.dat
rm x_beats.dat y_beats.dat metric.dat
