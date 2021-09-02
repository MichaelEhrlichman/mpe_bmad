mu = 0.02
wi = 0.001
f(x) = 3500.0 / (exp((x-mu)/wi)+1.0)

set xlabel 'delta_p'
set ylabel 'turns survived'

plot 'along_0p.dat' u 1:($2>0?$2:3500) w p lt 7,\
f(x)


pause -1

