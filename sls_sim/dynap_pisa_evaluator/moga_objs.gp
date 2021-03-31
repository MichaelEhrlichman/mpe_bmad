set terminal postscript eps color enhanced size 15,10 font "Arial" 32
set output 'report_objs.eps'

datafile='../moga_results.out'

countCommand='head -3 '.datafile.' | tail -1 | wc -w'
ncols = system(countCommand)
stats datafile

ngen = int(STATS_blocks-1)
stepsize = ngen/13
npop = int(STATS_records/(STATS_blocks-1))

fG = ngen%stepsize
fG = (fG==0)?stepsize:fG
nGenToPlot = ceil(1.0*ngen/stepsize)
o1 = int(ncols-3)
o2 = int(ncols-2)
o3 = int(ncols-1)
fe = int(ncols)
print ngen,o1,o2,o3,fe

set style line 20 linetype 7 linecolor rgb "black" ps 0.5
set style line 19 linetype 1 linecolor rgb "blue"
set style line 18 linetype 2 linecolor rgb "red"
set style line 17 linetype 4 linecolor rgb "green"
set style line 16 linetype 1 linecolor rgb "magenta"
set style line 15 linetype 2 linecolor rgb "orange"
set style line 14 linetype 4 linecolor rgb "grey"
set style line 13 linetype 1 linecolor rgb "blue"
set style line 12 linetype 2 linecolor rgb "red"
set style line 11 linetype 4 linecolor rgb "green"
set style line 10 linetype 1 linecolor rgb "magenta"
set style line  9 linetype 2 linecolor rgb "orange"
set style line  8 linetype 4 linecolor rgb "grey"

unset key
set size 1,1
#set multiplot layout 2,3 title 'db00g MOGA'
set multiplot layout 2,3

set log xy
set xlabel 'DA at 0%'
set ylabel 'DA at -3%'
set origin 0,0.5
set size 0.333,0.5
plot for [IDX=fG:ngen:stepsize] datafile i IDX-1 u (column(fe)==1?column(o1):1/0):(column(fe)==1?column(o2):1/0) w p ls 21-nGenToPlot+(IDX-fG)/stepsize

set xlabel 'DA at 0%'
set ylabel 'DA at 3%'
set origin 0.333,0.5
set size 0.333,0.5
plot for [IDX=fG:ngen:stepsize] datafile i IDX-1 u (column(fe)==1?column(o1):1/0):(column(fe)==1?column(o3):1/0) w p ls 21-nGenToPlot+(IDX-fG)/stepsize

set xlabel 'DA at -3%'
set ylabel 'DA at 3%'
set origin 0.666,0.5
set size 0.333,0.5
plot for [IDX=fG:ngen:stepsize] datafile i IDX-1 u (column(fe)==1?column(o2):1/0):(column(fe)==1?column(o3):1/0) w p ls 21-nGenToPlot+(IDX-fG)/stepsize
unset log xy

sedCommand(n) = sprintf("sed -n '/# Generation   %4d/,/^ *$/p' '../moga_results.out' | grep -v Generation > /dev/shm/objs.block.TEMP",n)
awkCommand(n1,n2,n3) = sprintf("cat /dev/shm/objs.block.TEMP | awk '{print $%2d, $%2d, $%2d}' > /dev/shm/objs.trans.TEMP",n1,n2,n3)

system(sedCommand(ngen))
system(awkCommand(o1,o2,o3))
!python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < /dev/shm/objs.trans.TEMP > /dev/shm/objs.ready.TEMP

unset xlabel
set ylabel 'DA metric'
set xrange[-0.25:2.25]
set xtics ('DA(0%%)' 0)
set xtics add ('DA(-3%%)' 1)
set xtics add ('DA(3%%)' 2)
set origin 0,0.
set size 0.666,0.5
plot for [COL=1:npop] '/dev/shm/objs.ready.TEMP' u COL w lp lw 4

unset label
set key
set key center
unset tics
unset border
unset xlabel
unset ylabel
set origin 0.666,0
set size 0.333,0.5
plot [0:1] [0:1] for [IDX=fG:ngen:stepsize] NaN t "generation ".IDX w p ls 21-nGenToPlot+(IDX-fG)/stepsize



