set terminal postscript eps color noenhanced size 15,30 font "Arial" 18
set output 'report_vars.eps'

datafile = '../moga_results.out'

fgen = 1

countCommand='head -3 '.datafile.' | tail -1 | wc -w'
ncols = system(countCommand)
stats datafile
ngen = int(STATS_blocks-1)

header = system('head -1 '.datafile)

unset key
set multiplot layout 12,3 title 'MOGA Variable Report'

set xlabel 'generation'

set auto xy
set xtics scale 0.1
set ytics scale 0.1
set xrange[fgen-1:ngen+1]
set grid

do for [jcol=2:(ncols-4)] {
  set ylabel word(header, (jcol+1))
  plot for [IDX=fgen-1:ngen-1] datafile i IDX u (IDX+1):jcol w d lt rgb "red"
}

unset multiplot
set output 'report_recent_vars.eps'
set multiplot layout 12,3 title 'MOGA Recent Variable Report'

fgen=(ngen-50>1?ngen-50:1)
set xrange[fgen-1:ngen+1]

do for [jcol=2:(ncols-4)] {
  set ylabel word(header, (jcol+1))
  plot for [IDX=fgen-1:ngen-1] datafile i IDX u (IDX+1):jcol w d lt rgb "red"
}
