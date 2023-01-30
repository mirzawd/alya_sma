# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>

#set terminal aqua dashed enhanced
#set terminal x11 dashed nopersist enhanced font "arial,15"
#set terminal wxt dashed nopersist enhanced font "arial,15"

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "One element (1-)"
set xlabel "Strain (-)"
set ylabel "Stress (MPa)"
set xr [0:1.0]
set yr [-100:2000]
set key right top
set grid

system("alya-sets sm154_d1C-boundary.sld.set 1")
   
#
# define line styles using explicit rgbcolor names
#
set for [i=1:3] linetype i dashtype i
set style line 1 lt 2 lc rgb "blue" lw 1.5 pt 1
set style line 2 lt 1 lc rgb "black" lw 1.5 pt 2
set style line 3 lt 1 lc rgb "red" lw 1.5 pt 2

plot "exp-XC.txt" using 1:2 title 'XC' with lines ls 1, \
     "the-XC.txt" using 1:2 title 'Theory (fXC = 0.1 & fGC = 0.85)' with lines ls 2, \
     "alya-sets-mu045-inf.txt" using (-$3/0.194):(-$4/0.194/0.194) title 'v23=0.45, INF' with lines, \
     "alya-sets-mu000-inf.txt" using (-$3/0.194):(-$4/0.194/0.194) title 'v23=0.0, INF' with lines, \
     "alya-sets-mu045-gre.txt" using (-$3/0.194):(-$4/0.194/0.194) title 'v23=0.45, GRE' with lines, \
     "alya-sets-mu000-gre.txt" using (-$3/0.194):(-$4/0.194/0.194) title 'v23=0.0, GRE' with lines, \
     "alya-sets.out" using (-$3/0.194):(-$4/0.194/0.194) title 'Curr' with lines ls 3



