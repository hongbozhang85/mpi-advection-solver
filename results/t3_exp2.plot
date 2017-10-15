set terminal eps truecolor font ",15"
#set output "t3_exp2.eps"
set output "t3_exp2_2.eps"
set multiplot
set key top left
# set key sample 1
# set size 0.6,0.7
#set xtics("100"100,"10^5"100000,"10^6"1000000)
#set ytics 0,10,50
set xlabel 'M'
set ylabel 'time'
y(x)=a*x+b
#fit y(x) "t3_exp2_new_plot.dat" using 1:2 via a,b
fit y(x) "t3_exp2_2.dat" using 1:2 via a,b
#plot "t3_exp2_new_plot.dat" using 1:2 with p ps 1.0 pt 3 lw 2 lt 1 t 'Exp', \
plot "t3_exp2_2.dat" using 1:2 with p ps 1.0 pt 3 lw 2 lt 1 t 'Exp', \
 y(x) with l lw 2 lt 2 t 'Fit' 
