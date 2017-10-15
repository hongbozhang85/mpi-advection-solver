set terminal eps truecolor font ",15"
set output "t3_exp1.eps"
set multiplot
set key top left
# set key sample 1
# set size 0.6,0.7
set xtics("1/64"0.015625, "1/16"0.0625, "1/8"0.125,"1/4"0.25,"1/2"0.5)
set ytics 0,0.1,0.5
set xlabel '1/Q'
set ylabel 'time'
y(x)=a*x+b
fit y(x) "t3_exp1_plot.dat" using (1/$1):2 via a,b
plot "t3_exp1_plot.dat" using (1/$1):2 with p ps 1.0 pt 3 lw 2 lt 1 t 'Exp', \
 y(x) with l lw 2 lt 2 t 'Fit' 
