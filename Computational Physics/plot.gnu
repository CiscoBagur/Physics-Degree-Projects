file1="posicions.dat"


set xlabel "t "
set ylabel "y (adim)"

set title "trajectòria y(t)"
plot file1 index 0 u 1:3 w l t"y(t)"
set term png
set output "Exa-jan-23-fig1.png"
replot

set xlabel "x (adim) "
set ylabel "z (adim)"

set title "Espai fàsic"
plot file1 index 0 u 2:4 w l t"z(x)" 
set term png
set output "Exa-jan-23-fig2.png"
replot
