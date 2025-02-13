file2="P6-22-23-res.dat"

set xlabel "N"
set ylabel "I1"

I1=pi**2*(2*pi**2+3./2.)

plot file2 index 0 u 1:2:3 w yerrorbars t"Error real", I1 t"valor teoric"

pause -1
set term png
set output "P6-22-23-c2-fig1.png"
replot

set ylabel "I2"

plot file2 index 1 u 1:2:3 w yerrorbars t"Error real"


pause -1
set term png
set output "P6-22-23-c2-fig2.png"
replot