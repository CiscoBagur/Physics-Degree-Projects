file1="P9-res2.dat"
file2="P9-dat.dat"

set xlabel "PAS"
set ylabel "emax (ºC)"

set key b r
set title "Convergencia T0: 15ºC"
plot file2 index 0 u 1:2 t"Jacobi",file2 index 1 u 1:2 t"Sobrerelaxació"
set term png
set output "P9-22-23-fig1.png"
replot

set key t r
set title "Convergencia T0: 220ºC"
plot file2 index 2 u 1:2 t"Jacobi",file2 index 3 u 1:2 t"Sobrerelaxació"
set term png
set output "P9-22-23-fig2.png"
replot

set title "Convergencia T0: 1280ºC"
plot file2 index 4 u 1:2 t"Jacobi",file2 index 5 u 1:2 t"Sobrerelaxació"
set term png
set output "P9-22-23-fig3.png"
replot

set ticslevel 0

set xlabel "x (cm)"
set ylabel "y (cm)"

set xrange [0:45.5]
set yrange [0:33.5]

set view 0,0

set colorbox vertical user origin .05, .2 size .04,.65
set title "Poisson 2D, Sobrerelaxació"
splot file1 index 0 u 1:2:3 w pm3d t"T (ºC)"
set term png
set output "P9-22-23-fig4.png"
replot

set title "Poisson 2D (no font), Sobrerelaxació"
splot file1 index 1 u 1:2:3 w pm3d t"T (ºC)"
set term png
set output "P9-22-23-fig5.png"
replot
