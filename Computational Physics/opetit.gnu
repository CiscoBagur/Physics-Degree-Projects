file1="P7-22-23-res1.dat"


set xlabel "t (s)"
set ylabel "angle (rad)"

set title "Oscilacions petites"
plot file1 index 0 u 1:2 w l t"Euler",file1 index 1 u 1:2 w l t"RK2"
set term png
set output "P7-22-23-fig1.png"
replot

set title "Oscilacions grans"
plot file1 index 2 u 1:2 w l t"Euler",file1 index 3 u 1:2 w l t"RK2"
set term png
set output "P7-22-23-fig2.png"
replot

set xlabel "angle (rad)"
set ylabel "angle punt (rad)"

set title "Espai fàsic"
plot file1 index 2 u 2:3 w l t"Euler",file1 index 3 u 2:3 w l t"RK2"
set term png
set output "P7-22-23-fig3.png"
replot

set title "Transició espai fàsic"
plot file1 index 6 u 2:3 w l t"suma",file1 index 7 u 2:3 w l t"resta"
set term png
set output "P7-22-23-fig5.png"
replot

set xlabel "t (s)"
set ylabel "Energia (J)"

set title "Energia cinètica vs total"
plot file1 index 4 u 1:2 w l t"EKINE Euler",file1 index 4 u 1:4 w l t"ETOT Euler",file1 index 5 u 1:2 w l t"EKINE RK2",file1 index 5 u 1:4 w l t"ETOT RK2"
set term png
set output "P7-22-23-fig4.png"
replot

set yrange[0:15]
set title "Energia total, convergència"
plot file1 index 8 u 1:2 w l t"N=300",file1 index 9 u 1:2 w l t"N=550",file1 index 10 u 1:2 w l t"N=1000", file1 index 11 u 1:2 w l t"N=20000"
set term png
set output "P7-22-23-fig6.png"
replot

#El primer punt de 20000 dona un valor loco i peta



set term gif size 1000,600 animate  delay 20 loop 0 
set output "anima2.gif"
datafile =file1


do for[i=1:1000:5]{

set multiplot

set size 0.5,0.8
set origin 0.0,0.0
set title "Como funcion del tiempo"
set size 0.5,0.8
set xrange[0:10]
set yrange[-2*pi:2*pi]
set xlabel "t (s)"
set ylabel "theta, thetap"
plot datafile index 3 every ::1::i with line linewidth 4 t"theta" ,datafile index 3 every ::1::i u 1:3 with line linewidth 4 t"thetap"

set origin 0.5,0
set size 0.5,0.8
set title "Diagrama de fases"
set yrange[-2*pi:2*pi]
set xrange[-pi:pi]
set xlabel "theta"
set ylabel "thetap"
if (i>200) { plot datafile index 3 every::i::i u 2:3 t"" ps 3,datafile index 3 every::i-100::i u 2:3 w l t"" } else {
plot datafile index 3 every::i::i u 2:3 t"" ps 3}
unset multiplot

}