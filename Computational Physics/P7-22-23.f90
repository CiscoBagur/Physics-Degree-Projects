!By Francesc Bagur
Program P7
implicit none
double precision angle1(20000),apunt1(20000),temps1(20000)
double precision g,l,pi,m,omega,Tn,a0,ap0 !a --> angle, p ---> punt(derivada)
double precision EKIN,EPOT
integer i, k
 common/pendol/l,g,m

pi=4*atan(1.d0)
l=1.07d0
g=10.44d0
m=0.98d0
omega=sqrt(g/l)
Tn=2*pi/omega
a0=0.025d0
ap0=0.d0

open(100,file="P7-22-23-res1.dat")
 call EulerPendol(1500,0.d0,7*Tn,angle1,apunt1,temps1,g,l,a0,ap0)
do i=1,1500
	write(100,*)temps1(i),angle1(i),apunt1(i)
end do
write(100,'(//A40)')"#metode RK2 Index1"
 call RK2(1500,0.d0,7*Tn,angle1,apunt1,temps1,g,l,a0,ap0)			
 do i=1,1500
	write(100,*)temps1(i),angle1(i),apunt1(i)
 end do

call EulerPendol(1500,0.d0,7*Tn,angle1,apunt1,temps1,g,l,pi-0.15d0,ap0)
write(100,'(//A40)')"#Grans oscilacions Euler Index2"
do i=1,1500
	write(100,*)temps1(i),angle1(i),apunt1(i)
end do 
write(100,'(//A40)')"#grans oscilacions RK2 Index3"
 call RK2(1500,0.d0,7*Tn,angle1,apunt1,temps1,g,l,pi-0.15d0,ap0)
do i=1,1500
	write(100,*)temps1(i),angle1(i),apunt1(i)
end do

write(100,'(//A40)')"#Energies Euler Index4"
call EulerPendol(1500,0.d0,7*Tn,angle1,apunt1,temps1,g,l,pi-0.025d0,0.12d0)
do i=1,1500
	write(100,*) temps1(i),EKIN(apunt1(i)),EPOT(angle1(i)),EKIN(apunt1(i))+EPOT(angle1(i))
end do
write(100,'(//A40)')"#Energies RK2 Index5"
call RK2(1500,0.d0,7*Tn,angle1,apunt1,temps1,g,l,pi-0.025d0,0.12d0)
do i=1,1500
	write(100,*) temps1(i),EKIN(apunt1(i)),EPOT(angle1(i)),EKIN(apunt1(i))+EPOT(angle1(i))
end do

write(100,'(//A40)')"#Transicio suma Index6"
!EXPLICACIO TRANSICIO
!la diferencia que es veu en els dos grafics es deu a que quan la velocitat inicial supera el 
!valor de 2*omega, el pendol comença a donar voltes sobre ell mateix, fent que l'angle augmenti de 
!forma indefinida

call RK2(6000,0.d0,15*Tn,angle1,apunt1,temps1,g,l,0.d0,2*omega+0.05d0)
do i=1,6000
	write(100,*)temps1(i),angle1(i),apunt1(i)
end do
write(100,'(//A40)')"#Transicio resta Index7"
call RK2(6000,0.d0,15*Tn,angle1,apunt1,temps1,g,l,0.d0,2*omega-0.05d0)
do i=1,6000
	write(100,*)temps1(i),angle1(i),apunt1(i)
end do

write(100,'(//A40)')"#Convergencia Index8"
call RK2(300,0.d0,11*Tn,angle1,apunt1,temps1,g,l,2.87d0,0.d0)
do i=1,300
	write(100,*)temps1(i),EKIN(apunt1(i))+EPOT(angle1(i))
end do
write(100,'(//A40)')"#Mes iteracions! Index9"
call RK2(550,0.d0,11*Tn,angle1,apunt1,temps1,g,l,2.87d0,0.d0)
do i=1,550
	write(100,*)temps1(i),EKIN(apunt1(i))+EPOT(angle1(i))
end do
write(100,'(//A40)')"#Mes iteracions! Index10"
call RK2(1000,0.d0,11*Tn,angle1,apunt1,temps1,g,l,2.87d0,0.d0)
do i=1,1000
	write(100,*)temps1(i),EKIN(apunt1(i))+EPOT(angle1(i))
end do
write(100,'(//A40)')"#Mes iteracions! Index11"
call RK2(20000,0.d0,11*Tn,angle1,apunt1,temps1,g,l,2.87d0,0.d0)
do i=1,20000
	write(100,*)temps1(i),EKIN(apunt1(i))+EPOT(angle1(i))
end do

close(100)
end program


subroutine RK2(N,a,b,th,thp,t,g,l,th0,thp0)
implicit none
!variables d'entrada
double precision a,b,th0,thp0,g,l
integer N  
!variables de sortida
double precision th(N),thp(N),t(N)
!varibles mudes
double precision dt
integer i
dt=(b-a)*1.d0/N
th(1)=th0
thp(1)=thp0
t(1)=a

do i=1,N 
	th(i+1)=th(i)+dt*(thp(i)-dt*g*sin(th(i))/(2*l))
	thp(i+1)=thp(i)-dt*g*sin(th(i)+dt*thp(i)/2.d0)/l
	t(i+1)=a+i*dt
end do 

return
end 



subroutine EulerPendol(N,a,b,theta,thetapunt,t,g,l,thzero,thpzero)
implicit none
double precision a,b,theta(N),thetapunt(N),g,l,thzero,thpzero,dt,t(N),t0
integer N,i
t0=a
dt=(b-a)*1.d0/N
theta(1)=thzero
thetapunt(1)=thpzero
t(1)=t0
!write(*,*)dt,a,b,N
do i=1,N-1 
	theta(i+1)=theta(i)+dt*thetapunt(i)
	thetapunt(i+1)=thetapunt(i)-dt*(g/l)*sin(theta(i))
	t(i+1)=t(1)+i*dt
end do

return
end


function EKIN(thp)
implicit none
double precision EKIN,thp,l,g,m
	common/pendol/l,g,m
	EKIN=m*thp**2*l**2/2.d0 
return
end  

function EPOT(th)
implicit none
double precision EPOT,th,l,g,m
	common/pendol/l,g,m
	EPOT=-m*g*l*cos(th)
return
end
