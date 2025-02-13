program practica6
implicit none

!variables primera integral
double precision pi,I1,sigmaI1
integer ndat,i,k
!variables per generar la primera distribucio
double precision rnums(1000000),L
!variables segona integral
double precision I2,sigmaI2,I3,sigmaI3

external prob, func1,func2
common/const/pi ,L   	


pi=4*atan(1.d0)
ndat=1000000
L=50.d0

open(100,file="P6-22-23-res.dat")
write(100,*)"#Per columnes: N, I1, sigmaI1"

do i=10000,ndat,10000
call MonteCru(0,2*pi,i,func1,I1,sigmaI1)
write(100,*)i,I1,sigmaI1
end do

call acceptrebuig(ndat,rnums,0,2*L,0.02,prob)
write(100, '(////A45)') "#Segona integral: N, I2, sigmaI2"
do i=10000,ndat,10000
call MonteImport(i,func2,rnums,I2,sigmaI2)
write(100,*)i,I2,sigmaI2
end do

write(100, '(////A45)') "#Tercera integral: N, I3, sigmaI3" !quant afegeixo aquesta subrutina el programa peta 
do i=10000,300000,1000000									!i deixa d'escriure en el fitxer, per aixo hi ha una copia
call ErradaEncert(0,2*pi,1/L**2,i,I3,sigmaI3)
write(100,*)i,I3,sigmaI3
end do
close(100)

end program

subroutine MonteCru(a,b,N,funci,vint,error) 			!a,b -> llimits integral esquerra i dreta. N -> nombre de punts
implicit none											!funci -> funcio a integrar. vint -> valor integral
double precision a, b, funci, vint,t,ht,error, Nr		!error -> estimacio error comes
integer i, N
vint=0
error=0
do i=1,N
call random_number(t)
ht=(b-a)*funci((b-a)*t+a)
vint=vint+ht
error=error+ht**2
end do
vint=vint/N
Nr=N*1.d0
error=sqrt(error/N-(vint)**2)/sqrt(Nr)
return
end

subroutine ErradaEncert(a,b,M,N,funci,vint,error) !la subrutina necessitara que aportis una cota maxima
implicit none
double precision a,b,funci,vint,error,rand1,rand2,x,p,conta,M
integer N,i
conta=0
do i=1,N 
call random_number(rand1)
call random_number(rand2)
x=a+rand1*(b-a)  			
p=rand2*M  					
if (funci(x).LE.p) then
conta=conta+1.d0
else
end if
end do 
vint=M*(b-a)*conta/N
error=M*(b-a)/sqrt(1.d0*N)*sqrt(conta/N*(1-conta/N))
return
end

subroutine acceptrebuig(ndat,xnums,a,b,M,fun)
implicit none
double precision xnums(ndat),a,b,M, x,p,aleatori1, aleatori2,fun
integer ndat,i, niter
i=0
niter=0
do while (i.LT.ndat)
niter=niter+1
call random_number(aleatori1)
call random_number(aleatori2)
x=a+aleatori1*(b-a)  			!escollim la base
p=aleatori2*M  					!escollim l'"altura"
if (fun(x).GE.p) then
i=i+1
xnums(i)=x
else
end if
end do
return
end

subroutine MonteImport(N,funci,distr,vint,error)
implicit none
double precision funci,distr(N),error,vint,xk,Nr
integer i ,N
vint=0
error=0
do i=1,N
xk=distr(i)
vint=vint+funci(xk)
error=error+(funci(xk))**2
end do
Nr=1.d0*N
vint=vint/N
error=sqrt(error/N-vint**2)/sqrt(Nr)

return
end

function func1(x)
implicit none
double precision x,func1
func1=x**3*cos(x)**2
return
end

function func2(x)
implicit none
double precision func2,x,L,pi 
common/const/pi,L 
func2=sin(pi*(x-2*L)/L)**2/L
return
end

function prob(x) !ja esta normalitzada
implicit none
double precision prob, x,L,pi
common/const/pi,L
prob=sin(pi*(x-2*L)/(2*L))**2/L
return
end
