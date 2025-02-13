program Final2023
implicit none

!	Variables problema 1
double precision sigma,ro,beta,t,dt,ini(3),fin1(3,10000)
integer i, passes

!	Variables problema 2
double precision integ,error,dades(1000000)
integer k,N
external DYY, RK41STEP,fun1,fun2,fun3
	common/lorentz/sigma,beta,ro

	sigma=10.d0
	ro=28.d0
	beta=8.d0/3.d0

!	Problema 1, apartat a)
!	Elegim un pas de temps

	passes=2000
	dt=20.d0/(passes*1.d0)

!	Establim les condicions inicials
	ini=(/0.d0,1.d0,0.d0/)
	open(100,file="Exa-jan-23-res1.dat")
	open(101,file="posicions.dat")
	write(100,'(//A40)')"#Resolucio equacions:"
	call RungeKutta(t,dt,ini,fin1,3,passes)
	write(101,'(//A40)')"#Resolucio equacions:"
	do i=1,passes+1
		write(101,*)(i-1)*dt,fin1(1,i),fin1(2,i),fin1(3,i)
	end do


!	Problema 2, apartat a)
	k=220 								!Nombre d'intervals necessari per tenir precisio vuit decimals
										!Calculat manualment, ja que l'error creix com h**4 i h=(b-a)/k
	call simpson(0.8d0,3.d0,k,fun1,integ)
	write(100,'(//A40)')"#Integral metode simpson:"
	write(100,*)integ

!	Problema 2, apartat b)
	N=1000000
	call acceptrebuig(N,dades,0.8d0,3.d0,1.2d0,fun2)

	
	call monteimport(N,fun3,dades,integ,error)
	write(100,'(//A40)')"#Integral Monte-Carlo amb sampleig:"
	write(100,*)integ,"+-",error


end program

!	Resoldrem el sistema no lineal amb el mètode de rungekutta

!	Copiam totes les subrutines relatives al rungekutta
!	Les modifiquem per resoldre el sistema d'equacions plantejat

SUBROUTINE RK41STEP(X,DX,YY,YYFINAL,NEQUA)
        IMPLICIT NONE
! NUMERO DE ECUACIONES
        INTEGER NEQUA,I,J
! VALOR PREVIO Y POSTERIOR
        DOUBLE PRECISION YY(NEQUA),YYTEMP(NEQUA),YYFINAL(NEQUA)
        DOUBLE PRECISION X,DX
        DOUBLE PRECISION  K1(NEQUA), K2(NEQUA), K3(NEQUA), K4(NEQUA)
!       EXTERNAL DYY

          CALL DYY(X,YY,K1,NEQUA)
          DO J=1,NEQUA
            YYTEMP(J)=YY(J)+DX*K1(J)/2.D0
          ENDDO

          CALL DYY(X+DX/2.D0,YYTEMP,K2,NEQUA)
          DO J=1,NEQUA
            YYTEMP(J)=YY(J)+DX*K2(J)/2.D0
          ENDDO

          CALL DYY(X+DX/2.D0,YYTEMP,K3,NEQUA)
          DO J=1,NEQUA
            YYTEMP(J)=YY(J)+DX*K3(J)
          ENDDO

          CALL DYY(X+DX,YYTEMP,K4,NEQUA)
          DO I=1,NEQUA
          YYFINAL(I)=YY(I)+1.D0/6.D0*(K1(I)+2.D0*K2(I)+2.D0*K3(I)+K4(I))*DX
          ENDDO

         END


!       Subrutina bucle per fer Runge Kutta
subroutine RungeKutta(X,DX,YY,resul,NEQUA,npas)
        IMPLICIT NONE
! NUMERO DE ECUACIONES
        INTEGER NEQUA,I,J,npas
! VALOR PREVIO Y POSTERIOR
        DOUBLE PRECISION YY(NEQUA),YYTEMP(NEQUA),YYFINAL(NEQUA),resul(NEQUA,npas+1)
        DOUBLE PRECISION X,DX
        DOUBLE PRECISION  K1(NEQUA), K2(NEQUA), K3(NEQUA), K4(NEQUA)
        X=0
        resul(1,1)=YY(1)
        resul(2,1)=YY(2)
        resul(3,1)=YY(3)
        do i=1,npas
            call RK41STEP(X,DX,YY,YYFINAL,NEQUA)
            X=DX*i
            YY=YYFINAL
            resul(1,i+1)=YYFINAL(1)
            resul(2,i+1)=YYFINAL(2)
            resul(3,i+1)=YYFINAL(3)
            if (mod(X-2.d0,4.d0).Eq.0.d0) write(100,*)x,YYFINAL(1),YYFINAL(2),YYFINAL(3)
        end do
        return
    end

!		Cream la subrutina per les derivades, on el potencial és una funció
SUBROUTINE DYY(x,yy,K,neq)
implicit none
double precision x,yy(neq),k(neq),sigma, beta, ro
integer neq
		common/lorentz/sigma,beta,ro
		k(1)=sigma*(yy(2)-yy(1))
		k(2)=yy(1)*(ro-yy(3))-yy(2)
		k(3)=yy(1)*yy(2)-beta*yy(3)
		return
end

!	Funcio a estudiar en el problema 2
function fun1(x)
implicit none
double precision x, fun1
	fun1=1.d0/(1+sinh(2.d0*x)*log(x)**2)
	return
end

!	Funcio distribucio problema 2
function fun2(x)
implicit none
double precision x, fun2
	fun2=exp(-x)/(exp(-0.8d0)-exp(-3.d0))
	return
end

!	Funcio per fer montecarlo amb sampleig d'importancia
function fun3(x)
implicit none
double precision x, fun1, fun2, fun3
	fun3=fun1(x)/fun2(x)
	return
end

! 		Subrutina que calcula integrals pel mètode simpson
subroutine simpson(x1,x2,N,funci,integral)
implicit none
double precision x1,x2,integral, xi,f1,f2,f3,h,simp_i,funci
integer N, part, i

	integral=0
	part=N
	h=(x2-x1)/part
	do i=0,part-2,2
		xi=x1+i*h
		f1=funci(xi)
		f2=funci(xi+h)
		f3=funci(xi+2*h)
		simp_i=h*(f1+4*f2+f3)/3.d0
		integral=integral+simp_i
	end do

	return
end


! subrutina que genera una llista de nombres aleatoris distribuits segons g(x)
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

!	Subrutina que fa integrals per montecarlo amb sampleig d'importància
subroutine monteimport(N,funci,distr,vint,error)
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
