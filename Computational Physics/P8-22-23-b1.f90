program Practica8
	implicit none
	! magnituds escalars
	double precision x,dx,x0,E,E1,E2,eps,hw,beta,sigma, prob
	!magnituds vectorials
	double precision ini(2),fin1(10000),fin2(10000),fin3(10000),energies(6),B(3)
	integer npas, i,j
	external dyy,RK41STEP,simpson,RungeKutta,V
	common/equacio/E,beta
	common/posi/dx,x0
	!definim algunes quantitats importants
	hw=2.76041d0
	energies=(/0.43d0*hw,0.48d0*hw,1.4d0*hw,1.45d0*hw,2.d0*hw,2.1d0*hw/)
	npas=500
	dx=16.d0/npas
	x0=-8.d0
	eps=1.d-8
	beta=0.d0

!	calculem les solucions per les energies de l'apartat 1, i les apuntam al document
	open(100,file="P8-dades1.dat")
	write(100,'(//A40)')"#E1:"
	E=energies(1)
	ini=(/0.d0,1.d-5/)
	call RungeKutta(x,dx,ini,fin1,2,npas)
	x=x0
	i=0
	do while (x.LE.4.d0)
		write(100,*) x,fin1(i)
		i=i+1
		x=x0+i*dx
	end do
	write(100,'(//A40)')"#E2:"
	E=energies(2)
	ini=(/0.d0,1.d-5/)
	call RungeKutta(x,dx,ini,fin1,2,npas)
	x=x0
	i=0
	do while (x.LE.4.d0)
		write(100,*) x,fin1(i)
		i=i+1
		x=x0+i*dx
	end do
	write(100,'(//A40)')"#E3:"
	E=energies(3)
	ini=(/0.d0,1.d-5/)
	call RungeKutta(x,dx,ini,fin1,2,npas)
	x=x0
	i=0
	do while (x.LE.4.d0)
		write(100,*) x,fin1(i)
		i=i+1
		x=x0+i*dx
	end do
	write(100,'(//A40)')"#E4:"
	E=energies(4)
	ini=(/0.d0,1.d-5/)
	call RungeKutta(x,dx,ini,fin1,2,npas)
	x=x0
	i=0
	do while (x.LE.4.d0)
		write(100,*) x,fin1(i)
		i=i+1
		x=x0+i*dx
	end do


!		Iniciam el metode de tir, obrim el fitxer per escriure la convergencia dels autovalors

	open(101,file="convergencia.dat")
	do j=0,4,2
		write(101,'(//A40)')"#convergencia:"
		E1=energies(j+1)
		E2=energies(j+2)
		E=E1
	    ini=(/0.d0,1.d-5/)
	    call RungeKutta(x,dx,ini,fin1,2,npas)
	    
	    E=E2
	    ini=(/0.d0,1.d-5/)
	    call RungeKutta(x,dx,ini,fin2,2,npas)

	    call metodetir(E1,E2,E,eps,fin1,fin2,fin3,npas)

	!   Normalitzem aquestes autofuncions
	    ini=(/0.d0,1.d-5/)
	    call  RungeKutta(x,dx,ini,fin3,2,npas)
	    call normalitzar(fin3,npas)

	!   Finalment les escrivim en un fitxer
	    write(100,'(//A40)')"#Autovector:"
	    do i=1,npas
          write(100,*)x0+i*dx, fin3(i)
        end do
    end do
    write(101,'(//A40)')"#basura:"


!		Per realitzar el darrer apartat definim el vector B(beta1,beta2,beta3)
!		Redefinim x0 i dx
	open(102,file="P8-22-23-b1-res1.dat")
	dx=10.d0/npas
	x0=-5.d0
	B=(/0.d0,0.1d0,0.25d0/)
	sigma=1.175d0
	do j=1,3  
		prob=0.d0
		beta=B(j)
		E1=energies(1)
		E2=energies(2)
		E=E1
	    ini=(/0.d0,1.d-5/)
	    call RungeKutta(x,dx,ini,fin1,2,npas)
	    
	    E=E2
	    ini=(/0.d0,1.d-5/)
	    call RungeKutta(x,dx,ini,fin2,2,npas)

	    call metodetir(E1,E2,E,eps,fin1,fin2,fin3,npas)

	!   Normalitzem aquestes autofuncions
	    ini=(/0.d0,1.d-5/)
	    call  RungeKutta(x,dx,ini,fin3,2,npas)
	    call normalitzar(fin3,npas)
!		Calculam la suma de la probabilitat
	    do i=1,npas
	    	if ((i.GT.191).AND.(i.LT.309)) prob=prob+fin3(i)
	    end do
	    write(100,'(//A40)')"#Probabilitat:"
	    write(102,*) prob
	!   Finalment les escrivim en un fitxer
	    write(100,'(//A40)')"#Autovector:"
	    do i=1,npas
          write(100,*)x0+i*dx, fin3(i)
        end do
    end do 


    close(100)
    close(101)
    close(102)
end program


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
        DOUBLE PRECISION YY(NEQUA),YYTEMP(NEQUA),YYFINAL(NEQUA),resul(npas)
        DOUBLE PRECISION X,DX,x0,basura
        DOUBLE PRECISION  K1(NEQUA), K2(NEQUA), K3(NEQUA), K4(NEQUA)
        common/posi/basura,x0
        X=x0
        do i=1,npas
            call RK41STEP(X,DX,YY,YYFINAL,2)
            X=x0+DX*i
            YY=YYFINAL
            resul(i)=YYFINAL(1)
        end do
        return
    end

!		Cream la subrutina per les derivades, on el potencial és una funció
SUBROUTINE DYY(x,yy,K,neq)
		implicit none
		double precision x,yy(neq),k(neq),V,E
		integer neq
		common/equacio/E
		k(1)=yy(2)
		k(2)=0.2624706361d0*(V(x)-E)*yy(1)
		return
	end

!	Cream la funció potencial, d'on podem modificar beta
function V(x)
	implicit none
	double precision x, beta,V,E 
	common/equacio/E,beta
	V=0.5d0*x**2+beta*x**4
	return
end

!       Subrutina integrals per simpson a partir d'una llista de valors
subroutine simpson(y,a,b,ndat,integral)
        implicit none
        double precision y(3,ndat),integral,h,int,a,b
        integer ndat,i
        h=(b-a)/ndat
        int=h*(y(2,1)+y(2,ndat))/3.d0

        do i=1,ndat-1
            if(mod(i,2).eq.0) then
                int=int+2*y(2,i)
            else
                int=int+4*y(2,i)
            endif
        end do
        integral=h/3.d0*int
    end 

!   Cream una subrutina a la que li dones una funcio i te la retorna normalitzada
subroutine normalitzar(func,npas)
    implicit none
    double precision func(npas),auto(3,npas),norm,dx,x0
    integer npas,i
    common/posi/dx,x0
    do i=1,npas
          auto(1,i)=x0+i*dx
          auto(2,i)=func(i)**2
          auto(3,i)=func(i)
    end do
    call simpson(auto,0.d0,1.d0,npas,norm)
    do i=1,npas
          func(i)=auto(3,i)/sqrt(norm)
    end do
    return
end


!   Programam una subrutina per fer el metode de tir
subroutine metodetir(E1,E2,E,eps,fin1,fin2,fin3,npas)
    implicit none
    double precision E1,E2,E,eps,x0
    double precision fin1(npas),fin2(npas),fin3(npas),ini(2),t,dx
    integer npas,i
    common/posi/dx,x0
    i=0
    fin3(npas)=1.d0
!   Comença el bucle que cerca els zeros: Metode de tir
    do while (abs(fin3(npas)).GT.eps)
      i=i+1
      ini=(/0.d0,0.25d0/)
      E=(E1*fin2(npas)-E2*fin1(npas))/(fin2(npas)-fin1(npas))
      
        call  RungeKutta(t,dx,ini,fin3,2,npas)
        write(101,*)i,E
!     Assignem les noves variables
      E1=E2
      E2=E
      fin1(npas)=fin2(npas)
      fin2(npas)=fin3(npas)
    end do
    return
end