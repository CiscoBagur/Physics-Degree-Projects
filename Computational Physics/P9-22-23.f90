program Practica9
	implicit none
	double precision t0x,t0y,tLy,tLx,Lx,Ly,h
	double precision T(0:182,0:134),Ti(0:182,0:134),eps,x,y
	integer i,j,Nx,Ny
	external rho,rho1,rho2,rho3,nofont,poisson2dpas
	common/inicials/t0y,t0x,tLy,tLx
	!dimensions de la caixa
	Lx=45.5d0
	Ly=33.5d0
	h=0.25d0
	Nx=int(Lx/h)
	Ny=int(Ly/h)

	!valors de les condicions de contorn
	t0y=0.5d0
	t0x=17.d0
	tLy=25.3d0
	tLx=11.2d0

!	Definim la precissio
	eps=1.d-5

!	inicialitzam els valors de Tij amb valors arbitraris
	do i=1,Nx-1
		do j=1,Ny-1
			T(i,j)=15.d0
		end do
	end do 
	!donam valors a la frontera
	do i=0,Nx
		T(i,0)=t0x
		T(i,Ny)=tLy
	end do
	do i=0,Ny
		T(0,i)=t0y
		T(Nx,i)=tLx
	end do 
	Ti=T
	open(101,file="P9-dat.dat")
!	Cridam la subbrutina i escrivim els valors de T(x,y) en un fitxer
	write(101,'(//A40)')"#Poisson 2D convergencia metode Jacobi:"
	call poisson2d(Nx,Ny,h,T,rho,0,eps)
	
	T=Ti
	write(101,'(//A40)')"#Poisson 2D convergencia metode Sobrerelaxacio:"
	call poisson2d(Nx,Ny,h,T,rho,1,eps)
	
	T=Ti
	open(100,file="P9-res2.dat")
	call poisson2dmap(Nx,Ny,h,T,rho,1,eps)
	write(100,'(//A40)')"#Poisson 2D metode Sobrerelaxacio:"
	do i=0,Nx
		x=i*h
		write(100,*)
		do j=0,Ny
		y=j*h
			write(100,*)x,y,T(i,j)
		end do 
	end do 

	call poisson2dmap(Nx,Ny,h,T,nofont,1,eps)
	write(100,'(//A40)')"#Poisson 2D metode Sobrerelaxacio (nofont):"
	do i=0,Nx
		x=i*h
		write(100,*)
		do j=0,Ny
		y=j*h
			write(100,*)x,y,T(i,j)
		end do 
	end do 
!	inicialitzam els valors de Tij amb valors diferents per l'apartat 3
	do i=1,Nx-1
		do j=1,Ny-1
			T(i,j)=220.d0
		end do
	end do 
	Ti=T
!	Cridam la subbrutina i escrivim els valors de T(x,y) en un fitxer
	write(101,'(//A40)')"#Poisson 2D convergencia metode Jacobi:"
	call poisson2d(Nx,Ny,h,T,rho,0,eps)
	
	T=Ti
	write(101,'(//A40)')"#Poisson 2D convergencia metode Sobrerelaxacio:"
	call poisson2d(Nx,Ny,h,T,rho,1,eps)
!	inicialitzam els valors de Tij amb valors diferents per l'apartat 3
	do i=1,Nx-1
		do j=1,Ny-1
			T(i,j)=1280.d0
		end do
	end do 
	Ti=T
!	Cridam la subbrutina i escrivim els valors de T(x,y) en un fitxer
	write(101,'(//A40)')"#Poisson 2D convergencia metode Jacobi:"
	call poisson2d(Nx,Ny,h,T,rho,0,eps)
	
	T=Ti
	write(101,'(//A40)')"#Poisson 2D convergencia metode Sobrerelaxacio:"
	call poisson2d(Nx,Ny,h,T,rho,1,eps)


end program


subroutine poisson2dpas(Nx,Ny,h,T,rho,ICONTROL,emax)
	implicit none
	double precision h,T(0:Nx,0:Ny),T_2(0:Nx,0:Ny),eps,y,x,rho
	double precision error,emax,t0x,t0y,tLx,tLy
	integer ICONTROL,i,j,Nx,Ny,PAS
	common/inicials/t0y,t0x,tLy,tLx
	T_2=T
!	Variable de control 0--> Jacobi
	if (ICONTROL.EQ.0) then
		do i=1,Nx-1 
			do j=1,Ny-1
				T(i,j)=(T_2(i+1,j)+T_2(i-1,j)+T_2(i,j+1)+T_2(i,j-1)+rho(i*h,j*h)*h**2)/4.d0
				!calculam l'error comes
				error=abs(T(i,j)-T_2(i,j))
				!write(*,*)error,emax
!				si aquest error es major que l'anterior error màxim que teniem guardat es substitueix el valor
				if (error.GT.emax) emax=error
			end do 
		end do
!	Variable de control 1--> Sobre-relaxacio
	else if (ICONTROL.EQ.1) then
		do i=1,Nx-1 
				do j=1,Ny-1
					T(i,j)=T(i,j)+1.45d0*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)-4.d0*T(i,j)+rho(i*h,j*h)*h**2)/4.d0
					!calculam l'error comes
					error=abs(T(i,j)-T_2(i,j))
	!				si aquest error es major que l'anterior error màxim que teniem guardat es substitueix el valor
					if (error.GT.emax) emax=error
				end do 
			end do
	else
	write(*,*)"l'has liada a s'ICONTROL xaval"
	stop
	end if
	end

subroutine poisson2d(Nx,Ny,h,T,rho,ICONTROL,eps)
	implicit none
	double precision h,T(0:Nx,0:Ny),T_2(0:Nx,0:Ny),eps,y,x,rho
	double precision error,emax,t0x,t0y,tLx,tLy
	integer ICONTROL,i,j,Nx,Ny,PAS
	common/inicials/t0y,t0x,tLy,tLx
	! Iniciam un contador de passos
	PAS=1
!	Iniciam bucle del metode,aplicam la condicio de convergencia
	emax=1.d0
	do while (abs(emax).GT.eps)
		write(101,*) PAS,T(102,54)
		emax=0.d0
		call poisson2dpas(Nx,Ny,h,T,rho,ICONTROL,emax)
		if (PAS.GT.60000) then
			write(*,*) "no ha convergit"
			exit
		else
			PAS=PAS+1
		end if
	end do
	write(*,*)"ha convergit amb pas:",PAS
	end

subroutine poisson2dmap(Nx,Ny,h,T,rho,ICONTROL,eps)
	implicit none
	double precision h,T(0:Nx,0:Ny),T_2(0:Nx,0:Ny),eps,y,x,rho
	double precision error,emax,t0x,t0y,tLx,tLy
	integer ICONTROL,i,j,Nx,Ny,PAS
	common/inicials/t0y,t0x,tLy,tLx
	! Iniciam un contador de passos
	PAS=1
!	Iniciam bucle del metode,aplicam la condicio de convergencia
	emax=1.d0
	do while (abs(emax).GT.eps)
		emax=0.d0
		call poisson2dpas(Nx,Ny,h,T,rho,ICONTROL,emax)
		if (PAS.GT.40000) then
			write(*,*) "no ha convergit"
			exit
		else
			PAS=PAS+1
		end if
	end do
	write(*,*)"ha convergit amb pas:",PAS
	end

! Definim la funcio rho
function rho1(x,y)
	implicit none
	double precision x,y,rho1,r 
	r=sqrt((x-22.5d0)**2+(y-8.d0)**2)
	rho1=10.d0*exp(-(r-4.d0)**2/0.7d0**2)
	return
end

function rho2(x,y)
	implicit none
	double precision x,y,rho2
	if ((x.GT.29.d0).and.(x.LT.35.d0).and.(y.GT.18.d0).and.(y.LT.22.d0)) then
		rho2=7.d0
	else
		rho2=0.d0 
	end if
end

function rho3(x,y)
	implicit none
	double precision x,y,rho3,r 
	r=sqrt((x-10.5d0)**2+(y-22.d0)**2)
	rho3=5.5d0*exp(-(r-5.d0)**2/1.2d0**2)
	return
end

function rho(x,y)
	implicit none
	double precision x,y,rho1,rho2,rho3,rho
	rho=rho1(x,y)+rho2(x,y)+rho3(x,y)
	return
end

! Si volem llevar les fonts empram aquesta funcio
function nofont(x,y)
	implicit none
	double precision nofont,x,y
	nofont=0
	return
end