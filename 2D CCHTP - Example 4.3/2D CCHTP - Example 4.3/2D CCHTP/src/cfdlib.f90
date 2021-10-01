module cfdlib
	use user
	implicit none


contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					start (T)						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine start (T,Ti,Tref)

	implicit none
	integer:: i,j
	real(KIND = DP)  :: Tref
	real(KIND = DP), dimension(Nx,Ny) :: T,Ti

	Tref=(Ta+Tb)/2.0D+00
	
	DO i=2,Nx
	
		T(i,1)=Ta
	
	END DO
	 
	DO i=2,Nx
	
		Ti(i,1)=Ta
	
	END DO	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	T(1,Ny)=Tb
	
	DO i=2,Nx-1

		T(i,Ny)=Tb
	
	END DO
			
	T(Nx,Ny)=Ta

	Ti(1,Ny)=Tb
	
	DO i=2,Nx-1

		Ti(i,Ny)=Tb
	
	END DO
			
	Ti(Nx,Ny)=Ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO j=1,Ny-1

		T(1,j)=Tb
	
	END DO

	DO j=1,Ny-1

		Ti(1,j)=Tb
	
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO j=2,Ny-1

		T(Nx,j)=Ta
	
	END DO
	
	DO j=2,Ny-1

		Ti(Nx,j)=Ta
	
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Do j=2,Ny-1
		Do i=2,Nx-1
			
			T(i,j)=Tref
			Ti(i,j)=Tref
		
		end do
	end do

	!write(*,*) Ti
	!write(*,*) T

    return
end subroutine start

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_T 					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_T(DXP,DYP,DX,DY,X,Y,		&
	         ap,ae,aw,an,as,b,Fep,Fnp)

	implicit none

	integer:: i,j

	real(KIND = DP)  :: 			&
		G,DX,DY,Fe,Fw,Fs,Fn,De,	&
		Dw,Ds,Dn,Pee,Pew,Pes,Pen,	&
		xu1,xu2,xu3,yv1,yv2,yv3	

	real(KIND = DP),intent(in) :: 	&	
		Fep,Fnp

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b

!-----------------------------------------------------------------------------------------------!
!					varaibles reset					    !
!-----------------------------------------------------------------------------------------------!

	do j=1,Ny
		do i=1,Nx

			ae(i,j)=0.0D+00

			aw(i,j)=0.0D+00

			an(i,j)=0.0D+00

			as(i,j)=0.0D+00			
			
			ap(i,j)=0.0D+00

			b(i,j)=0.0D+00

		end do
	end do

	do j=2,Ny-1
		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

!      Fe=Fep

!      Fw=Fep
      
!      Fs=Fnp
      
!      Fn=Fnp

      Fe=Fep*DY

      Fw=Fep*DY

      Fn=Fnp*DX

      Fs=Fnp*DX

!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

!	 Calculo da gamma
	G=Lambda/Cp

!	Calculo do flux difusivo

	De=G/DXP(I+1)*DY

	Dw=G/DXP(I)*DY
	
	Dn=G/DYP(J+1)*DX
	
	Ds=G/DYP(J)*DX
	
	!De=G*DYP(j)/( (x(i+1)-x(i)) )

	!Dw=G*DYP(j)/( (x(i)-x(i-1)) )

	!Dn=G*DXP(i)/( (y(j+1)-y(j)) )

	!Ds=G*DXP(i)/( (y(j)-y(j-1)) )

!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!

	Pee=Fe/De

      	Pew=Fw/Dw
      	
      	Pen=Fn/Dn
      	
      	Pes=Fs/Ds

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!	

	ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)
	
	an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

	as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)

      	b(i,j)=0	
      	
      	
      	!if (j==3 .and. i==2) then
      	!write(*,*) Dw
      	!write(*,*) De
      	!write (*,*) Dn*A(n,Pen)
      	!write (*,*) Ds*A(n,Pes)
      	!write(*,*) Fn
      	!write(*,*) Fs
      	!write(*,*) Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)
      	!write(*,*) Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)
      	!write(*,*) Pew
      	!write(*,*) Pee	
      	!end if

		end do
	end do
!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!

!					SOUTH					!
	Do i=2,Nx-1
		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=0.0D+00
		as(i,1)=0.0D+00
		b(i,1)=Ta
	end do

!					NORTH					!
	Do i=2,Nx-1
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00
		b(i,Ny)=Tb
		!if (i==1 .or. i==Nx) b(i,Ny)=0.0D+00
		!if (i/=1 .and. i/=Nx) b(i,Ny)=100.0D+00
	end do


!					WEST					!
	Do j=1,Ny
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j)=Tb
	end do


!					EAST					!
	Do j=1,Ny
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=0.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j)=Ta
	END DO
	
end Subroutine coef_T

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						converegence					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function converegence(Res_T,iter)


	real(KIND = DP):: Res_T,Rmax
	logical:: converegence
	integer:: iter

	Rmax = Res_T
	if (Rmax .gt. 1.0D+20) stop 'A friendly message: diverged'

	converegence=(Res_T.LE.epsilonT)

	return
end

!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!

subroutine residual(T,ap,ae,aw,an,as,b,Res_T,Rmax_T)

	real(KIND = DP)  :: 			&
		RES_T,resi,summ,Rmax_T

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,ap,an,as,b
	
	real(KIND = DP), dimension(Nx,Ny) :: T

	summ=0.0D+00
	Rmax_T=0.0D+00

	do  j=2,Ny-1
		do i=2,Nx-1

			Resi=dabs(ap(i,j)*T(i,j)-ae(i,j)*T(i+1,j)-aw(i,j)*T(i-1,j)-as(i,j)*T(i,j-1)-an(i,j)*T(i,j+1)-b(i,j))
			summ=summ+(Resi*Resi)

			if (Resi .gt. Rmax_T) Rmax_T=Resi

		end do
	end do

	Res_T=dsqrt(summ)

	return	
end subroutine residual

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Mass_flux				                           !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine Mass_Flux(u,v,Rho,Fep,Fnp)

	implicit none
	real(KIND = DP),intent(in) :: u,v,Rho
	real(KIND = DP),intent(out) :: Fep,Fnp

	Fep=u*Rho
	Fnp=v*Rho

    return
end subroutine Mass_Flux

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						 A(n,Pe)					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function A(n,Pe)

	real(KIND = DP)  :: Pe,A
	integer :: n


	if (n.EQ.1) A=1.0D+00-0.5D+00*DABS(Pe)
	if (n.EQ.2) A=1.0D+00
	if (n.EQ.3) A=DMAX1( 0.0D+00,1.0D+00-( 5.0D-01*DABS(Pe)) )
	if (n.EQ.4) A=DMAX1(0.0D+00,(1.0D+00-( 1.0D-01*DABS(Pe) ))**5.0D+0)
	if (n.EQ.5) A=DABS(Pe)/(DEXP(ABS(Pe))-1.0D+00)

	return
end

end module cfdlib
