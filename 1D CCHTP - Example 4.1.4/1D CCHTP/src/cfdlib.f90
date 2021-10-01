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
	integer:: i
	real(KIND = DP)  :: Tref
	real(KIND = DP), dimension(Nx) :: T,Ti

	Tref=(Ta+Tb)/2.0D+00

	T(1)=Ta
	T(Nx)=Tb

	Ti(1)=Ta
	Ti(Nx)=Tb

	Do i=2,Nx-1
		
		Ti(i)=Tref
	
	end do

    return
end subroutine start

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_T 					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_T(DXP,DX,X,			&
	         ap,ae,aw,b,Fep)

	implicit none

	integer:: i

	real(KIND = DP)  :: 			&
		DX,G,Fe,Fw,De,			&
		Dw,Pee,Pew
		
	real(KIND = DP),intent(in) :: 	&	
		Fep

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X

	real(KIND = DP), dimension(Nx) :: 	&
		ae,aw,ap,b

!-----------------------------------------------------------------------------------------------!
!					varaibles reset					    !
!-----------------------------------------------------------------------------------------------!


		do i=1,Nx

			ae(i)=0.0D+00

			aw(i)=0.0D+00

			ap(i)=0.0D+00

			b(i)=0.0D+00

		end do


		do  i=2,Nx-1
		
!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

      Fe=Fep

      Fw=Fep

!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

!	 Calculo da gamma
	G=Lambda/Cp

!	Calculo do flux difusivo

	De=G/DXP(I+1)

	Dw=G/DXP(I)

!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!

	Pee=Fe/De

      	Pew=Fw/Dw

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!	

	ae(i)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

     	ap(i)=ae(i)+aw(i)

      	b(i)=0	
      	
      	write(*,*) Dw
      	write(*,*) De
      	write(*,*) Fw
      	write(*,*) Fe
      	write(*,*) Pew
      	write(*,*) Pee
      
		end do

!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!

!					WEST					!

		ae(1)=0.0D+00
		aw(1)=0.0D+00
		ap(1)=1.0D+00
		b(1)=Ta

!					EAST					!
		ae(Nx)=0.0D+00
		aw(Nx)=0.0D+00
		ap(Nx)=1.0D+00
		b(Nx)=Tb
	
end Subroutine coef_T

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						converegence					    !
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
!					RESIDUAL						    !
!-----------------------------------------------------------------------------------------------!

subroutine residual(T,ap,ae,aw,b,Res_T,Rmax_T)

	real(KIND = DP)  :: 			&
		RES_T,resi,summ,Rmax_T

	real(KIND = DP), dimension(Nx) :: 	&
		ae,aw,ap,b
	
	real(KIND = DP), dimension(Nx) :: T

	summ=0.0D+00
	Rmax_T=0.0D+00

	do i=2,Nx-1

		Resi=dabs(ap(i)*T(i)-ae(i)*T(i+1)-aw(i)*T(i-1)-b(i))
		summ=summ+(Resi*Resi)

		if (Resi .gt. Rmax_T) Rmax_T=Resi

		end do

	Res_T=dsqrt(summ)

	return	
end subroutine residual

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Mass_flux				                           !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine Mass_Flux(u,Rho,Fep)

	implicit none
	real(KIND = DP),intent(in) :: u,Rho
	real(KIND = DP),intent(out) :: Fep

	Fep=u*Rho

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
