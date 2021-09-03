module cfdlib
	use user
	implicit none


contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					start (T)						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine start (T,X)

	implicit none
	integer:: i
	real(KIND = DP), dimension(Nx) :: X
	real(KIND = DP), dimension(Nx,itermax,itermax_t) :: T

	T(1,:,:)=0
	T(Nx,:,:)=0

	Do i=2,Nx-1
		T(i,1,1)=tin*SIN(PI*X(i))
		!print *,T(i,1,1)
	
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
	         ap,apo,ae,aw,b,T,iter_t)

	implicit none

	integer:: i,iter_t

	real(KIND = DP)  :: 			&
		G,DX	

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X

	real(KIND = DP), dimension(Nx) :: 	&
		ae,aw,apo,ap,b
		
	real(KIND = DP), dimension(Nx,itermax,itermax_t) :: T	

!-----------------------------------------------------------------------------------------------!
!					varaibles reset					    !
!-----------------------------------------------------------------------------------------------!


		do i=1,Nx

			ae(i)=0.0D+00

			aw(i)=0.0D+00

			ap(i)=0.0D+00
			
			apo(i)=0.0D+00

			b(i)=0.0D+00

		end do


		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!

!	 Calculo da gamma
	G=Lambda/Cp
	

	ae(i)=G/DXP(I+1)

	aw(i)=G/DXP(I)
     	 	
     	apo(i)=rho*DX/dt

	ap(i)=apo(i)+ae(i)+aw(i)

	if (iter_t==1 .or. iter_t==2) then

      	b(i)=apo(i)*T(i,1,1)

	else 
	
	b(i)=apo(i)*T(i,1,iter_t)

	end if
		end do

!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!


!					WEST					!

		ae(1)=0.0D+00
		aw(1)=0.0D+00
		ap(1)=1.0D+00
		b(1)=0.0D+00

!					EAST					!
		ae(Nx)=0.0D+00
		aw(Nx)=0.0D+00
		ap(Nx)=1.0D+00
		b(Nx)=0.0D+00
	
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

subroutine residual(T,ap,apo,ae,aw,b,Res_T,Rmax_T)

	real(KIND = DP)  :: 			&
		RES_T,resi,summ,Rmax_T

	real(KIND = DP), dimension(Nx) :: 	&
		ae,aw,ap,apo,b
	
	real(KIND = DP), dimension(Nx,itermax,itermax_t) :: T

	summ=0.0D+00
	Rmax_T=0.0D+00

	do i=2,Nx-1

		Resi=dabs(ap(i)*T(i,iter,iter_t)-ae(i)*T(i+1,iter,iter_t)-aw(i)*T(i-1,iter,iter_t)-b(i))
		summ=summ+(Resi*Resi)

		if (Resi .gt. Rmax_T) Rmax_T=Resi

		end do

	Res_T=dsqrt(summ)

	return	
end subroutine residual

end module cfdlib
