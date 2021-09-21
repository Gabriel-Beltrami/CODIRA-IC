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

	T(1)=0
	T(Nx)=100

	Ti(1)=0
	Ti(Nx)=100

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
	         ap,ae,aw,b)

	implicit none

	integer:: i,P2

	real(KIND = DP)  :: 			&
		G,G1,G2,DX,P1,P	

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


	P2=int((Hx1-DXP(2))/DX)+2
	
	write(*,*) P2


		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!

!	 Calculo da gamma
	G1=Lambda1/Cp1
	G2=Lambda2/Cp2
	
	if (i<=P2) then

		ae(i)=G1/DXP(I+1)

		aw(i)=G1/DXP(I)

	     	ap(i)=ae(i)+aw(i)

	      	b(i)=0
      	
      	end if
            	
      	if (i>P2) then
      	
	      	ae(i)=G2/DXP(I+1)

		aw(i)=G2/DXP(I)

	     	ap(i)=ae(i)+aw(i)

	      	b(i)=0      	

	end if
	     	
		end do
	
		P1=real(int((Hx1-DXP(2))/DX)+0.5)
	
		P=Hx1-P1*DX
	 	
	 	write (*,*) DX
	 	write (*,*) Hx1
	 	write (*,*) P1
	 	write (*,*) P
	 	
	 	G=(G1*G2*DX)/(G2*P+G1*(DX-P))
	 	
	 	!G=G1+P/DX*(G2-G1)
	 	
	 	write (*,*) G1
	 	write (*,*) G2
	 	write (*,*) G
	 	write (*,*) DXP(4)
      	
      		ae(P2)=G/DXP(P2+1)
      		
      		write (*,*) ae(P2)

		aw(P2+1)=G/DXP(P2+1)

	     	ap(P2)=ae(P2)+aw(P2)

	     	ap(P2+1)=ae(P2+1)+aw(P2+1)

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

end module cfdlib
