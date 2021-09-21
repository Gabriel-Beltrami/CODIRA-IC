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
	
	T(:,1)=Ta
	
	Ti(:,1)=Ta
	
	T(1,Ny)=Ta
	
	DO i=2,Nx-1

		T(i,Ny)=Tb
	
	END DO
			
	T(Nx,Ny)=Ta

	Ti(1,Ny)=Ta
	
	DO i=2,Nx-1

		Ti(i,Ny)=Tb
	
	END DO
			
	Ti(Nx,Ny)=Ta

!	Ti(1,5)=Ta
!	Ti(2,5)=Tb
!	Ti(3,5)=Tb
!	Ti(4,5)=Tb
!	Ti(5,5)=Ta

	DO j=2,Ny-1

		T(1,j)=Ta
	
	END DO

	DO j=2,Ny-1

		Ti(1,j)=Ta
	
	END DO
		
!	Ti(1,2)=Ta
!	Ti(1,3)=Ta
!	Ti(1,4)=Ta

	DO j=2,Ny-1

		T(Nx,j)=Ta
	
	END DO
	
	DO j=2,Ny-1

		Ti(Nx,j)=Ta
	
	END DO

!	Ti(5,2)=Ta
!	Ti(5,3)=Ta
!	Ti(5,4)=Ta

	Do j=2,Ny-1
		Do i=2,Nx-1
			
			Ti(i,j)=Tref
		
		end do
	end do

    return
end subroutine start

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_T 					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_T(DXP,DYP,DX,DY,X,Y,		&
	         ap,ae,aw,an,as,b)

	implicit none

	integer:: i,j

	real(KIND = DP)  :: 			&
		G,DX,DY	

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
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!

!	 Calculo da gamma
	G=Lambda/Cp
	

	ae(i,j)=G/DXP(I+1)*DY

	aw(i,j)=G/DXP(I)*DY

	an(i,j)=G/DYP(J+1)*DX

	as(i,j)=G/DYP(J)*DX

     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)

      	b(i,j)=0

		end do
	end do
!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!

!					SOUTH					!
	Do i=1,(Nx)
		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=0.0D+00
		as(i,1)=0.0D+00
		b(i,1)=0.0D+00
	end do

!					NORTH					!
	Do i=1,(Nx)
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00
		if (i==1 .or. i==Nx) b(i,Ny)=0.0D+00
		if (i/=1 .and. i/=Nx) b(i,Ny)=100.0D+00
	end do


!					WEST					!
	Do j=2,Ny-1
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j)=0.0D+00
	end do


!					EAST					!
	Do j=2,Ny-1
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=0.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j)=0.0D+00
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

end module cfdlib
