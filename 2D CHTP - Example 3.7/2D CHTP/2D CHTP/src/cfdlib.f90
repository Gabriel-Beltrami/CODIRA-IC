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
	real(KIND = DP)  :: Tref,DX
	real(KIND = DP), dimension(Nx,Ny) :: T,Ti

	Tref=(Ta+Tb)/2.0D+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	DO i=2,Nx-1

		T(i,1)=0
	
	END DO

!!!!!!!!!!!!!!!!!!

	DO i=2,Nx-1

		Ti(i,1)=0
	
	END DO
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO j=1,Ny

		T(1,j)=Ta
	
	END DO

!!!!!!!!!!!!!!!!!!

	DO j=1,Ny

		Ti(1,j)=Ta
	
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	DO j=1,Ny

		T(Nx,j)=Tb
	
	END DO

!!!!!!!!!!!!!!!!!!
	
	DO j=1,Ny

		Ti(Nx,j)=Tb
	
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Do j=2,Ny
		Do i=2,Nx-1
			
			Ti(i,j)=50.0D+00
			T(i,j)=50.0D+00
			
		end do
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!DO i=2,Nx-1

	!	if (i<k) then
		
	!		T(i,Ny)=T(i,Ny-1)-Q/(Lambda1*DX)*DYP(Ny)
	
	!	else if (i==k) then
		
	!		T(i,Ny)=((T(i,Ny-1)-Q/(Lambda1*DX)*DYP(Ny))+T(i,Ny-1)-Q/(Lambda2*DX)*DYP(Ny))/2
		
	!	else
		
	!		T(i,Ny)=T(i,Ny-1)-Q/(Lambda2*DX)*DYP(Ny)
	!	end if
		
	!END DO

!!!!!!!!!!!!!!!!!!
	
	!DO i=2,Nx-1

	!	if (i<k) then
		
	!		Ti(i,Ny)=Ti(i,Ny-1)-Q/(Lambda1*DX)*DYP(Ny)
	
	!	else if (i==k) then
		
	!		Ti(i,Ny)=((Ti(i,Ny-1)-Q/(Lambda1*DX)*DYP(Ny))+Ti(i,Ny-1)-Q/(Lambda2*DX)*DYP(Ny))/2
		
	!	else
		
	!		Ti(i,Ny)=Ti(i,Ny-1)-Q/(Lambda2*DX)*DYP(Ny)
	!	end if
		
	!END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
		G,G1,G2,DX,DY	

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
	G1=Lambda1/Cp1
	G2=Lambda2/Cp2
	
	if (i<k) then

		ae(i,j)=G1/DXP(I+1)*DY

		aw(i,j)=G1/DXP(I)*DY

		an(i,j)=G1/DYP(J+1)*DX

		as(i,j)=G1/DYP(J)*DX

	     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)

	      	b(i,j)=0

	end if

	if (i>k) then
	
		ae(i,j)=G2/DXP(I+1)*DY

		aw(i,j)=G2/DXP(I)*DY

		an(i,j)=G2/DYP(J+1)*DX

		as(i,j)=G2/DYP(J)*DX

	     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)

	      	b(i,j)=0

	end if
	
		end do
	end do
	
	
		G=(G1*G2*2.0D+00)/(G2+G1)
	 	
	! 	write (*,*) G1
	! 	write (*,*) G2
	! 	write (*,*) G
      	
      	!do j=2,Ny-1
      	
      	!	if (j/=k) then
      	
         !   		ae(k,j)=G2/DXP(k+1)*DY

      	!		aw(k,j)=G1/DXP(k)*DY

	!		if (j==k+1) then

	!			an(k,j)=G/DYP(J+1)*DX
	!			as(k,j)=G1/DYP(J)*DX
	
	!		else if (j==k-1) then
	!		
	!			an(k,j)=G1/DYP(J+1)*DX
	!			as(k,j)=G/DYP(J)*DX		
	!						
	!		else
			
	!			an(k,j)=G/DYP(J+1)*DX
	!			as(k,j)=G/DYP(J)*DX

	!		end if

	!	     	ap(k,j)=ae(k,j)+aw(k,j)+an(k,j)+as(k,j)
	
	!	      	b(k,j)=0
	      
	!	end if
	      	
	!end do

	!ae(k,k)=G/DXP(k+1)*DY

	!aw(k,k)=G1/DXP(k)*DY
	
	!an(k,k)=G1/DYP(k+1)*DX
	
	!as(k,k)=G1/DYP(k)*DX
	
	!b(k,k)=0

	!ap(k,k)=ae(k,k)+aw(k,k)+an(k,k)+as(k,k)
	
	!aw(k+1,k)=G/DXP(k+1)*DY
	
	!ap(k+1,k)=ae(k+1,k)+aw(k+1,k)+an(k+1,k)+as(k+1,k)
	
	
	do j=2,Ny-1
      	
    		ae(k,j)=G2/DXP(k+1)*DY

		aw(k,j)=G1/DXP(k)*DY	
			
		an(k,j)=G/DYP(J+1)*DX
	
		as(k,j)=G/DYP(J)*DX

	    	ap(k,j)=ae(k,j)+aw(k,j)+an(k,j)+as(k,j)

	      	b(k,j)=0
		      	      	
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
		b(i,1)=0.0D+00
	end do

!					NORTH					!
	Do i=2,Nx-1
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=1.0D+00
		b(i,Ny)=0.0D+00
		
		!if (i<k) then
		
		!	b(i,Ny)=Q/(Lambda1*DX)*DYP(Ny)
	
		!else if (i==k) then
		
		!	b(i,Ny)=((Q/(Lambda1*DX)*DYP(Ny))+(Q/(Lambda2*DX)*DYP(Ny)))/2
		
		!else
		
		!	b(i,Ny)=Q/(Lambda2*DX)*DYP(Ny)
			
		!end if
		
	end do


!					WEST					!
	Do j=1,Ny
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j)=Ta
	end do


!					EAST					!
	Do j=1,Ny
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=0.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j)=Tb
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
