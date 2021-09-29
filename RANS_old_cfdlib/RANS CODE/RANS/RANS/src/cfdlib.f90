!--------------------------------------------------------------------------------------------------
!     MODULE   Cfd Lib
!>    @ingroup RANS
!--------------------------------------------------------------------------------------------------
module cfdlib
	use user
	implicit none


contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					start (U,V,P)						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE U, V and P Start
!>    @brief     Define the initial solution of all unknown variables for the iterative process.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
subroutine start (U,V,T,T_ref, Told)

	implicit none
	integer:: i,j
	real(KIND = DP)  :: T_ref
	real(KIND = DP), dimension(Nx,Ny) :: U,V,T, Told


	Do j=1,Ny
		Do i=1,(Nx-1)

			U(i,j)=U_const
		end do
	end do


	Do j=1,(Ny-1)
		Do i=1,Nx
			V(i,j)=V_const
		end do
	end do

	Do j=1,Ny
		Do i=1,Nx
			T(i,j)=T_ref
			Told(i, j) = T_ref
		end do
	end do


    return
end subroutine start




!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					update					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE U Update
!>    @brief     Subroutine to copy the values of a matrix U into a matrix Ux, both of size Nx-1 * Ny.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
subroutine update_U (Ux,U)

	implicit none
	integer:: i,j
	real(KIND = DP), dimension(Nx,Ny) :: Ux,U

	Do j=1,Ny
		Do i=1,(Nx-1)
			Ux(i,j)=U(i,j)
		end do
	end do

    return
end subroutine update_U
!					update_V (Ux,U)						!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE V Update
!>    @brief     Subroutine to copy the values of a matrix V into a matrix Vx, both of size Nx * Ny-1.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
subroutine update_V (Ux,U)

	implicit none
	integer:: i,j
	real(KIND = DP), dimension(Nx,Ny) :: Ux,U

	Do j=1,(Ny-1)
		Do i=1,Nx
			Ux(i,j)=U(i,j)
		end do
	end do

    return
end subroutine update_V

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_T 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE T Coefficient
!>    @brief     Compute the coefficients to solve the energy transport equation for each node in the domain.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
subroutine COEF_T(DXP,DYP,DXU,DYV,X,Y,XU,YV,			&
		CONTERo,Cpo,RHOo,U,V,T,ap,ae,aw,as,an,b,Res_T,Rmax_T, dt, T_s, T_n, T_e, T_w)

	implicit none

	integer:: i,j,npas

	real(KIND = DP)  :: 			&
		Fe,Fw,Fn,Fs,			&
		De,Dw,Dn,Ds,			&
		RES_T,resi,summ,Rmax_T,		&
		Pes,Pen,Pee,Pew,		&
		dt,				&
		GIE,GIW,GIN,GIS,		&
		GIEL,GIWL,GINL,GISL,		&
		GIET,GIWT,GINT,GIST, gamma, contero,Cpo,rhoo, T_W, T_E, T_S, T_N

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X,XU,DXU

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y,YV,DYV

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b,apo,		&
		U,V,T,		&
		Fep,Fnp,MUt

!-----------------------------------------------------------------------------------------------!
!					varaibles reset 		!
!-----------------------------------------------------------------------------------------------!


	do j=1,Ny
		do i=1,Nx

			ae(i,j)=0.0D+00

			aw(i,j)=0.0D+00

			an(i,j)=0.0D+00

			as(i,j)=0.0D+00

			ap(i,j)=0.0D+00

			b(i,j)=0.0D+00

			apo(i,j)=0.0D+00

		end do
	end do



	do  j=2,Ny-1
		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

!      Fe=Fep(i,j)*DYP(j)

!      Fw=Fep(i-1,j)*DYP(j)

!      Fn=Fnp(i,j)*DXP(i)

!      Fs=Fnp(i,j-1)*DXP(i)

      Fe=U(i+1,j)*RHOO*DYP(i)

      Fw=U(i,j)*RHOO*DYP(i)



      Fn=V(i,j+1)*Rhoo*DXP(j)

      Fs=V(i,j)*Rhoo*DXP(j)

!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

 Gamma = CONTERo/Cpo

!	Calculo do flux difusivo

	De=Gamma*DYP(j)/( (x(i+1)-x(i)) )

	Dw=Gamma*DYP(j)/( (x(i)-x(i-1)) )

	Dn=Gamma*DXP(i)/( (y(j+1)-y(j)) )

	Ds=Gamma*DXP(i)/( (y(j)-y(j-1)) )

!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


	Pee=Fe/De

      	Pew=Fw/Dw

	Pen=Fn/Dn

	Pes=Fs/Ds


!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!



	ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

	an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

	as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

	apo(i,j)=(RHOo*DXP(i)*DYP(j))/dt
	apo(i,j) = 0
     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j)
      	b(i,j)=apo(i,j)*T(i,j)

!if ((i==12).and.(j==10)) then

!write(*,*) "Fw"
!write(*,*) Fw
!write(*,*) "DMAX1(Fw,0.0D+00)"
!write(*,*) DMAX1(Fw,0.0D+00)

!stop
!end if

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
		b(i,1)=T_S
	end do

!					NORTH					!
	Do i=1,(Nx)
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00
		b(i,Ny)=T_N
	end do


!					WEST					!
	Do j=2,Ny-1
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j)=T_W
	end do


!					EAST					!
	Do j=2,Ny-1
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=0.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j)=T_E
	END DO


!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


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

end Subroutine coef_T




!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						reset(Pc)					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Pc Reset
!>    @brief     Filling.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
Subroutine reset(Pc)

	real(KIND = DP), dimension(Nx,Ny) :: Pc
	integer :: i,j

     	Do j=1,Ny
		Do i=1,Nx

			Pc(i,j)=1.0D-30
		end do
	end do

	return
end Subroutine reset

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						converegence					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     FUNCTION   Convergence
!>    @brief     Check convergence of all equation.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
function converegence(Res_T,iter)


	real(KIND = DP):: Res_T, Rmax
	logical:: converegence
	integer:: iter


	Rmax = Res_T
	if (Rmax .gt. 1.0D+20) stop 'A friendly message: diverged'

	converegence=(Res_T.LE.epsilonT)
	return
end

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
