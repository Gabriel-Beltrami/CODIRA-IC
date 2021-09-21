!--------------------------------------------------------------------------------------------------
!     MODULE   Visual
!>    @ingroup RANS
!!    @authors Jan Mateu
!!    @date    22/08/2017
!--------------------------------------------------------------------------------------------------
module visual
	use user
	implicit none

contains


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!	Calculo del Nusselt,temperaturas y velocidades adimensionales				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Comparison with the Literature
!>    @brief     Calculates and writes the respective velocities in the mean horizontal plan and left side Nusselt numbers, in other to compare them with the ones found in other papers.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     In Simple iteration, x-coordinate, both convection coefficients, time delta, dimensionless time, both velocities, temperature, mean pressure, initial pressure, CONTERo, CONTER, Mu and MuT.
!!    @param     Out Vertical velocity in the mean horizontal plan, simple interpolated mean Nusselt number, double interpolated Nusselt number, Max Mut, Max Mut adim, maximum vertical velocity, maximum horizontal velocity and maximum Nusselt.
!--------------------------------------------------------------------------------------------------
subroutine literatura_compare(iter_t,X,Hy,Hx,DeltaT,tadim,U,V,T,Pmed,Po,CONTERo,CONTER,Mu,MuT)

	implicit none

	integer :: i,j,iter_t
	real(KIND = DP) :: Hy,Hx, Numedio, Numedio2, DeltaT,tadim,Pmed,Po,dtdx,CONTERo,MuT_max, &
				u_max,V_max, U_adim, nu_max
	real(KIND = DP), dimension(Nx) :: VV,TT,X,Nu,Nu2
	real(KIND = DP), dimension(Nx,Ny) :: U,V,T,CONTER,Mu,MuT

! Calculo de la velocidad vertical en el plano horizontal medio

!	if ((iter_t.eq.100).or.(iter_t.eq.1000).or.(iter_t.eq.2000)) then

	U_adim = sqrt(Hy*DeltaT*g*bbeta )
		Do i=1,Nx
			VV(i)= V(i,Ny/2) / U_adim
		end do

! Calculo de la temperatura en el plano horizontal medio

!		Do i=1,Nx
!			TT(i)= (T(i,Ny/2)-Tc )/( Th-Tc )
!		end do
		Write(2,*)
		Write(2,*)
! Write T and VV
!		Write(2,*) tadim

		do i=1,Nx
			Write(2,*) (X(i)/Hx),VV(i)
	   	end do
!	end if

! Calculo de Nusselt medio en la pared izquierda

	NUmedio=0.0D+00
	Numedio2=0.0D+00


	Do j=1,Ny-2
		NU(j)=( CONTER(2,j+1)*(T(2,j+1)-Th)*Hx )/(-DeltaT*CONTERo*( X(2)-X(1) ))

		NUmedio=NUmedio+NU(j)*DYP(j+1)

		dtdx=Th*( -X(2)-X(3) )/( X(2)*X(3) )+ T(2,j+1)*(-X(3))/( X(2)*(X(2)-X(3)) ) &
			                            + T(3,j+1)*(-X(2))/( X(3)*(X(3)-X(2)) )
		NU2(j)=dtdx*Hx/(-DeltaT)
!		NU2(j)=( CONTER(1,j+1)*Hx )/(-DeltaT*CONTERo ) * (T(2,j+1)-Th) / ( X(2)-X(1) )
		Numedio2=Numedio2+Nu2(j)*DYP(j+1)
   	end do

	nu_max = maxval(NU)

	NUmedio=NUmedio/Hx
	NUmedio2=Numedio2/Hx

	Mut_max =  maxval(maxval(Mut,2))

	v_max =  maxval(maxval(v/U_adim,2))
	u_max =  maxval(maxval(u/U_adim,2))

	write(3,*) "NUmedio intepolacion simple =", NUmedio, "	|	","NUmedio intepolacion doble =", Numedio2
	write(3,*) "Max Mut", Mut_max
	write(3,*) "Max Mut adim", Mut_max / Muo

	write(3,*) "Max v", v_max
!	write(3,*) "Max Mut adim", Mut_max / Muo

	write(3,*) "Max u", u_max
!	write(3,*) "Max Mut adim", Mut_max / Muo
	write(3,*) "nu_max", nu_max
end subroutine literatura_compare


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing 					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Variable Fields Data Writing
!>    @brief     Write output data in the VTK format, to be opened by Paraview or Visit (for example).
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Horizontal velocity, vertical velocity, pressure, temperature, b, time iteration, Tk, e and MuT.
!!    @param     [out] X-coordinate, y-coordinate, z-coordinate, cell data, scalar pressure, temperature, e, Tk, MuT and str in table form and vector velocity.
!--------------------------------------------------------------------------------------------------
subroutine write_vtk(U,V,P,T,b,iter_t,TK,e,MuT)

	implicit none

	character(len = 27) :: str, tmp
	integer :: i,j,iter_t
	real(KIND = DP), dimension(Nx+1) :: x_print
	real(KIND = DP), dimension(Ny+1) :: Y_print
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,b,PSIP,TK,e,inter_U,inter_V,MuT

	x_print(1)=XU(1)
	do i = 1,Nx-1
		x_print(i+1)=XU(i)
	end do
	x_print(Nx+1)=XU(Nx-1)+1.0D-5

	y_print(1)=yv(1)
	do j = 1,Ny-1
		y_print(j+1)=yv(j)
	end do
	y_print(Ny+1)=yv(Ny-1)+1.0D-5
	x_print(2)=XU(1)+1.0D-5
	y_print(2)=yv(1)+1.0D-5


	write (tmp,'(I5.5)') iter_t
	str = 'outputdata/output_'//trim(tmp)//'.vtk'
	open(unit=1, file = str)

	write ( 1,'(a)' )'# vtk DataFile Version 3.0'
	write ( 1,'(a)' ) 'Non-uniform Rectilinear - Rectilinear Grid'
	write ( 1,'(a)' ) 'ASCII'
	write ( 1,'(a)' ) 'DATASET RECTILINEAR_GRID'
	write ( 1,'(a,3i4)' )'DIMENSIONS', nx+1,ny+1,1
	write ( 1,'(a13,1i4,a6)' ) 'X_COORDINATES',nx+1,' float'
	write ( 1,'(f22.14)' ) (x_print(i),i=1,nx+1)
	write ( 1,'(a13,1i4,a6)' ) 'Y_COORDINATES',ny+1,' float'
	write ( 1,'(f22.14)' ) (y_print(j),j=1,ny+1)
	write ( 1,'(a13,1i4,a6)' ) 'Z_COORDINATES',1,' float'
	write ( 1,'(f22.14)' ) 1.0

	write ( 1, '(a, I16)' ) 'CELL_DATA ', (nx)*(ny)*1

	write ( 1, '(a)' ) 'SCALARS P float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((P(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS T float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((T(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS e float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((e(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS Tk float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((TK(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS MuT float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((MuT(i,j),i=1,(nx)),j=1,(ny))

      	CALL STREAMLINES(DYP,U,PSIP,RHO)

	write ( 1, '(a)' ) 'SCALARS Str float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((PSIP(i,j),i=1,(nx)),j=1,(ny))

	CALL interpol_UV(U,V,inter_U,inter_V)

	write ( 1, '(a)' ) 'VECTORS velocity float'
	do j=1,ny
		do i=1,nx
			write ( 1, "(f22.14,f22.14,f22.14)" ) inter_U(i,j), inter_V(i,j), 0.0D+00
		end do
	end do

!	write ( 1, '(a)' ) 'SCALARS u float'
!	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
!	write ( 1, '(f22.14)' ) ((U(i,j),i=1,(nx-1)),j=1,(ny-1))

!	write ( 1, '(a)' ) 'SCALARS v float'
!	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
!	write ( 1, '(f22.14)' ) ((V(i,j),i=1,(nx-1)),j=1,(ny-1))



	close(1)

!	TIMECODE=TIMEF()

end subroutine write_vtk

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					print_Residual						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Print Residual
!>    @brief     Subroutine to print during program execution the residuals of every equation in order to monitor the convergence.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Simple iteration, time iteration, global errors of horizontal velocity, vertical velocity, pressure, temperature, Tk and e and maximum local erros of horizontal velocity, vertical velocity, pressure, temperature, Tk and e.
!!    @param     [out] Printed simple iteration, global errors of horizontal velocity, vertical velocity, pressure, temperature, Tk and e and maximum local erros of horizontal velocity, vertical velocity, pressure, temperature, Tk and e.
!--------------------------------------------------------------------------------------------------
Subroutine print_Residual(iter,iter_t,Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Rmax_T,	&
					Res_Tk,Rmax_Tk,Res_e,Rmax_e)

	real(KIND = DP) :: Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Rmax_T,Res_Tk,Rmax_Tk,Res_e,Rmax_e
	integer :: iter,iter_t

!	write(*,*) 'Iteracao no tempo:',iter_t
	write(*,*) "______________________________________________"
	write(*,*) 'Iteracao no SIMPLE:',iter
	write(*,*)
	write(*,*) 'Global error of U:',Res_U,"  |  ",'MAX local error of U:',Rmax_U
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of V:',Res_V,"  |  ",'MAX local error of V:',Rmax_V
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of P:',Res_P,"  |  ",'MAX local error of P:',Rmax_P
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of T:',Res_T,"  |  ",'MAX local error of T:',Rmax_T
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of TK:',Res_TK,"  |  ",'MAX local error of TK:',Rmax_TK
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of e:',Res_e,"  |  ",'MAX local error of e:',Rmax_e
	return

end Subroutine print_Residual

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing for Matlab				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE interpol_UV
!>    @brief     Interpolate the velocity values (u and v) at the Principal mesh.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Horizontal velocity and vertical velocity.
!!    @param     [out] Horizontal velocity interpolation and vertical velocity interpolation.
!--------------------------------------------------------------------------------------------------
subroutine interpol_UV(U,V,inter_U,inter_V)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,V,inter_U,inter_V


	do j=1,Ny
		do i=2,Nx-1
			inter_U(i,j) = (U(i,j)+U(i-1,j)) /2.0D+0
		end do
	end do
	do j=1,Ny
		inter_U(1,j) = U(1,j)
		inter_U(Nx,j) = U(Nx-1,j)
	end do

	do j=2,Ny-1
		do i=1,Nx
			inter_V(i,j) = (V(i,j)+V(i,j-1)) /2.0D+0
		end do
	end do

	do i=1,nx
		inter_V(i,1) = V(i,1)
		inter_V(i,Ny) = V(i,Ny-1)
	end do

end subroutine interpol_UV
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing for Matlab				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Other Variable Fields Data Writing for Matlab
!>    @brief     Data output Writing for Matlab.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
subroutine write_matlab(U,V,P,T,b,DYP,RHO)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,b,PSIP,RHO,inter_U,inter_V

	open(unit=1, file ='./outputdata/x.dat')
	open(unit=11, file ='./outputdata/y.dat')
	open(unit=9, file ='./outputdata/temperature.dat')
	open(unit=19, file ='./outputdata/stream.dat')
	open(unit=29, file ='./outputdata/u_vel.dat')

      	CALL STREAMLINES(DYP,U,PSIP,RHO)

	do i=1,Nx
		write(1,*) X(i)
	end do

	do j=1,Ny
		write(11,*) y(j)
	end do

	do j=1,Ny
		write(9,*) T(:,j)
	end do

	do j=1,Ny
		write(19,*) PSIP(:,j)
	end do

	CALL interpol_UV(U,V,inter_U,inter_V)

	do j=1,Ny
		write(29,*) inter_U(:,j)
	end do

	do j=1,Ny
		write(29,*) inter_V(:,j)
	end do

	close(1)
        close(9)
	close(11)
        close(19)
        close(29)
end subroutine write_matlab

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing for TECPLOT				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Variable Fields Data Writing for TecPlot
!>    @brief     Output data to be read by tecplot.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] In Filling.
!!    @param     [out] Out Filling.
!--------------------------------------------------------------------------------------------------
subroutine write_tecplot(U,V,P,T,RHO,MuT,TK,e)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,PSIP,RHO,MuT,TK,e,inter_U,inter_V

	open(unit=100, file ='./outputdata/tecplot.dat')


      	CALL STREAMLINES(DYP,U,PSIP,RHO)
	CALL interpol_UV(U,V,inter_U,inter_V)

	write(100,*)"Title= RANS"
 	write(100,*)"Variables=X,Y,u,v,PSIP,T,P,TK,Mut"
	write(100,*)"Zone I=",Nx, "J=" ,Ny, "F=point"
	write(100,*)

	do  j=1,Ny
		do i=1,Nx

		Write(100,"(200F25.15)") X(i),Y(j),inter_U(i,j),inter_V(i,j),PSIP(i,j),T(i,j),	&
		P(i,j),TK(i,j),e(i,j),Mut(i,j)

		end do
	end do


	close(100)

end subroutine write_tecplot


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			STREAM LINES (PSIP = LINEAS DE CORRIENTE) 				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Streamlines
!>    @brief     Computes the streamlines in both velocity and main meshes.
!!    @authors   Jan Mateu
!!    @date      08/22/2017
!!    @param     [in] Filling.
!!    @param     [out] Filling.
!--------------------------------------------------------------------------------------------------
SUBROUTINE STREAMLINES(DYP,U,PSIP,RHO)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,RHO,PSI,PSIP

	DO I=1,Nx
	DO J=1,Ny
		PSI(I,J)=0.0D+00
		PSIP(I,J)=0.0D+00
	end do
	end do

!  NOTA: CALCULO DE PSI SOBRE LA MALLA DE LA VELOCIDAD "V"

        DO I=2,NX-1
        DO J=2,NY-2

              PSI(I,J)=(((U(I,J)+U(I-1,J))/2.0D+00)*DYP(J) * RHO(i,j)*1000 )+PSI(I,J-1)

	END DO
	END DO

!  NOTA: CALCULO DE PSI SOBRE LA MALLA PRINCIPAL (T,P)

        DO I=2,NX-1
        DO J=3,NY-2

                  PSIP(I,J)=(PSI(I,J)+PSI(I,J-1))/2.0D+00

	END DO
	END DO

!     LADO SUR
        DO I=2,NX-1
               J=2
                  PSIP(I,J)=(PSI(I,J)+PSI(I,J-1))/2.0D+00
	END DO

!     LADO NORTE
        DO I=2,NX-1
               J=NY-1
                  PSIP(I,J)=(PSI(I,J)+PSI(I,J-1))/2.0D+00
	END DO

      RETURN
END SUBROUTINE STREAMLINES


end module visual
