module visual
	use user
	implicit none

contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing 					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine write_vtk(U,V,T)

	implicit none

	character(len = 27) :: str, tmp
	integer :: i,j
	real(KIND = DP), dimension(Nx+1) :: x_print
	real(KIND = DP), dimension(Ny+1) :: Y_print
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,b,PSIP,TK,e,inter_U,inter_V,MuT,	&
										 Pk,epsi,div_T,DTK_Dt,vorticity_z, P_ent

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

	str = 'outputdata/output.vtk'
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

	write ( 1, '(a)' ) 'SCALARS T float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((T(i,j),i=1,(nx)),j=1,(ny))

    CALL STREAMLINES(DYP,U,PSIP,RHOo)

	write ( 1, '(a)' ) 'SCALARS Str float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((PSIP(i,j),i=1,(nx)),j=1,(ny))
	CALL interpol_UV(U,V,inter_U,inter_V)

	write ( 1, '(a)' ) 'SCALARS u float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.8)' ) ((inter_U(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS v float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.8)' ) ((inter_V(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'VECTORS velocity float'
	do j=1,ny
		do i=1,nx
			write ( 1, "(f22.14,f22.14,f22.14)" ) inter_U(i,j), inter_V(i,j), 0.0D+00
		end do
	end do


	close(1)

!	TIMECODE=TIMEF()

end subroutine write_vtk

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					print_Residual						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine print_Residual(iter,Res_T,Rmax_T)

	real(KIND = DP) :: Res_T,Rmax_T
	integer :: iter

!	write(*,*) 'Iteracao no tempo:',iter_t
	write(*,*) "______________________________________________"
	write(*,*) 'Iteracao no SIMPLE:',iter
	write(*,*)
	write(*,*) 'Global error of T:',Res_T,"  |  ",'MAX local error of T:',Rmax_T
	write(*,*) "-----------------------------------------------"
	return

end Subroutine print_Residual

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!!							Interpoling stagared values to the main mesh			!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine interpol_UV(U,V,inter_U,inter_V)

	implicit none

	integer :: i,j
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

SUBROUTINE STREAMLINES(DYP,U,PSIP,RHOo)

	implicit none

	integer :: i,j
	real(kind = dp) :: rhoo
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,PSI,PSIP

	DO I=1,Nx
	DO J=1,Ny
		PSI(I,J)=0.0D+00
		PSIP(I,J)=0.0D+00
	end do
	end do

!  NOTA: CALCULO DE PSI SOBRE LA MALLA DE LA VELOCIDAD "V"

        DO I=2,NX-1
        DO J=2,NY-2

              PSI(I,J)=(((U(I,J)+U(I-1,J))/2.0D+00)*DYP(J) * RHOo*1000 )+PSI(I,J-1)

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
