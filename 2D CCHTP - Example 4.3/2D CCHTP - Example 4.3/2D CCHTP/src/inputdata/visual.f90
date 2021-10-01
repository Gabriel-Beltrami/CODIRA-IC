module visual
	use user
	implicit none
	
contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing 					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine write_vtk(T,b)

	implicit none

	character(len = 27) :: str, tmp
	integer :: i,j,k
	real(KIND = DP), dimension(Nx) :: x_print
	real(KIND = DP), dimension(Ny) :: y_print
	real(KIND = DP), dimension(Nx,Ny) :: b
	real(KIND = DP), dimension(Nx,Ny) :: T

	x_print(1)=X(1)
	do i = 2,Nx-1
		x_print(i)=X(i)
	end do
	x_print(Nx)=X(Nx) 	

	y_print(1)=Y(1)
	do j = 2,Ny-1
		y_print(j)=Y(j)
	end do
	y_print(Ny)=Y(Ny) 

	str = 'outputdata/output.vtk'
	open(unit=1, file = str)

	write ( 1,'(a)' )'# vtk DataFile Version 3.0'
	write ( 1,'(a)' ) 'Uniform Rectilinear - Rectilinear Grid'
	write ( 1,'(a)' ) 'ASCII'
	write ( 1,'(a)' ) 'DATASET RECTILINEAR_GRID'
	write ( 1,'(a,3i4)' )'DIMENSIONS', nx+1,ny+1,1
	write ( 1,'(a13,1i4,a6)' ) 'X_COORDINATES',nx,' float'
	write ( 1,'(f22.14)' ) (x_print(i),i=1,nx)
	write ( 1,'(a13,1i4,a6)' ) 'Y_COORDINATES',ny,' float'
	write ( 1,'(f22.14)' ) (y_print(j),j=1,ny)
	write ( 1,'(a13,1i4,a6)' ) 'Z_COORDINATES',1,' float'
	write ( 1,'(f22.14)' ) 1.0
	!write ( 1,'(a13,1i4,a6)' ) 'COEFFICIENTS_DXP',nx,' float'
	!write ( 1,'(f22.14)' ) (DXP(i),i=1,nx)
	!write ( 1,'(a13,1i4,a6)' ) 'COEFFICIENTS_DYP',ny,' float'
	!write ( 1,'(f22.14)' ) (DYP(j),j=1,ny)
	
	write ( 1, '(a, I16)' ) 'CELL_DATA ', (nx)*(ny)*1
	
	!write ( 1, '(a)' ) 'SCALARS ae float'
	!write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	!write ( 1, '(f22.14)' ) ((ae(i,j),i=1,(nx)),j=1,(ny))

	!write ( 1, '(a)' ) 'SCALARS aw float'
	!write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	!write ( 1, '(f22.14)' ) ((aw(i,j),i=1,(nx)),j=1,(ny))
	
	!write ( 1, '(a)' ) 'SCALARS an float'
	!write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	!write ( 1, '(f22.14)' ) ((an(i,j),i=1,(nx)),j=1,(ny))
	
	!write ( 1, '(a)' ) 'SCALARS as float'
	!write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	!write ( 1, '(f22.14)' ) ((as(i,j),i=1,(nx)),j=1,(ny))
	
	!write ( 1, '(a)' ) 'SCALARS ap float'
	!write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	!write ( 1, '(f22.14)' ) ((ap(i,j),i=1,(nx)),j=1,(ny))
	
	!write ( 1, '(a)' ) 'SCALARS b float'
	!write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	!write ( 1, '(f22.14)' ) ((b(i,j),i=1,(nx)),j=1,(ny))
	
	write ( 1, '(a)' ) 'SCALARS T float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((T(i,j),i=1,(nx)),j=1,(ny))

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
	write(*,*) 'Iteracao no Jacobi:',iter	
	write(*,*)
	write(*,*) 'Global error of T:',Res_T,"  |  ",'MAX local error of T:',Rmax_T
	return

end Subroutine print_Residual
 
end module visual
