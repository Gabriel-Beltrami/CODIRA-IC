module visual
	use user
	use grids
	implicit none
	
contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!                            Varaible Fields Data Writing                                       !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine write_vtk(u_final,v_final)

	implicit none

	character(len = 27) :: str
	integer :: i,j
	real(KIND = DP), dimension(n_points) :: x_print,y_print
	real(KIND = DP), dimension(n_points,n_points) :: u_final,v_final

	do i = 1,n_points
		x_print(i)=x(i)
	end do

	do j = 1,n_points
		y_print(j)=y(j)
	end do

	str = 'outputdata/output.vtk'
	open(unit=1, file = str)

	write ( 1,'(a)' )'# vtk DataFile Version 3.0'
	write ( 1,'(a)' ) 'Uniform Rectilinear - Rectilinear Grid'
	write ( 1,'(a)' ) 'ASCII'
	write ( 1,'(a)' ) 'DATASET RECTILINEAR_GRID'
	write ( 1,'(a,3i4)' )'DIMENSIONS', n_points,n_points,1
	write ( 1,'(a13,1i4,a6)' ) 'X_COORDINATES',n_points,' float'
	write ( 1,'(f22.14)' ) (x_print(i),i=1,n_points)
	write ( 1,'(a13,1i4,a6)' ) 'Y_COORDINATES',n_points,' float'
	write ( 1,'(f22.14)' ) (y_print(j),j=1,n_points)
	write ( 1,'(a13,1i4,a6)' ) 'Z_COORDINATES',1,' float'
	write ( 1,'(f22.14)' ) 1.0
	
	write ( 1, '(a, I16)' ) 'POINT_DATA ', (n_points)*(n_points)*1
	
	write ( 1, '(a)' ) 'SCALARS u_final float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((u_final(i,j),i=1,(n_points)),j=1,(n_points))

	write ( 1, '(a)' ) 'SCALARS v_final float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((v_final(i,j),i=1,(n_points)),j=1,(n_points))

	close(1)

end subroutine write_vtk
 
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					print_Residual						    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine print_Residual(iter,error)

	real(KIND = DP) :: error
	integer :: iter

	write(*,*) "______________________________________________"
	write(*,*) 'SIMPLE Iteration:',iter	
	write(*,*)
	write(*,*) 'Error:',error
	return

end Subroutine print_Residual

end module visual
