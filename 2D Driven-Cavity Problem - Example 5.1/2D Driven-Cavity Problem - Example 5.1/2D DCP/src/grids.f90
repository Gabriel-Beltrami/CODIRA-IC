
module grids
	use user
	implicit none

contains

subroutine GRID2d (x,y,h,dom_length)

	INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

	real(KIND = DP)  :: h,dom_length 

	real(KIND = DP), dimension(n_points) :: x,y 

	integer :: i,j

!-----------------------------------------------------------------------------------------------!
!				Defining the problem domain                                      !
!-----------------------------------------------------------------------------------------------!
      
      h=dom_length/dble(n_points-1)

      x(1)=0.0D+00
      x(n_points)=dom_length

      DO I=2,N_POINTS-1
                  X(I)=h+X(I-1)	!X domain span
      END DO
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      y(1)=0.0D+00
      y(n_points)=dom_length

      DO J=2,N_POINTS-1
                  Y(J)=h+Y(J-1)	!Y domain span
      END DO

      return
end subroutine GRID2d

end module grids
