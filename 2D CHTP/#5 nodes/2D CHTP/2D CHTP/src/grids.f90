
module grids
	use user
	implicit none

contains

subroutine GRID2d (DXP,DYP,DX,DY,X,Y,Hx,Hy)

	INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

	real(KIND = DP)  :: Hx,Hy,DX,DY

	real(KIND = DP), dimension(Nx) :: DXP,X
	
	real(KIND = DP), dimension(Ny) :: DYP,Y

	integer :: i,j

!-----------------------------------------------------------------------------------------------!
!				NUMBER OF LINES BETWEEN CV					!
!-----------------------------------------------------------------------------------------------!

      X(1)=0.0D+00
      X(Nx)=Hx
      
      DX=Hx/dble(Nx-2)

      DO I=2,NX-1
                  X(I)=dble(I-2)*DX+DX/2.0D+00
      END DO

      DXP(1)=0.0D+00
      DXP(Nx)=1.0D+00/6.0D+00

      DO I=2,NX-1
                 DXP(I)=dble(X(I)-X(I-1))
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      Y(1)=0.0D+00
      Y(Ny)=Hy
      
      DY=Hy/dble(Ny-2)

      DO J=2,NY-1
                  Y(J)=dble(J-2)*DY+DY/2.0D+00
      END DO

      DYP(1)=0.0D+00
      DYP(Ny)=1.0D+00/6.0D+00

      DO J=2,NY-1
                 DYP(J)=dble(Y(J)-Y(J-1))
      END DO

      return
end subroutine GRID2d

end module grids
