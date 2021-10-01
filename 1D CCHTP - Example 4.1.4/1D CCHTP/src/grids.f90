
module grids
	use user
	implicit none

contains

subroutine GRID1d (DXP,DX,X,Hx)

	INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

	real(KIND = DP)  :: Hx,DX

	real(KIND = DP), dimension(Nx) :: DXP,X

	integer :: i

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
      DXP(Nx)=DX/2.0D+00

      DO I=2,NX-1
                 DXP(I)=dble(X(I)-X(I-1))
      END DO

      return
end subroutine GRID1d

end module grids
