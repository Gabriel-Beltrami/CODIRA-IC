!--------------------------------------------------------------------------------------------------
!     MODULE   Grids
!>    @ingroup CommComb_tabulation_grp
!!    @brief   Filling.
!!    @authors Jan Mateu
!!    @date    22/08/2017
!--------------------------------------------------------------------------------------------------
module grids
	use user
	implicit none

contains
!--------------------------------------------------------------------------------------------------
!     SUBROUTINE Grid 2D
!>    @brief     Filling.
!!    @authors   Jan Mateu
!!    @date      22/08/2017
!!    @param     [in] In Filling
!!    @param     [out] Out Filling
!--------------------------------------------------------------------------------------------------
subroutine GRID2d (DXP,DYP,DXU,DYV,X,Y,XU,YV,Hx,Hy)

	INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

	real(KIND = DP)  :: TETA1,TETA2,ALFA1,FUNCION1,FUNCION2,Hy,Hx

	real(KIND = DP), dimension(Nx) :: DXP,X,DXU,XU

	real(KIND = DP), dimension(Ny) :: DYP,Y,DYV,YV

	integer :: j,i,NLX,NLY,AI


!-----------------------------------------------------------------------------------------------!
!				COUNTING TIME STARTING						!
!-----------------------------------------------------------------------------------------------!

!	TIMECODE=TIMEF()

!-----------------------------------------------------------------------------------------------!
!				ADJUST FACTOR BETWEEN NODES					!
!-----------------------------------------------------------------------------------------------!

	ALFA1=2.0D+00*FACTORX

!-----------------------------------------------------------------------------------------------!
!				NUMBER OF LINES BETWEEN CV					!
!-----------------------------------------------------------------------------------------------!

	NLX=NX-1
	NLY=NY-1

      DO AI=1,NLX

            TETA1=ALFA1*(((dble(AI-1))/dble(NLX-1))-0.5D+00)
            TETA2=(ALFA1/2D+00)
            FUNCION1=DTANH(TETA1)
            FUNCION2=DTANH(TETA2)

        XU(AI)=(HX/2.0D+00)*(dble(1)+(FUNCION1/FUNCION2))

      END DO

      ALFA1=2.0D+00*FACTORY

      DO AI=1,NLY

            TETA1=ALFA1*((dble(AI-1)/dble(NLY-1))-0.5D+00)
            TETA2=(ALFA1/2.0D+00)
            FUNCION1=DTANH(TETA1)
            FUNCION2=DTANH(TETA2)

        YV(AI)=(HY/2.0D+00)*(1+(FUNCION1/FUNCION2))
      END DO

      X(1)=XU(1)
      X(NX)=XU(NX-1)

      DO I=2,NX-1
                  X(I)=((XU(I-1)+XU(I))/2.0D+00)
      END DO

      Y(1)=YV(1)
      Y(NY)=YV(NY-1)

      DO J=2,NY-1
                  Y(J)=((YV(J-1)+YV(J))/2.0D+00)
      END DO

      DXP(1)=0.0D+00
      DXP(NX)=DXP(1)

      DO I=2,NX-1
                 DXP(I)=XU(I)-XU(I-1)
      END DO

      DYP(1)=0.0D+00
      DYP(NY)=DYP(1)

      DO J=2,NY-1
                  DYP(J)=YV(J)-YV(J-1)
      END DO

      DXU(1)= DXP(2)/2.0D+00
      DXU(NX-1)=DXP(NX-1)/2.0D+00

      DO I=2,NX-2
                 DXU(I)=X(I+1)-X(I)
      END DO


      DYV(1)= DYP(2)/2.0D+00
      DYV(NY-1)=DYP(NY-1)/2.0D+00

      DO J=2,NY-2
                 DYV(J)=Y(J+1)-Y(J)
      END DO

	return
end subroutine GRID2d

end module grids
