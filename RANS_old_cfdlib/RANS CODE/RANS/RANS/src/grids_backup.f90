module grids
	use user
	implicit none
	
contains

!-----------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
subroutine GRID_half_jet(DXP,DYP,DXU,DYV,X,Y,XU,YV)
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------!

  implicit none

  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14) 

  integer :: n_core_y, variable_zone_y

  real(KIND = DP)  :: delta, delta_x, delta_y, raison_x, raison_y

  real(KIND = DP), dimension(Nx) :: DXP,X,DXU,XU
    
  real(KIND = DP), dimension(Ny) :: DYP,Y,DYV,YV,y_med
      
  integer :: j,i,NLX,NLY,AI


!    delta_x = h/5D+0  ! 101x61
!    raison_x = 1.0118D+0
!    delta_y = h/5.0D+0    
!    raison_y = 1.008D+0

   delta_x = h/20D+0  ! 141x81
   raison_x = 1.00D+0
   delta_y = h/3000.0D+0    
   raison_y = 1.05D+0

!    delta_x = h/11D+0  ! 181x101
!    raison_x = 1.0084D+0
!    delta_y = h/11.0D+0    
!    raison_y = 1.0096D+0   
!-----------------------------------------------------------------------------------------------!
!       NUMBER OF LINES BETWEEN CV          !
!-----------------------------------------------------------------------------------------------!

  NLX=NX-1   
  NLY=NY-1

  !n_core_y = NLY/1.22
  n_core_y = 1 


    YV(1) = 0.0D+0

!      ! 1st Region of the sponge layer in y direction 
!      delta = delta_y 
!      do j = 2, n_core_y

!        YV(j) = yv(j-1) + delta
 
!      enddo

    do j =  n_core_y + 1, NLY
        delta_y  = delta_y   * raison_y 
        YV(j) = YV(j-1) + delta_y 
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
  delta = delta_x
  XU(1) = 0.0D+0

      DO AI=2,NLX

        delta = delta * raison_x

        XU(AI)=XU(AI-1) + delta

      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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



      Hx = x(Nx)-x(1) 
      Hy = y(Ny)-y(1)

      write(*,*) "---------------------------------"
      write(*,*) "Hx = " , Hx
      write(*,*) "Hy = " , Hy 
      write(*,*) "---------------------------------"     
  return
end subroutine GRID_half_jet




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















!-----------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
subroutine GRID_jet(DXP,DYP,DXU,DYV,X,Y,XU,YV)
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------!



  implicit none

  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14) 

  integer :: n_core_y, variable_zone_y

  real(KIND = DP)  :: delta, delta_x, delta_y, raison_x, raison_y

  real(KIND = DP), dimension(Nx) :: DXP,X,DXU,XU
    
  real(KIND = DP), dimension(Ny) :: DYP,Y,DYV,YV,y_med
      
  integer :: j,i,NLX,NLY,AI



  delta_x = h/10.0D+0
  raison_x = 1.010D+0

  delta_y = h/15.0D+0
  raison_y = 1.03D+0



!-----------------------------------------------------------------------------------------------!
!       NUMBER OF LINES BETWEEN CV          !
!-----------------------------------------------------------------------------------------------!

  NLX=NX-1   
  NLY=NY-1

  n_core_y = NLY/3


                  delta = delta_y 
                  variable_zone_y = NLY/2 - n_core_y/2 
                  do j = n_core_y/2 + 1 , NLY/2
                      delta = delta  * raison_y
                  end do 

    YV(1) = 0.0D+0

    ! 1st Region of the sponge layer in y direction 

    do j = 2, variable_zone_y

      YV(j) = yv(j-1) + delta
      delta = delta  / raison_y      
    enddo

    delta = delta_y 

    do j = variable_zone_y + 1, variable_zone_y + n_core_y
        YV(j) = YV(j-1) + delta 
    enddo 

    ! 2nd Incrising zone in y direction 

    do j = variable_zone_y + n_core_y + 1, NLY
        delta = delta  * raison_y 
        YV(j) = YV(j-1) + delta 
    enddo


  do j=1,NLY
  y_med(j) = (0.0D+0 + yv(NLY) ) / 2.0D+0
  end do
  

   YV(:) = YV(:) - (y_med )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
  delta = delta_x
  XU(1) = 0.0D+0

      DO AI=2,NLX

        delta = delta * raison_x

        XU(AI)=XU(AI-1) + delta

      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





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

      write(*,*) "Hx = " , x(Nx)-x(1) 
      write(*,*) "Hy = " , y(Ny)-y(1) 

      Hx = x(Nx)-x(1) 
      Hy = y(Ny)-y(1)

      write(*,*) "Hx = " , Hx, x(Nx)-x(1) 
      write(*,*) "Hy = " , Hy , y(Ny)-y(1) 
  return
end subroutine GRID_jet

end module grids
