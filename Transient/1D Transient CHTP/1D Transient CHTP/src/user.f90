module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!				GRID, DOMAIN SIZE AND SOME NUMERICAL PARAMETERS       	    !
!-----------------------------------------------------------------------------------------------!
		
		Nx=7,				& ! Nodes in the x-direction
		itermax=100			  ! Maximum number of iterations
	
	real(KIND = DP), PARAMETER ::		&
		
		tmax=1.0D-01,               	& ! Maximum time
		dt =1.0D-02			  ! Time step
	
	real(KIND = DP)  ::			&

		Hx=1.0D+00,			& ! Domain size [m]
		DX,				&
		
!-----------------------------------------------------------------------------------------------!
!					OTHER NUMERICAL PARAMETERS				    !
!-----------------------------------------------------------------------------------------------!

		epsilonT=1.0D-10,		& ! Residual tolerance for T
		
		to=0.0D+00,			& ! Initial time

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				    !
!-----------------------------------------------------------------------------------------------!

		Ta=0D+00,			& ! Temperature at A [ºC]
		Tb=0D+00,			& ! Temperature at B [ºC]
		tin=20D+00,			& ! Temperature at 0 [ºC]
		Lambda=2.0D+00,		& ! Thermal conductivity [W/m.ºC]
		Rho=1.0D+00,			& ! Density [kg/m³]
		Cp=1.0D+00,			& ! Specific heat [J/kg.ºC]
		
		Res_T,Rmax_T, &
		time_start,time_finish

!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						   !
!-----------------------------------------------------------------------------------------------!
	
	integer::				&

		iter,i,iter_t
	
	integer, PARAMETER ::			&
	
		itermax_t=int(tmax/dt+2.0D+00)

!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					   !
!-----------------------------------------------------------------------------------------------!
	
	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X, 	&
		ae,aw,apo,	&
		ap,b	
		
	real(KIND = DP), dimension(Nx,itermax,itermax_t) :: 	&
		T			

	real(KIND = DP), parameter :: PI = 4*ATAN(1.0D+00)

end module user
