module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
		Nx=7,			&	 ! Nodes in the x-direction
		itermax=2000			 ! Maximum number of iterations
		
		  
	real(KIND = DP)  ::			&

		Hx=1.0D+00,			& ! Domain size [m]
		DX,				&
		
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!

		epsilonT=1.0D-10,		& ! Residual tolerance for T

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!

		Ta=0D+00,			& ! Temperature at A [ºC]
		Tb=100D+00,			& ! Temperature at B [ºC]
		Lambda=2.0D+00,		& ! Thermal conductivity [W/m.ºC]
		Cp=1.0D+00,			& ! Specific heat [J/kg.ºC]
		Rho=1.0D+00,			& ! Density [kg/m³]
		u=50.0D+00,			& ! x-direction speed [m/s]
		
		Res_T,Rmax_T, &
		Tref, &
		time_start,time_finish

!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						!
!-----------------------------------------------------------------------------------------------!
	
	integer::				&

		iter,i,			&

		n=2				! Interpolation scheme: 1:Central, 2:Upwind, 3:Hybrid, 4:Power law, 5:Exponential
!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					!
!-----------------------------------------------------------------------------------------------!
	
	real(KIND = DP) ::			&
		Fep
	
	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X, 	&
		ae,aw,ap,b			
		
	real(KIND = DP), dimension(Nx) :: 	&
		T,Ti			

end module user
