module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
		Nx=101,			&	 ! Nodes in the x-direction
		Ny=101,			&	 ! Nodes in the y-direction
		itermax=20000			 ! Maximum number of iterations
		
		  
	real(KIND = DP)  ::			&

		Hx=1.0D+00,			& ! Domain size [m]
		Hy=1.0D+00,			& ! Domain size [m]
		DX,DY,				&
		
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!

		epsilonT=1.0D-10,		& ! Residual tolerance for T
		relaxT=0.5D+00,		& ! Relax factor for T

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!

		Ta=0D+00,			& ! Temperature at A [ºC]
		Tb=100D+00,			& ! Temperature at B [ºC]
		Lambda=237D+00,		& ! Thermal conductivity [W/m.ºC]
		Cp=903D+00,			& ! Specific heat [J/kg.ºC]
		Rho=2702D+00,			& ! Density [kg/m³]
		u=2.0D+00,			& ! x-direction speed [m/s]
		v=2.0D+00,			& ! y-direction speed [m/s]
		
		Res_T,Rmax_T, &
		Tref, &
		time_start,time_finish

!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						!
!-----------------------------------------------------------------------------------------------!
	
	integer::				&

		iter,i,j,			&

		n=2				! Interpolation scheme: 1:Central, 2:Upwind, 3:Hybrid, 4:Power law, 5:Exponential

!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					!
!-----------------------------------------------------------------------------------------------!
	
	real(KIND = DP) ::			&
		Fep,Fnp
		
	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X
	
	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y
	
	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b			
		
	real(KIND = DP), dimension(Nx,Ny) :: 	&
		T,Ti			

end module user
