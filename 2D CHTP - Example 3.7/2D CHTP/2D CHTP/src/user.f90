module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
		Nx=51,			&	 ! Nodes in the x-direction
		Ny=51,			&	 ! Nodes in the y-direction
		itermax=200000,	&	 ! Maximum number of iterations
		k=(Nx+1)/2
		  
	real(KIND = DP)  ::			&

		Hx=1.0D+00,			& ! Domain size [m]
		Hy=5.0D-01,			& ! Domain size [m]
		DX,DY,				&
		Hx1=5.0D-01,			&
		
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!

		epsilonT=1.0D-10,		& ! Residual tolerance for T

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!

		Ta=600D+00,			& ! Temperature at A [ºC]
		Tb=100D+00,			& ! Temperature at B [ºC]
		Q=0.0D+00,			& ! Imposed Heat flux [W]
		Lambda1=6.0D-02,		& ! Thermal conductivity #1 [W/m.ºC]
		Cp1=1D+00,			& ! Specific heat #1 [J/kg.ºC]
		Lambda2=1.0D-03,		& ! Thermal conductivity #2 [W/m.ºC]
		Cp2=1D+00,			& ! Specific heat #2 [J/kg.ºC]
		
		Res_T,Rmax_T, &
		Tref, &
		time_start,time_finish

!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						!
!-----------------------------------------------------------------------------------------------!
	
	integer::				&

		iter,i,j

!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					!
!-----------------------------------------------------------------------------------------------!
	
	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X
	
	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y
	
	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b			
		
	real(KIND = DP), dimension(Nx,Ny) :: 	&
		T,Ti			

end module user
