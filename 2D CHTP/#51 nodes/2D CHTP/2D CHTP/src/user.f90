module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
		Nx=51,			&	 ! Nodes in the x-direction
		Ny=51,			&	 ! Nodes in the y-direction
		itermax=80000			 ! Maximum number of iterations
		
		  
	real(KIND = DP)  ::			&

		Hx=1.0D+00,			& ! Domain size [m]
		Hy=1.0D+00,			& ! Domain size [m]
		DX,DY,				&
		
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!

		epsilonT=1.0D-10,		& ! Residual tolerance for T

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!

		Ta=0D+00,			& ! Temperature at A [ºC]
		Tb=100D+00,			& ! Temperature at B [ºC]
		Lambda=35D+00,			& ! Thermal conductivity [W/m.ºC]
		Cp=130D+00,			& ! Specific heat [J/kg.ºC]
		
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
		
	real(KIND = DP), dimension(Nx,Ny,itermax) :: 	&
		T			

end module user
