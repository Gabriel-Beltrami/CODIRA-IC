module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					    !
!-----------------------------------------------------------------------------------------------!
		
		n_points=101			  ! No_of_points
	
	real(KIND = DP)  ::			&

		dom_length=1.0D+00,		& ! Domain size [m]
		h,				&
		
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					   !
!-----------------------------------------------------------------------------------------------!

		error_req=1.0D-07, 		& ! Final required error residual

		alpha=8.0D-01,			& ! Under-relaxation factor for u,v
		alpha_p=8.0D-01,		& ! Under-relaxation factor for p

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				    !
!-----------------------------------------------------------------------------------------------!

		Re=100.0D+00,			& ! Reynolds number
		nu,				& ! Viscosity
		u_E,u_W,			&
		v_N,v_S,			&
		a_E1,a_E2,a_E3,		& ! A_e,a_E,a_e
		a_W,a_S,			&
		a_N1,a_N2,a_N3,		& ! A_n,a_N,a_n
		a_P,				&
				
		error,time_start,time_finish

!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						   !
!-----------------------------------------------------------------------------------------------!
	
	integer::				&

		iterations,i,j			

!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					    !
!-----------------------------------------------------------------------------------------------!
		
	real(KIND = DP), dimension(n_points) ::	&
		x,y
	
	real(KIND = DP), dimension(n_points,n_points) ::	&
		u_final,v_final,p_final	! Initializing the variables: Final collocated variables
		
	real(KIND = DP), dimension(n_points+1,n_points) :: 	&
		u,u_star,d_e,u_new		! Initializing the variables: Staggered variables

	real(KIND = DP), dimension(n_points,n_points+1) :: 	&
		v,v_star,d_n,v_new		! Initializing the variables: Staggered variables
	
	real(KIND = DP), dimension(n_points+1,n_points+1) ::	&
		p,p_star,pc,b,p_new		! Initializing the variables: Staggered variables

end module user
