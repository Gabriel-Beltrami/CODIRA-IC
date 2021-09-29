module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
		Nx=101,				& ! Nodes in the x-direction
		Ny=101				! Nodes in the y-directionc
		!Ny=69 to have Ly=12h and Ny 91 to have Ly=18h


	real(KIND = DP)  ::			&

		Hx = 1, &
		Hy = 1,					& ! Domain size
		Factorx=1.0D-010,		& ! x-Screcth grid factor
		Factory=1.0D-010,		& ! y-Screcth grid factor
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!

		epsilonT=1.0D-10,		& ! Residual tolerance for T

		relaxT=0.5D+00,			& ! Relax factor for T

		dt =2.0D-01,			& ! Time step
		to = 0,			 &! Initial time

!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!
                Ra=1.0D+10,			&
		MUo,				& ! Viscoity
		RHOO = 2702,				& ! Density
		Pr=0.71D+0,			& ! Prantl number
		Th=100D+0,			& ! Hot temperature
		Tc=0D+00,			& ! Cold temperature
		g=9.8D+00,			& ! gravity
		CONTERo = 237,			& ! Thermal conductivity
		Cpo	= 903,		& ! Specific heat
		ts=0.0D+0,			&


		U_const = 2, &
		V_const = 0, &

		T_ref = 50, &
		T_S = 0, &
		T_N = 100, &
		T_E = 0, &
		T_W = 100, &

		Res_T, Rmax_T, &
		tadim,DeltaT,dt_adim,Tref,PSIP,RHOmed,Pmed, &
		time_start,time_finish








!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						!
!-----------------------------------------------------------------------------------------------!

	integer::					&
!		itermax_t=5,			& ! Maximum number of iterations

		itermax=100000000,			& ! Maximum number of iterations
		iter,iter_t,i,j,		&

		n=2 ! Interpolation scheme: 1:Central, 2:Upwind, 3:Hybrid, 4:Power law, 5:Exponential


!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					!
!-----------------------------------------------------------------------------------------------!
	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X,DXU,XU

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y,DYV,YV,U_in,Tk_in,e_in,T_in

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		apu,apv,ae,aw,an,as,ap,b,du,dv,	&
		Ux,U,Uold,			&
		Vx,V,Vold,			&
		Px,P,Pc,T,Told,			&
		mu,CONTER,cp,rho,RHOold,	&
		TxxL, TyyL, TxyL,		&
		TxxT, TyyT, TxyT,		&
		MuT, TK,Tkx, e,			&
		Fep,Fnp,Pk,Gk,APUNB,APVNB,Rij_T,Rij_u,Rij_v,Sdc_ij_v, &
		P_rad

end module user
