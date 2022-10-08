module user

	implicit none

	INTEGER, PARAMETER :: 			&
		DP = SELECTED_REAL_KIND(14),	&

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
		Nx=81,				& ! Nodes in the x-direction
		Ny=81,				& ! Nodes in the y-directionc
		!Ny=69 to have Ly=12h and Ny 91 to have Ly=18h

		ALGORITMO=1	!**    (ALGORITMO = SIMPLE (1), SIMPLEC (2))

	real(KIND = DP)  ::			&

		Hx,Hy,					& ! Domain size
		Factorx=1.0D+0,		& ! x-Screcth grid factor
		Factory=1.0D+0,		& ! y-Screcth grid factor
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!
		epsilonU=1.0D-9,		& ! Residual tolerance for U
		epsilonV=1.0D-9,		& ! Residual tolerance for V
		epsilonP=1.0D-8,		& ! Residual tolerance for P
		epsilonT=1.0D-9,		& ! Residual tolerance for T
		epsilonTK=1.0D-9,		& ! Residual tolerance for TK
		epsilone=1.0D-9,		& ! Residual tolerance for e

		relaxU=0.4D+00,			& ! Relax factor for U
		relaxV=0.4D+00,			& ! Relax factor for V
		relaxP=1.0D+00,			& ! Relax factor for P
		relaxT=0.3D+00,			& ! Relax factor for T
		relaxTK=0.2D+00,		& ! Relax factor for TK
		relaxe=0.2D+00,			& ! Relax factor for e

		dt =5.0D-02,			& ! Time step
		to=0.0D+00			 ! Initial time

!-----------------------------------------------------------------------------------------------!
!					RANS PARAMATERS						!
!-----------------------------------------------------------------------------------------------!
	integer::															&
		rans_model = 0 !  0:Laminar |1:k-epsilon HH |2:k-epsilon JL |3: k-epsilon LS | 4:k-w PHD | 5:k-w WX 1994

	real(KIND = DP)  ::													&
		Cmu,C1e,C2e,Theta_k,Theta_e,Theta_t,C_k,C3e,Crw,Cw1,Cw2,		&
!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!

		Ra=1.0D+07,			  & ! Rayleigh Number
		Bbeta=3.322D-03,		  & ! Expansion Coefficient for the Boussinesq Approximation
		Pr=0.71D+0,			  & ! Prantl Number
		Th=300D+0,			  & ! Hot Temperature
		Tc=288D+00,			  & ! Cold Temperature
		G=9.8D+00,			  & ! Gravity
		R=287.0D+00,			  & ! Gas Constant for Air
		Po=101325.0D+00,		  & ! Initial Pressure			

		CONTERo,							& ! Thermal conductivity
		MUo,ratio_mu,ratio_conter,			& ! Viscoity
		RHOo,								& ! Density
		Cpo	,								& ! Specific heat	
		Res_U,Res_V,Res_P,Res_T,Res_TK,Res_e,Rmax_U,Rmax_V,Rmax_P,Rmax_T,Rmax_e,Rmax_TK, &
		tadim,DeltaT,dt_adim,Tref,PSIP,RHOmed,Pmed, &
		time_start,time_finish

	integer::				&
		case_jet,			&
		Boussinesq=1,	  	&	! If Boussinesq=1 applies Bouss. Approx.. If Boussinesq=2 DO NOT APPLY Bouss. Appr. and density is variable
		CttProperties=1,	&  	! If CttProperties=1 properties are constant. if CttProperties=2 mu and conter are tmeprature dependent by Sutherland's law	
		ite_coupling=0
!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						!
!-----------------------------------------------------------------------------------------------!
	
	integer::					&
!		itermax_t=5000,			& ! Maximum number of iterations

		itermax=5000,			& ! Maximum number of iterations
		iter,iter_t,i,j,		&

		npas_P=15,			& ! Time to apply the P-solver subroutine per iteration 
		npas_U=1,			& ! Time to apply the U-solver subroutine per iteration 
		npas_V=1,			& ! Time to apply the V-solver subroutine per iteration
		npas_T=1,			& ! Time to apply the T-solver subroutine per iteration 

		n=2,				& ! Interpolation scheme: 1:Central, 2:Upwind, 3:Hybrid, 4:Power law, 5:Exponential


		n_high_order_scheme = 0 	! No HOS: 		n_high_order_scheme = 0	
									! Van Leer: 	n_high_order_scheme = 1
									! Van Albada: 	n_high_order_scheme = 2
									! Min-Mod: 		n_high_order_scheme = 3
									! SUPERBEE : 	n_high_order_scheme = 4
									! QUICK :		n_high_order_scheme = 5
									! UMIST :		n_high_order_scheme = 6

	logical :: high_order_on_velocities_u=.false., high_order_on_velocities_v=.false., high_order_on_scalars=.false.


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
