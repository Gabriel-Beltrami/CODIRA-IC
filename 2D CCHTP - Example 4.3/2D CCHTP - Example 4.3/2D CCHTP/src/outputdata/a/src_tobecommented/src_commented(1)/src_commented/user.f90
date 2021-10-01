!--------------------------------------------------------------------------------------------------
!     MODULE    User
!>    @ingroup  RANS
!!    @brief    Variables to be specified by the user, change in every set up problem.
!!    @authors	Jan Mateu
!!    @date     22/08/2017
!!    @param    [in] Filling.
!!    @param    [out] Filling.
!--------------------------------------------------------------------------------------------------
module user

implicit none

	INTEGER, DIMENSION(2) ::						&
		Z

	DOUBLE PRECISION, DIMENSION(35) ::	&
		ZZ

	INTEGER, PARAMETER ::								&
		DP = SELECTED_REAL_KIND(14)

!-----------------------------------------------------------------------------------------------!
!					GRID AND DOMAIN SIZE					!
!-----------------------------------------------------------------------------------------------!
	INTEGER ::                              &
		Nx=61,				& ! Nodes in the x-direction
		Ny=61,				& ! Nodes in the y-directionc

		ALGORITMO=2	!**    (ALGORITMO = SIMPLE (1), SIMPLEC (2))

	real(KIND = DP)  ::			&
!		Hx=0.045841D+00,		& ! Domain width
!		Hy=0.045841D+00,		& ! Domain higth
		Hx,Hy,				&
		Factorx=4.0D+00,		& ! x-Screcth grid factor
		Factory=3.0D+00,		& ! y-Screcth grid factor
!-----------------------------------------------------------------------------------------------!
!					NUMERICAL PARAMETERS					!
!-----------------------------------------------------------------------------------------------!
		epsilonU=1.0D-10,		& ! Residual tolerance for U
		epsilonV=1.0D-10,		& ! Residual tolerance for V
		epsilonP=1.0D-10,		& ! Residual tolerance for P
		epsilonT=1.0D-10,		& ! Residual tolerance for T

		relaxU=0.9D+00,		& ! Relax factor for U
		relaxV=0.9D+00,		& ! Relax factor for V
		relaxP=1.0D+00,			& ! Relax factor for P
		relaxT=1.0D+00,			& ! Relax factor for T

		dt =2.0D-01,			& ! Time step
		to=0.0D+00			 ! Initial time
!-----------------------------------------------------------------------------------------------!
!					RANS PARAMATERS						!
!-----------------------------------------------------------------------------------------------!
	integer::				&
		rans_model = 4 !  1:k-epsilon HH |2:k-epsilon JL |3: k-epsilon LS | 4:k-w PHD | 5:k-w WX 1994

	real(KIND = DP)  ::			&
		Cmu,				&
		C1e,				&
		C2e,				&
		Theta_k,			&
		Theta_e,			&
		Theta_t,			&
		C_k,				&
		C3e,Crw,Cw1,Cw2,				&

		relaxTK=0.4D+00,	& ! Relax factor for TK
		relaxe=0.4D+00,		& ! Relax factor for e

		epsilonTK=1.0D-10,		& ! Residual tolerance for TK
		epsilone=1.0D-10,		& ! Residual tolerance for e
!-----------------------------------------------------------------------------------------------!
!					THERMOPHYSICAL PROPERTIES				!
!-----------------------------------------------------------------------------------------------!

		Ra=1.0D+10,			&
		Bbeta=3.322D-03,		& ! Expansion coefficient for the Boussinesq approximation
		MUo,				& ! Viscoity
		RHOo,				& ! Density
		Pr=0.71D+0,			& ! Prantl number
		Th=300D+0,			& ! Hot temperature
		Tc=288D+00,			& ! Cold temperature
		g=9.8D+00,			& ! gravity
		R=287.0D+00,			& ! gas constant for air
		CONTERo,			& ! Thermal conductivity
		Po=101325.0D+00,		& ! Initial Pressure
		Cpo	,			& ! Specific heat
		ts=0.0D+0,			&

		Res_U,Res_V,Res_P,Res_T,Res_TK,Res_e,Rmax_U,Rmax_V,Rmax_P,Rmax_T,Rmax_e,Rmax_TK, &
		tadim,DeltaT,dt_adim,Tref,PSIP,RHOmed,Pmed, &
		time_start,time_finish

	integer::				&

		Boussinesq=1,	  & ! If Boussinesq=1 applies Bouss. Approx.. If Boussinesq=2 DO NOT APPLY Bouss. Appr. and density is variable
		CttProperties=1	    ! If CttProperties=1 properties are constant. if CttProperties=2 mu and conter are tmeprature dependent by Sutherland's law

!-----------------------------------------------------------------------------------------------!
!					SOLVER OPTIONS						!
!-----------------------------------------------------------------------------------------------!

	integer::				&
!		itermax_t=5000,			& ! Maximum number of iterations
		itermax=800000,			& ! Maximum number of iterations
		iter,iter_t,i,j,		&

		npas_P=10,			& ! Time to apply the P-solver subroutine per iteration
		npas_U=1,			& ! Time to apply the U-solver subroutine per iteration
		npas_V=1,			& ! Time to apply the V-solver subroutine per iteration
		npas_T=1,			& ! Time to apply the T-solver subroutine per iteration

		n=2				! Interpolation scheme: 1:Central, 2:Upwind, 3:Hybrid, 4:Power law, 5:Exponential

!-----------------------------------------------------------------------------------------------!
!					MATRIX DIMENSIONS					!
!-----------------------------------------------------------------------------------------------!
	real(KIND = DP), dimension(:), allocatable ::		&
		DXP,X,DXU,XU

	real(KIND = DP), dimension(:), allocatable ::		&
		DYP,Y,DYV,YV

	real(KIND = DP), dimension(:,:), allocatable ::		&
		apu,apv,ae,aw,an,as,ap,b,du,dv,	&
		Ux,U,Uold,			&
		Vx,V,Vold,			&
		Px,P,Pc,T,Told,			&
		mu,CONTER,cp,rho,RHOold,	&
		TxxL, TyyL, TxyL,		&
		TxxT, TyyT, TxyT,		&
		MuT, TK,Tkx, e,			&
		Fep,Fnp,Pk,Gk,APUNB,APVNB

	CONTAINS

	subroutine readinterface()

	implicit none

		INTEGER :: 																&
			ios,nn,ii,XX,YY,GG,k,q,vv,h

		INTEGER, PARAMETER ::											&
			read_unit = 99

		CHARACTER(len=200), ALLOCATABLE ::				&
			command(:)

		CHARACTER(len=200) ::				&
			line

		DOUBLE PRECISION :: 				&
			c

		CHARACTER(len=500) ::				&
			sig,a,char_a,char_b,char_c

	  OPEN(UNIT=read_unit, FILE='interface.txt', IOSTAT=ios)
			IF (ios/=0) STOP "Error opening file data.dat"
				nn=0
			DO
			READ(read_unit, '(A)', IOSTAT=ios) line
				IF (ios /= 0) EXIT
				nn = nn + 1
			END DO
			ALLOCATE(command(nn))
			REWIND(read_unit)
			DO ii = 1,nn
				read(read_unit, '(A)') command(ii)
			END DO
		CLOSE(read_unit)

		h=1

		DO ii=1,nn
			sig=command(ii)
			XX=1
			YY=1
			GG=0
			q=0
			vv=0
			DO WHILE (GG/=1 .AND. YY<=LEN(sig))
				IF (sig(YY:YY)=='=') THEN
					vv=1
					YY=YY+1
					a=sig(YY:YY)
					YY=YY+1
					DO WHILE (sig(YY:YY)/=',')
						q=1
						char_a=a
						char_b=sig(YY:YY)
						char_c=trim(adjustl(char_a))//trim(adjustl(char_b))
						a=char_c
						YY=YY+1
					END DO
					IF (q/=1) THEN
						READ(unit=a,fmt=*)c
					END IF
					GG=1
				END IF
				YY=YY+1
			END DO
			IF (q==1) THEN
				READ(unit=char_c,fmt=*)c
			END IF
			IF (vv==1) THEN
				IF (h<3) THEN
					Z(h)=c
					h=h+1
				ELSE
					ZZ(h-2)=c
					h=h+1
			  END IF
			END IF
		END DO

		ALLOCATE(DXP(Z(1)))
		ALLOCATE(X(Z(1)))
		ALLOCATE(DXU(Z(1)))
		ALLOCATE(XU(Z(1)))
		ALLOCATE(DYP(Z(2)))
		ALLOCATE(Y(Z(2)))
		ALLOCATE(DYV(Z(2)))
		ALLOCATE(YV(Z(2)))
		ALLOCATE(apu(Z(1),Z(2)))
		ALLOCATE(apv(Z(1),Z(2)))
		ALLOCATE(ae(Z(1),Z(2)))
		ALLOCATE(aw(Z(1),Z(2)))
		ALLOCATE(an(Z(1),Z(2)))
		ALLOCATE(as(Z(1),Z(2)))
		ALLOCATE(ap(Z(1),Z(2)))
		ALLOCATE(b(Z(1),Z(2)))
		ALLOCATE(du(Z(1),Z(2)))
		ALLOCATE(dv(Z(1),Z(2)))
		ALLOCATE(Ux(Z(1),Z(2)))
		ALLOCATE(U(Z(1),Z(2)))
		ALLOCATE(Uold(Z(1),Z(2)))
		ALLOCATE(Vx(Z(1),Z(2)))
		ALLOCATE(V(Z(1),Z(2)))
		ALLOCATE(Vold(Z(1),Z(2)))
		ALLOCATE(Px(Z(1),Z(2)))
		ALLOCATE(P(Z(1),Z(2)))
		ALLOCATE(Pc(Z(1),Z(2)))
		ALLOCATE(T(Z(1),Z(2)))
		ALLOCATE(Told(Z(1),Z(2)))
		ALLOCATE(mu(Z(1),Z(2)))
		ALLOCATE(CONTER(Z(1),Z(2)))
		ALLOCATE(cp(Z(1),Z(2)))
		ALLOCATE(rho(Z(1),Z(2)))
		ALLOCATE(RHOold(Z(1),Z(2)))
		ALLOCATE(TxxL(Z(1),Z(2)))
		ALLOCATE(TyyL(Z(1),Z(2)))
		ALLOCATE(TxyL(Z(1),Z(2)))
		ALLOCATE(TxxT(Z(1),Z(2)))
		ALLOCATE(TyyT(Z(1),Z(2)))
		ALLOCATE(TxyT(Z(1),Z(2)))
		ALLOCATE(MuT(Z(1),Z(2)))
		ALLOCATE(TK(Z(1),Z(2)))
		ALLOCATE(Tkx(Z(1),Z(2)))
		ALLOCATE(e(Z(1),Z(2)))
		ALLOCATE(Fep(Z(1),Z(2)))
		ALLOCATE(Fnp(Z(1),Z(2)))
		ALLOCATE(Pk(Z(1),Z(2)))
		ALLOCATE(Gk(Z(1),Z(2)))
		ALLOCATE(APUNB(Z(1),Z(2)))
		ALLOCATE(APVNB(Z(1),Z(2)))

		Nx=Z(1)
		Ny=Z(2)
		ALGORITMO=ZZ(1)
		Factorx=ZZ(2)
		Factory=ZZ(3)
		EpsilonU=ZZ(4)
		EpsilonV=ZZ(5)
		EpsilonP=ZZ(6)
		EpsilonT=ZZ(7)
		RelaxU=ZZ(8)
		RelaxV=ZZ(9)
		RelaxP=ZZ(10)
		RelaxT=ZZ(11)
		Dt=ZZ(12)
		To=ZZ(13)
		Rans_model = ZZ(14)
		RelaxTK=ZZ(15)
		Relaxe=ZZ(16)
		EpsilonTK=ZZ(17)
		Epsilone=ZZ(18)
		Ra=ZZ(19)
		Bbeta=ZZ(20)
		Pr=ZZ(21)
		Th=ZZ(22)
		Tc=ZZ(23)
		G=ZZ(24)
		R=ZZ(25)
		Po=ZZ(26)
		Ts=ZZ(27)
		Boussinesq=ZZ(28)
		CttProperties=ZZ(29)
		Itermax=ZZ(30)
		Npas_P=ZZ(31)
		Npas_U=ZZ(32)
		Npas_V=ZZ(33)
		Npas_T=ZZ(34)
		N=ZZ(35)

		CONTAINS

		FUNCTION is_numeric(string)
		IMPLICIT NONE

			CHARACTER(len=*), INTENT(IN) ::				&
				string

			LOGICAL ::				&
				is_numeric

			REAL ::					  &
				x

			INTEGER ::				&
				e

			READ(string,*,IOSTAT=e) x
			is_numeric = e == 0
		END FUNCTION is_numeric

		INTEGER FUNCTION logic2dbl(a)
		IMPLICIT NONE

			LOGICAL, INTENT(in) ::				&
				 a

			IF (a) THEN
				logic2dbl = 1
			ELSE
				logic2dbl = 0
			END IF
		END FUNCTION logic2dbl

		INTEGER FUNCTION numlines()
		IMPLICIT NONE

			INTEGER ::				&
				nlines

			nlines = 0

			OPEN(UNIT=1,FILE='interface.txt',FORM="FORMATTED",STATUS="OLD",ACTION="READ")
				DO
					READ (1,*, END=10)
					nlines = nlines + 1
				END DO
			10 CLOSE (1)

			numlines=nlines
		end function numlines
	end subroutine readinterface
end module user
