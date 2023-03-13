module cfdlib
	use user
	use rans_models
	implicit none

	
contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					Compute dependent variables				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine dependent_var(DeltaT,Hx,Hy,dt,Tref,CONTERo,Po,Muo,U_bg)
			

	implicit none
	integer:: i,j
	real(KIND = DP)  :: DeltaT,Tref,Hx,Hy,dt,CONTERo,Po,Muo,D_hidraulic,U_bg

	Tref=(Th+Tc)/2.0D+00
	DeltaT=Th-Tc
!	RHOo=Po/(R*Tref);
	RHOo=1.1774D+0
	Cpo = 1006D+0
	contero=0.02617D+0
	Muo=1.847D-5

!	Muo= 1.68D-5*( (Tref/273.0D+0)**(3.0D+0/2.0D+0) )*( 383.5D+0/(Tref+110.5) ) !SUTHERLAND's Law
!	CONTERo=Muo*401.8D+0/2.84D-1
!	Cpo=Pr*CONTERo/MUo

	if (Boussinesq.eq.1) then
!		Hx=( (Ra*CONTERo*MUo)/(BBETA*g*DeltaT*Cpo*RHOo*RHOo) )**(1.0D+00/3.0D+00)
		!Hx=(RA*(Muo/RHOo)*(Muo/RHOo)/(PR*G*BBETA*(DeltaT)))**(1.0D+00/3.0D+00)
		Hx=BD+2.0D+0*BW
	else
		Hx=( (Ra*MUo*MUo*Tref)/(Pr*g*RHOo*RHOo*DeltaT) )**(1.0D+00/3.0D+00)
	end if
!write(*,*) "Hx",Hx,"RHOo",RHOo,"Cpo",Cpo,"contero",contero,"Muo",Muo
!stop
	Hy=B1+ABD
!	dt=dt_adim*Hx*Hx*RHOo*Cpo/(CONTERo)
!	write(*,*) dt
!	stop

	U_bg = V10*cos(4D+00*atan(1D+00)/180D+00*(90D+00-zeta))-U10*cos(4D+00*atan(1D+00)/180D+00*zeta)	
	!write(*,*) U_bg
    return
end subroutine dependent_var



!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			stressLT(Mu,U,V,DXP,DYP,DYV,DXU,TxxL,TyyL,Txy)					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine stressLT(Mu,U,V,DXP,DYP,DYV,DXU,TxxL,TyyL,TxyL, TxxT, TyyT, TxyT, Tk, MuT, RHO)	

	implicit none
	integer:: i,j
	real(KIND = DP)  :: MUn,Mus,Mup
	real(KIND = DP), dimension(Nx,Ny), intent(in) :: Mu,U,V, Tk, MuT, RHO
	real(KIND = DP), dimension(Nx,Ny), intent(out) :: TxxL, TyyL, TxyL, TxxT, TyyT, TxyT
	real(KIND = DP), dimension(Nx), intent(in) :: DXP,DXU
	real(KIND = DP), dimension(Ny), intent(in) :: DYP,DYV

!-----------------------------------------------------------------------------------------------!
!				Laminar stresses						!
!-----------------------------------------------------------------------------------------------!

	Do j=2,Ny-1
		Do i=2,(Nx-1)

			TxxL(i,j)=2.0D+0 * mu(i,j) * (U(i,j)-U(i-1,j)) / DXP(i)
		end do
	end do

	Do j=2,Ny-1
		Do i=1,(Nx)

			TyyL(i,j)=2.0D+0 * mu(i,j) * (V(i,j)-V(i,j-1)) / DYP(j)
		end do
	end do
   
	Do j=2,(Ny-2)
		Do i=2,(Nx-2)
			Mun = flux(mu(i,j+1),mu(i+1,j+1),dxu(i),DXP(i+1)/2.0D+0 )
			Mus = flux(mu(i,j),mu(i+1,j),dxu(i),DXP(i+1)/2.0D+0 )
			Mup = flux(mus,mun,dyv(j),DYP(j+1)/2.0D+0 )
			TxyL(i,j)=MUp * ((U(i,j+1)-U(i,j)) / DYV(j) + (V(i+1,j)-V(i,j)) / DXU(i)) 
		end do
	end do


!-----------------------------------------------------------------------------------------------!
!				Laminar Boundary stresses 					!	

	Do j = 2, Ny-1
		! West Boundary
		TxxL(1,j)  = 2.0D+0 * mu(1,j)  * (U(2,j)   - U(1   ,j))    / DXP(2)
		! East Boundary
		TxxL(Nx,j) = 2.0D+0 * mu(Nx,j) * (U(Nx-1,j)- U(Nx-2,j))    / DXP(Nx-1)
	end do

	Do i=2,(Nx-1)
		! South boundary
		TyyL(i,1)=2.0D+0 * mu(i,1) * (V(i,2)-V(i,1)) / DYP(2)
		! North boundary
		TyyL(i,Ny)=2.0D+0 * mu(i,Ny) * (V(i,Ny-1)-V(i,Ny-2)) / DYP(Ny-1)
	end do


	Do i=2,(Nx-1)

		MUn= (mu(i,1)+mu(i+1,1)) /2.0D+0 ! Strictly speaking mu in i=1 and Nx shouldnt be interpolated
		TxyL(i,1)=MUn * ((U(i,2)-U(i,1)) / DYV(1) + (V(i+1,1)-V(i,1)) / DXU(i)) 

		MUn= (mu(i,Ny)+mu(i+1,Ny)) /2.0D+0 ! Strictly speaking mu in i=1 and Nx shouldnt be interpolated
		TxyL(i,Ny-1)=MUn * ((U(i,Ny)-U(i,Ny-1)) / DYV(Ny-1) + (V(i+1,Ny-1)-V(i,Ny-1)) / DXU(i))


	!write(*,*) TxyL(i,Ny-1),MUn ,(U(i,Ny)-U(i,Ny-1)) / DYV(Ny-1) , (V(i+1,Ny-1)-V(i,Ny-1)) / DXU(i)


	end do



	Do j=2,(Ny-1)
		MUn= (mu(1,j)+mu(1,j+1)) /2.0D+0 ! Strictly speaking mu in i=1 and Nx shouldnt be interpolated
		TxyL(1,j)=MUn * ((U(1,j+1)-U(1,j)) / DYV(j) + (V(2,j)-V(1,j)) / DXU(1)) 

		MUn= (mu(Nx,j)+mu(Nx,j+1)) /2.0D+0 ! Strictly speaking mu in i=1 and Nx shouldnt be interpolated
		TxyL(Nx-1,j)=MUn * ((U(Nx-1,j+1)-U(Nx-1,j)) / DYV(j) + (V(Nx,j)-V(Nx-1,j)) / DXU(Nx-1)) 


	end do



!-----------------------------------------------------------------------------------------------!
!				Turbulent stresses 						!
!-----------------------------------------------------------------------------------------------!
	Do j=2,Ny-1
		Do i=2,(Nx-1)

			TxxT(i,j) = (-2.0D+0/3.0D+0)*RHO(i,j)*TK(i,j) +	&
					 2.0D+0 * muT(i,j) * (U(i,j)-U(i-1,j)) / DXP(i)
		end do
	end do

	Do j=2,Ny-1
		Do i=2,(Nx-1)

			TyyT(i,j)=(-2.0D+0/3.0D+0)*RHO(i,j)*TK(i,j) +	&
					2.0D+0 * muT(i,j) * (V(i,j)-V(i,j-1)) / DYP(j)
		end do
	end do
   
	Do j=2,(Ny-2)
		Do i=2,(Nx-2)
			Mun = flux(mut(i,j+1),mut(i+1,j+1),dxu(i),DXP(i+1)/2.0D+0 )
			Mus = flux(mut(i,j),mut(i+1,j),dxu(i),DXP(i+1)/2.0D+0 )
			Mup = flux(mus,mun,dyv(j),DYP(j+1)/2.0D+0 )

			TxyT(i,j)=MUp * ((U(i,j+1)-U(i,j)) / DYV(j) + (V(i+1,j)-V(i,j)) / DXU(i)) 
		end do
	end do

!-----------------------------------------------------------------------------------------------!
!				Turbulent Boundary stresses 					!
!-----------------------------------------------------------------------------------------------!
	

	Do j = 2, Ny-1
		! West Boundary
		TxxT(1,j)  = (-2.0D+0/3.0D+0)*RHO(1,j)*TK(1,j) +	&			!Fluid Boundary	
			2.0D+0 * muT(1,j)  * (U(2,j)   - U(1   ,j))    / DXP(2)			!Fluid Boundary	
		! East Boundary
		TxxT(Nx,j) = (-2.0D+0/3.0D+0)*RHO(Nx,j)*TK(Nx,j) +	&			!Fluid Boundary	
				 2.0D+0 * muT(Nx,j) * (U(Nx-1,j)- U(Nx-2,j))    / DXP(Nx-1)	!Fluid Boundary	
		! West Boundary
!		TxxT(1,j) = 0.0D+0	!Solid Boundary
!		TyyT(1,j) = 0.0D+0	!Solid Boundary
		! East Boundary
!		TxxT(Nx,j) = 0.0D+0	!Solid Boundary
!		TyyT(Nx,j) = 0.0D+0	!Solid Boundary
	end do





	Do i=2,(Nx-1)
		! South boundary
		TyyT(i,1)= (-2.0D+0/3.0D+0)*RHO(i,1)*TK(i,1) +	&				!Fluid Boundary	
				2.0D+0 * muT(i,1) * (V(i,2)-V(i,1)) / DYP(2)			!Fluid Boundary	
		! North boundary
		TyyT(i,Ny)= (-2.0D+0/3.0D+0)*RHO(i,Ny)*TK(i,Ny) +	&			!Fluid Boundary	
				2.0D+0 * muT(i,Ny) * (V(i,Ny-1)-V(i,Ny-2)) / DYP(Ny-1)		!Fluid Boundary	
		! South boundary
!		TyyT(i,1) = 0.0D+0	!Solid Boundary
!		TxxT(i,1) = 0.0D+0	!Solid Boundary
		! North boundary
!		TyyT(i,Ny) = 0.0D+0	!Solid Boundary
!		TxxT(i,Ny) = 0.0D+0	!Solid Boundary
	end do


	Do i=2,(Nx-2)
		! South boundary
		MUn= (muT(i,1)+muT(i+1,1)) /2.0D+0  !FLUX OJOJO					!Fluid Boundary
		TxyT(i,1)=MUn * ((U(i,2)-U(i,1)) / DYV(1) + (V(i+1,1)-V(i,1)) / DXU(i))		!Fluid Boundary
		! North boundary
		MUn= (muT(i,Ny)+muT(i+1,Ny)) /2.0D+0 !FLUX OJOJO					!Fluid Boundary
		TxyT(i,Ny-1)=MUn * ((U(i,Ny)-U(i,Ny-1)) / DYV(Ny-1) + (V(i+1,Ny-1)-V(i,Ny-1)) / DXU(i)) !Fluid Boundary

		! South boundary
!		TxyT(i,1) = 0.0D+0	!Solid Boundary
		! North boundary
!		TxyT(i,Ny-1) = 0.0D+0	!Solid Boundary

	end do

	Do j=2,(Ny-2)
		! West Boundary
		MUn= flux( muT(1,j), muT(1,j+1), DYV(j), DYP(j+1)/2.0D+0 )			!Fluid Boundary
		TxyT(1,j)=MUn * ((U(1,j+1)-U(1,j)) / DYV(j) + (V(2,j)-V(1,j)) / DXU(1))		!Fluid Boundary 



		! East Boundary		
		MUn= flux( muT(Nx,j), muT(Nx,j+1), DYV(j), DYP(j+1) /2.0D+0 ) 					        !Fluid Boundary
		TxyT(Nx-1,j)=MUn * ((U(Nx-1,j+1)-U(Nx-1,j)) / DYV(j) + (V(Nx,j)-V(Nx-1,j)) / DXU(Nx-1))	!Fluid Boundary

		! West Boundary
!		TxyT(1,j) = 0.0D+0	!Solid Boundary
		! East Boundary
!		TxyT(Nx-1,j) = 0.0D+0	!Solid Boundary

	end do

!if (iter .eq. 2) then
!write(*,*) Mut(:,Ny-1) 
!write(*,*)
!write(*,*) Mut(:,Ny-2) 
!write(*,*)
!write(*,*) Mut(:,Ny-3) 
!write(*,*)
!write(*,*) Txyl(Nx-2,Ny-2) , 'dudy:', (U(Nx-2,Ny-2+1)-U(Nx-2,Ny-2)) / DYV(Ny-2) , 'dvdx:', (V(Nx-2+1,Ny-2)-V(Nx-2,Ny-2)) / DXU(Nx-2)
!write(*,*)
!write(*,*) TxyT(Nx-2,Ny-2) , 'dvdx:', (U(Nx-2,Ny-2+1)-U(Nx-2,Ny-2)) / DYV(Ny-2) , 'dvdx:' ,(V(Nx-2+1,Ny-2)-V(Nx-2,Ny-2)) / DXU(Nx-2)
!stop
!end if

    return
end subroutine stressLT


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					start (U,V,P)						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine start (U,V,P,T,C,Tref,rho,Fep,Fnp,TK,e)	

	implicit none
	integer:: i,j,stat
	real(KIND = DP)  :: Tref
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,C,RHO
	real(KIND = DP), dimension(Nx,Ny),intent(out) :: Fep,Fnp,TK,e
	character(len=8) :: case_jet_string

	Do j=1,Ny
		Do i=1,(Nx-1)

!			U(i,j)=x(i)*y(j)*20.0D+0
			U(i,j)=U_bg
		end do
	end do


	Do j=1,(Ny-1)
		Do i=1,Nx

!			V(i,j)=x(i)*y(j)*50.0D+0
			V(i,j)=0.0D+00
		end do
	end do

!write(*,*) Tref ;stop


	Do j=1,Ny
		Do i=1,Nx

!			T(i,j)=x(i)*y(j)*Tref
			T(i,j)=Th !Th
		end do
	end do

	Do j=1,Ny
		Do i=1,Nx

!			C(i,j)=x(i)*y(j)*Tref
			C(i,j)=C_coflow
		end do
	end do

	Do j=1,Ny
		Do i=1,Nx
!			P(i,j)=x(i)*y(j)*20.0D+0
			P(i,j)= Po
		end do
	end do

	Do j=1,Ny
		Do i=1,Nx
			RHO(i,j)=RHOo
		end do
	end do
	call MASSFLUXES(RHOo,U,V,FEP,FNP,iter)
!	call Mass_Flux(U,V,RHO,Fep,Fnp)

	Do j=1,Ny
		Do i=1,Nx
!			TK(i,j)=x(i)*y(j)*10.0D+0
			TK(i,j)=0.0D+0
		end do
	end do

	Do j=1,Ny
		Do i=1,Nx
!			e(i,j)=x(i)*y(j)*20.0D+0
			e(i,j)=0.0D+0
		end do
	end do

	CALL init_rans_constants()

return
end subroutine start




!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					update					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!


subroutine update_U (Ux,U)	

	implicit none
	integer:: i,j
	real(KIND = DP), dimension(Nx,Ny) :: Ux,U

	Do j=1,Ny
		Do i=1,(Nx-1)
			Ux(i,j)=U(i,j)
		end do
	end do

    return
end subroutine update_U
!					update_V (Ux,U)						!
subroutine update_V (Vx,V)	

	implicit none
	integer:: i,j
	real(KIND = DP), dimension(Nx,Ny) :: Vx,V

	Do j=1,(Ny-1)
		Do i=1,Nx
			Vx(i,j)=V(i,j)
		end do
	end do

    return
end subroutine update_V
!					update_V (Ux,U)						!

subroutine update_P (Ux,U)	

	implicit none
	integer:: i,j
	real(KIND = DP), dimension(Nx,Ny) :: Ux,U

	Do j=1,Ny
		Do i=1,Nx
			Ux(i,j)=U(i,j)
		end do
	end do

    return
end subroutine update_P


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_T 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_T(DXP,DYP,DXU,DYV,X,Y,XU,YV,			&
		CONTER,Cp,RHO,U,V,T,ap,ae,aw,as,an,b,Res_T,Rmax_T,dt,Fep,Fnp,MUt)

	implicit none
		   
	integer:: i,j,npas

	real(KIND = DP)  :: 			&
		Fe,Fw,Fn,Fs,			&
		De,Dw,Dn,Ds,			&
		RES_T,resi,summ,Rmax_T,		&
		Pes,Pen,Pee,Pew,		&
		dt,				&
		GIE,GIW,GIN,GIS,		&
		GIEL,GIWL,GINL,GISL,		&
		GIET,GIWT,GINT,GIST,sc,sp,Sdc		

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X,XU,DXU

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y,YV,DYV

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b,apo,		&
		U,V,T,conter,Cp,rho,		&
		Fep,Fnp,MUt

!-----------------------------------------------------------------------------------------------!
!					varaibles reset 		!
!-----------------------------------------------------------------------------------------------!
																	 
																	 
	do j=1,Ny
		do i=1,Nx							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			ap(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00

			apo(i,j)=0.0D+00
																	 
		end do
	end do 



	do  j=2,Ny-1
		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

      Fe=Fep(i,j)*DYP(j)
	
      Fw=Fep(i-1,j)*DYP(j)
	
      Fn=Fnp(i,j)*DXP(i)
	
      Fs=Fnp(i,j-1)*DXP(i)


!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

	Pr = Cpo * mu(i,j) / conter(i,j)


!	 Calculo da gamma
	GIEL=flux(Mu(i,j),Mu(i+1,j),DXU(i),(DXP(i+1)/2.0D+0) )/Pr

	GIWL=flux(Mu(i,j),Mu(i-1,j),DXU(i-1),(DXP(i-1)/2.0D+0))/Pr
	
	GINL=flux(Mu(i,j),Mu(i,j+1),DYV(j),(DYP(j+1)/2.0D+0))/Pr
	
	GISL=flux(Mu(i,j),Mu(i,j-1),DYV(j-1),(DYP(j-1)/2D+0))/Pr


	   
	GIET=flux(MUT(i,j),MUT(i+1,j),DXU(i),(DXP(i+1)/2.0D+0))	/ Theta_t

	GIWT=flux(MUT(i,j),MUT(i-1,j),DXU(i-1),(DXP(i-1)/2.0D+0)) / Theta_t
	
	GINT=flux(MUT(i,j),MUT(i,j+1),DYV(j),(DYP(j+1)/2.0D+0)) / Theta_t
	
	GIST=flux(MUT(i,j),MUT(i,j-1),DYV(j-1),(DYP(j-1)/2D+0)) / Theta_t


	GIE = GIEL + GIET
	GIW = GIWL + GIWT
	GIN = GINL + GINT
	GIS = GISL + GIST

	if (i==2) 	GIW = Mu(i-1,j)/Pr      + MuT(i-1,j)/Theta_t    ! fluid boundary
	if (i==Nx-1) 	GIE = Mu(i+1,j)/Pr  + MuT(i+1,j)/Theta_t    ! fluid boundary
	if (j==2) 	GIS = Mu(i,j-1)/Pr      + MuT(i,j-1)/Theta_t   	!  fluid boundary
	if (j==Ny-1) 	GIN = Mu(i,j+1)/Pr  + MuT(i,j+1)/Theta_t    !  fluid boundary

	
!if ( (i.eq.3) .and. (j.eq.3) ) then
!	write(*,*) GIE, GIEL 
!	stop
!end if


!	Calculo do flux difusivo
	
	De=GIE*DYP(j)/( (x(i+1)-x(i)) )
	
	Dw=GIW*DYP(j)/( (x(i)-x(i-1)) )
	
	Dn=GIN*DXP(i)/( (y(j+1)-y(j)) )
	
	Ds=GIS*DXP(i)/( (y(j)-y(j-1)) )


!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


	Pee=Fe/De

    Pew=Fw/Dw

	Pen=Fn/Dn

	Pes=Fs/Ds


!-----------------------------------------------------------------------------------------------!
!					High_order_scheme				!
!-----------------------------------------------------------------------------------------------!
	sc = 0.0D+0
	sp = 0.0D+0
	Sdc = 0.0D+0

	if ( (i.gt. 2).and.(i.lt. Nx-1).and.(j.gt. 2).and.(j.lt. Ny-1).and.(n_high_order_scheme.ne. 0) &
		.and. high_order_on_scalars ) then

		call high_order_scheme( i,j,Fe,Fw,Fn,Fs,					&
			T(i-2,j),T(i-1,j),T(i,j),T(i+1,j),T(i+2,j),		&
			T(i,j-2),T(i,j-1),T(i,j+1),T(i,j+2),			&
			dxu(i-2),dxu(i-1),dxu(i),dxu(i+1),			&
			dyv(j-2),dyv(j-1),dyv(j),dyv(j+1),Sdc			)
		

		sc = MAX(Sdc,0.0D+0)
		sp = Min(Sdc / (T(i,j)+1.0D-25),0.0D+0)

	end if

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!

	

	ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

	an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

	as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

	apo(i,j)=(RHO(i,j)*DXP(i)*DYP(j))/dt
	
        ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j) !- sp

        b(i,j)=apo(i,j)*T(i,j) !+ sc + P_rad(i,j)*DXP(i)*DYP(j)/cp(i,j)

	end do
	end do 



!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!
     

!					SOUTH					!
	Do i=1,(Nx)
		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=0.0D+00
		as(i,1)=0.0D+00
		b(i,1) =Th
	end do
   
!					NORTH					!
	Do i=1,(Nx)
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00!0.0D+00
		b(i,Ny) =Tc !T_coflow
	end do
 
!					WEST					!
	Do j=1,Ny
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j) =Th
	end do

!if (street_canyon) then
	if (.TRUE.) then

	!				WEST					!
		Do j=1,Ny
			outside_canyon = Y(j).gt.B1
			if (outside_canyon) then
			ae(1,j)		=	0.0D+00
			aw(1,j)		=	0.0D+00
			ap(1,j)	=	1.0D+00
			an(1,j)		=	0.0D+00
			as(1,j)		=	0.0D+00
			b(1,j)		=       Tc
			end if  	
		end do


	end if

!write(*,*) b(1,:)
!stop

!					EAST					!
	Do j=2,Ny-1
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=0.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j) =Th
	END DO

!if (street_canyon) then
	if (.TRUE.) then

	!				EAST					!
		Do j=2,Ny-1
			outside_canyon = Y(j).gt.B2
			if (outside_canyon) then
			ae(Nx,j)		=	0.0D+00
			aw(Nx,j)		=	0.0D+00
			ap(Nx,j)	=	1.0D+00
			an(Nx,j)		=	0.0D+00
			as(Nx,j)		=	0.0D+00
			b(Nx,j)		=       Tc
			end if  	
		end do


	end if

!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


	summ=0.0D+00
	Rmax_T=0.0D+00

if (street_canyon .eqv. .false.) then 

	do  j=2,Ny-1
		do i=2,Nx-1

		Resi=dabs(ap(i,j)*T(i,j)-ae(i,j)*T(i+1,j)-aw(i,j)*T(i-1,j)-as(i,j)*T(i,j-1)-an(i,j)*T(i,j+1)-b(i,j))
		Rij_T(i,j) = Resi
		summ=summ+(Resi*Resi)

		if (Resi .gt. Rmax_T) Rmax_T=Resi

		end do
	end do

else

	do  j=2,Ny-1
	do i=2,Nx-1
		outside_canyon =(X(i).ge.((BD+BW)-1.0D+00)).and.(Y(j).le.B2+1.0D+00).or.&
		(X(i).le.BW+1.0D+00).and.(Y(j).le.B1+1.0D+00)
			if (outside_canyon) then
			else
				Resi=dabs(ap(i,j)*T(i,j)-ae(i,j)*T(i+1,j)-aw(i,j)*T(i-1,j)-as(i,j)*T(i,j-1)-an(i,j)*T(i,j+1)-b(i,j))
				Rij_T(i,j) = Resi
				summ=summ+(Resi*Resi)

				if (Resi .gt. Rmax_T) Rmax_T=Resi
			end if
	end do
	end do
end if

	Res_T=dsqrt(summ)

	
	return

end Subroutine coef_T

!**************************************************************************
!*******      FUNCIÓN FLUX (RESUELVE EL FLUJO EN LA INTERFACE)      *******        
!**************************************************************************
!C        
      FUNCTION  FLUX2(F1,F2,D12,D)
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!C
!CC    CALCULO DEL FLUJO EN LA INTERFACE
!C
      FLUX2=F1+(((F2-F1)/D12)*(D/2.0D+00))
!C
      RETURN  
      END
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_C 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_C(DXP,DYP,DXU,DYV,X,Y,XU,YV,			&
		CONTER,Cp,RHO,U,V,C,ap,ae,aw,as,an,b,Res_C,Rmax_C,dt,Fep,Fnp,MUt)

	implicit none
		   
	integer:: i,j,npas

	real(KIND = DP)  :: 			&
		Fe,Fw,Fn,Fs,			&
		De,Dw,Dn,Ds,			&
		RES_C,resi,summ,Rmax_C,		&
		Pes,Pen,Pee,Pew,		&
		dt,				&
		GIE,GIW,GIN,GIS,		&
		GIEL,GIWL,GINL,GISL,		&
		GIET,GIWT,GINT,GIST,sc,sp,Sdc		

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X,XU,DXU

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y,YV,DYV

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b,apo,		&
		U,V,C,conter,Cp,rho,		&
		Fep,Fnp,MUt

!-----------------------------------------------------------------------------------------------!
!					varaibles reset 		!
!-----------------------------------------------------------------------------------------------!
																	 
																	 
	do j=1,Ny
		do i=1,Nx							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			ap(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00

			apo(i,j)=0.0D+00
																	 
		end do
	end do 



	do  j=2,Ny-1
		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

      Fe=Fep(i,j)*DYP(j)
	
      Fw=Fep(i-1,j)*DYP(j)
	
      Fn=Fnp(i,j)*DXP(i)
	
      Fs=Fnp(i,j-1)*DXP(i)


!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

	Pr = Sch


!	 Calculo da gamma
	GIEL=flux(Mu(i,j),Mu(i+1,j),DXU(i),(DXP(i+1)/2.0D+0) )/Pr

	GIWL=flux(Mu(i,j),Mu(i-1,j),DXU(i-1),(DXP(i-1)/2.0D+0))/Pr
	
	GINL=flux(Mu(i,j),Mu(i,j+1),DYV(j),(DYP(j+1)/2.0D+0))/Pr
	
	GISL=flux(Mu(i,j),Mu(i,j-1),DYV(j-1),(DYP(j-1)/2D+0))/Pr


	   
	GIET=flux(MUT(i,j),MUT(i+1,j),DXU(i),(DXP(i+1)/2.0D+0))	/ Theta_t

	GIWT=flux(MUT(i,j),MUT(i-1,j),DXU(i-1),(DXP(i-1)/2.0D+0)) / Theta_t
	
	GINT=flux(MUT(i,j),MUT(i,j+1),DYV(j),(DYP(j+1)/2.0D+0)) / Theta_t
	
	GIST=flux(MUT(i,j),MUT(i,j-1),DYV(j-1),(DYP(j-1)/2D+0)) / Theta_t


	GIE = GIEL + GIET
	GIW = GIWL + GIWT
	GIN = GINL + GINT
	GIS = GISL + GIST

	if (i==2) 	GIW = Mu(i-1,j)/Pr      + MuT(i-1,j)/Theta_t    ! fluid boundary
	if (i==Nx-1) 	GIE = Mu(i+1,j)/Pr  + MuT(i+1,j)/Theta_t    ! fluid boundary
	if (j==2) 	GIS = Mu(i,j-1)/Pr      + MuT(i,j-1)/Theta_t   	!  fluid boundary
	if (j==Ny-1) 	GIN = Mu(i,j+1)/Pr  + MuT(i,j+1)/Theta_t    !  fluid boundary

	
!if ( (i.eq.3) .and. (j.eq.3) ) then
!	write(*,*) GIE, GIEL 
!	stop
!end if


!	Calculo do flux difusivo
	
	De=GIE*DYP(j)/( (x(i+1)-x(i)) )
	
	Dw=GIW*DYP(j)/( (x(i)-x(i-1)) )
	
	Dn=GIN*DXP(i)/( (y(j+1)-y(j)) )
	
	Ds=GIS*DXP(i)/( (y(j)-y(j-1)) )


!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


	Pee=Fe/De

    Pew=Fw/Dw

	Pen=Fn/Dn

	Pes=Fs/Ds


!-----------------------------------------------------------------------------------------------!
!					High_order_scheme				!
!-----------------------------------------------------------------------------------------------!
	sc = 0.0D+0
	sp = 0.0D+0
	Sdc = 0.0D+0

	if ( (i.gt. 2).and.(i.lt. Nx-1).and.(j.gt. 2).and.(j.lt. Ny-1).and.(n_high_order_scheme.ne. 0) &
		.and. high_order_on_scalars ) then

		call high_order_scheme( i,j,Fe,Fw,Fn,Fs,					&
			C(i-2,j),C(i-1,j),C(i,j),C(i+1,j),C(i+2,j),		&
			C(i,j-2),C(i,j-1),C(i,j+1),C(i,j+2),			&
			dxu(i-2),dxu(i-1),dxu(i),dxu(i+1),			&
			dyv(j-2),dyv(j-1),dyv(j),dyv(j+1),Sdc			)
		

		sc = MAX(Sdc,0.0D+0)
		sp = Min(Sdc / (C(i,j)+1.0D-25),0.0D+0)

	end if

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!

	

	ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

	an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

	as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

	apo(i,j)=(RHO(i,j)*DXP(i)*DYP(j))/dt
	
        ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j) !- sp

        b(i,j)=apo(i,j)*C(i,j) !+ sc

!How to parameterize
	if ( (i.gt.Nx/2).and.(i.le.(Nx/2)+1).and.(j==3)) then
		b(i,j)=b(i,j)+1.0D-09
	end if

	end do
	end do 



!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!
!Is it right?     

!					SOUTH					!
	Do i=1,(Nx)
		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=1.0D+00
		as(i,1)=0.0D+00
		b(i,1) =0.0D+00
	end do
   
!					NORTH					!
	Do i=1,(Nx)
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00!0.0D+00
		b(i,Ny) =C_coflow
	end do
 
!					WEST					!
	Do j=1,Ny
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j) =C_coflow
	end do

!write(*,*) b(1,:)
!stop

!					EAST					!
	Do j=2,Ny-1
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=1.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j) =0.0D+00
	END DO



!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


	summ=0.0D+00
	Rmax_C=0.0D+00

if (street_canyon .eqv. .false.) then 

	do  j=2,Ny-1
		do i=2,Nx-1

		Resi=dabs(ap(i,j)*C(i,j)-ae(i,j)*C(i+1,j)-aw(i,j)*C(i-1,j)-as(i,j)*C(i,j-1)-an(i,j)*C(i,j+1)-b(i,j))
		Rij_C(i,j) = Resi
		summ=summ+(Resi*Resi)

		if (Resi .gt. Rmax_C) Rmax_C=Resi

		end do
	end do

else

	do  j=2,Ny-1
	do i=2,Nx-1
		outside_canyon =(X(i).ge.((BD+BW)-1.0D+00)).and.(Y(j).le.B2+1.0D+00).or.&
		(X(i).le.BW+1.0D+00).and.(Y(j).le.B1+1.0D+00)
			if (outside_canyon) then
			else
				Resi=dabs(ap(i,j)*C(i,j)-ae(i,j)*C(i+1,j)-aw(i,j)*C(i-1,j)-as(i,j)*C(i,j-1)-an(i,j)*C(i,j+1)-b(i,j))
				Rij_C(i,j) = Resi
				summ=summ+(Resi*Resi)

				if (Resi .gt. Rmax_C) Rmax_C=Resi
			end if
	end do
	end do

end if

	Res_C=dsqrt(summ)
	
	return

end Subroutine coef_C

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_TK 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_TK(DXP,DYP,DXU,DYV,X,Y,XU,YV,								&
		Mu,Mut,RHO,U,V,TKx,e,ap,ae,aw,as,an,b,Res_Tk,Rmax_Tk,dt,Fep,Fnp,T,Txxt,Tyyt,Txyt,Pk,Gk)

	implicit none
		   
	integer:: i,j,npas

	real(KIND = DP)  :: 			&
		Fe,Fw,Fn,Fs,			&
		De,Dw,Dn,Ds,			&
		RES_Tk,resi,summ,Rmax_Tk,	&
		Pes,Pen,Pee,Pew,		&
		dt,dtdy,			&
		GIE,GIW,GIN,GIS,Txyt_P,		&
		Un, Up, Us, Vw, Vp, Ve, 	&
		dudx, dudy, dvdx, dvdy,sp,sc,	&
		GIEL,GIWL,GINL,GISL,		&
		GIET,GIWT,GINT,GIST,sp_hos,sc_hos,Sdc

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X,DXU,XU

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y,DYV,YV

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b,apo,		&
		U,V,Tkx,mu,MuT,e,rho,		&
		Fep,Fnp,Pk,Gk,T,Txxt,Tyyt,Txyt


!-----------------------------------------------------------------------------------------------!
!					varaibles reset 		!
!-----------------------------------------------------------------------------------------------!
																	 
																	 
	do j=1,Ny
		do i=1,Nx							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			ap(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00

			apo(i,j)=0.0D+00
																	 
		end do
	end do 



	do  j=2,Ny-1
		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

      Fe=Fep(i,j)*DYP(j)
	
      Fw=Fep(i-1,j)*DYP(j)
	
      Fn=Fnp(i,j)*DXP(i)
	
      Fs=Fnp(i,j-1)*DXP(i)




!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!


!	 Calculo da gamma
      
	GIEL = flux(Mu(i,j),	Mu(i+1,j),	DXU(i),(DXP(i+1)/2.0D+0) )
	GIET = flux(MuT(i,j),	MuT(i+1,j),	DXU(i),(DXP(i+1)/2.0D+0))  	/Theta_k

	GIWL = flux(Mu(i,j)	,Mu(i-1,j)	,DXU(i-1),(DXP(i-1)/2.0D+0))
	GIWT = flux(MuT(i,j)	,MuT(i-1,j)	,DXU(i-1),(DXP(i-1)/2.0D+0)) 	/Theta_k
	
	GINL = flux(Mu(i,j),	Mu(i,j+1)	,DYV(j),(DYP(j+1)/2.0D+0))
	GINT = flux(MuT(i,j),	MuT(i,j+1)	,DYV(j),(DYP(j+1)/2.0D+0)) 	/Theta_k
	
	GISL = flux(Mu(i,j)	,Mu(i,j-1)	,DYV(j-1),(DYP(j-1)/2D+0))
	GIST = flux(MuT(i,j)	,MuT(i,j-1)	,DYV(j-1),(DYP(j-1)/2.0D+0)) 	/Theta_k


	GIE = GIEL + GIET
	GIW = GIWL + GIWT
	GIN = GINL + GINT
	GIS = GISL + GIST

	if (i==2) 	GIW = Mu(i-1,j)  	 + 	MuT(i-1,j)/Theta_k	! fluid boundary
	if (i==Nx-1) 	GIE = Mu(i+1,j)  +  MuT(i+1,j)/Theta_k	! fluid boundary
	if (j==2) 	GIS = Mu(i,j-1)      +  MuT(i,j-1)/Theta_k   			! Add turbulent viscosity if fluid boundary
	if (j==Ny-1) 	GIN = Mu(i,j+1)  +  MuT(i,j+1)/Theta_k      			! Add turbulent viscosity if fluid boundary


!	Calculo do flux difusivo
	
	De=GIE*DYP(j)/( (x(i+1)-x(i)) )
	
	Dw=GIW*DYP(j)/( (x(i)-x(i-1)) )
	
	Dn=GIN*DXP(i)/( (y(j+1)-y(j)) )
	
	Ds=GIS*DXP(i)/( (y(j)-y(j-1)) )


!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


	Pee=Fe/De

    Pew=Fw/Dw

	Pen=Fn/Dn

	Pes=Fs/Ds


!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------!
!					High_order_scheme				!
!-----------------------------------------------------------------------------------------------!
	sc_hos = 0.0D+0
	sp_hos = 0.0D+0
	Sdc=0.0D+0

	if ( (i.gt. 2).and.(i.lt. Nx-1).and.(j.gt. 2).and.(j.lt. Ny-1).and.(n_high_order_scheme.ne. 0) &
		.and. high_order_on_scalars ) then

		call high_order_scheme( i,j,Fe,Fw,Fn,Fs,						&
			TKx(i-2,j),TKx(i-1,j),TKx(i,j),TKx(i+1,j),TKx(i+2,j),		&
			TKx(i,j-2),TKx(i,j-1),TKx(i,j+1),TKx(i,j+2),				&
			dxu(i-2),dxu(i-1),dxu(i),dxu(i+1),							&
			dyv(j-2),dyv(j-1),dyv(j),dyv(j+1),Sdc			)
		

		sc_hos = MAX(Sdc,0.0D+0)
		sp_hos = Min(Sdc / (Tkx(i,j)+1.0D-25),0.0D+0)

	end if
!-----------------------------------------------!
!	Source terms				!

	Call rans_source_terms_TK(i,j,u,v,Mut,mu,T,rho,e,TKx,sc,sp)

!-----------------------------------------------!
!	Coeficients				!	

	ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

	an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

	as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

	apo(i,j)=(RHO(i,j)*DXP(i)*DYP(j))/dt
	
     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j) - Sp * DXP(i)*DYP(j) - sp_hos

      	b(i,j)=apo(i,j)*TKx(i,j) + Sc  * DXP(i)*DYP(j) + sc_hos

		end do
	end do 



!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!
     

!					SOUTH					!
	Do i=1,(Nx)
		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=0.0D+00
		as(i,1)=0.0D+00
		b(i,1)=0.0D+00
	end do
   
!					NORTH					!
	Do i=1,(Nx)
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=1.0D+00
		b(i,Ny)=0.0D+00
	end do
	
 
!					WEST					!
	Do j=2,Ny-1
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j) = 0.0D+00
	end do


!					EAST					!
	Do j=2,Ny-1
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=1.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j) =0.0D+00
	END DO



!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


	summ=0.0D+00
	Rmax_Tk=0.0D+00

	do  j=3,Ny-2
		do i=3,Nx-2

		Resi=dabs(ap(i,j)*Tk(i,j)-ae(i,j)*Tk(i+1,j)-aw(i,j)*Tk(i-1,j)-as(i,j)*Tk(i,j-1)-an(i,j)*Tk(i,j+1)-b(i,j))
		summ=summ+(Resi*Resi)

		if (Resi .gt. Rmax_Tk) Rmax_Tk=Resi

		end do
	end do

	Res_Tk=dsqrt(summ)
	
	return

end Subroutine coef_TK


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_e 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine COEF_e(DXP,DYP,DXU,DYV,X,Y,XU,YV,			&
		Mu,MuT,RHO,U,V,Tk,e,ap,ae,aw,as,an,b,Res_e,Rmax_e,dt,Fep,Fnp,T,Txxt,Tyyt,Txyt,Pk,Gk)

	implicit none
		   
	integer:: i,j,npas

	real(KIND = DP)  :: 			&
		Fe,Fw,Fn,Fs,			&
		De,Dw,Dn,Ds,			&
		RES_e,resi,summ,Rmax_e,		&
		Pes,Pen,Pee,Pew,		&
		dt,bb,				&
		GIE,GIW,GIN,GIS,		&
		C3e, Dis,Txyt_P,		&
		Up,Vp,sp,sc,			&
		GIEL,GIWL,GINL,GISL,		&
		GIET,GIWT,GINT,GIST,sp_hos,sc_hos,Sdc

	real(KIND = DP), dimension(Nx) :: 	&
		DXP,X,XU,DXU

	real(KIND = DP), dimension(Ny) :: 	&
		DYP,Y,YV,DYV

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ae,aw,an,as,ap,b,apo,		&
		U,V,Tk,e,mu,mut,rho,		&
		Fep,Fnp,T,Txxt,Tyyt,Txyt,Pk,Gk


!-----------------------------------------------------------------------------------------------!
!					varaibles reset 		!
!-----------------------------------------------------------------------------------------------!
																	 
	do j=1,Ny
		do i=1,Nx							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			ap(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00

			apo(i,j)=0.0D+00
																	 
		end do
	end do 



	do  j=2,Ny-1
		do  i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

      Fe=Fep(i,j)*DYP(j)
	
      Fw=Fep(i-1,j)*DYP(j)
	
      Fn=Fnp(i,j)*DXP(i)
	
      Fs=Fnp(i,j-1)*DXP(i)




!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!
!	 Calculo da gamma
      
	GIEL = flux(Mu(i,j) ,	Mu(i+1,j)	,DXU(i),(DXP(i+1)/2.0D+0) )
	GIET = flux(MuT(i,j),	MuT(i+1,j)	,DXU(i),(DXP(i+1)/2.0D+0)) 		/Theta_e

	GIWL = flux(Mu(i,j),	Mu(i-1,j)	,DXU(i-1),(DXP(i-1)/2.0D+0))
	GIWT = flux(MuT(i,j),	MuT(i-1,j)	,DXU(i-1),(DXP(i-1)/2.0D+0)) 		/Theta_e
	
	GINL = flux(Mu(i,j)	,Mu(i,j+1),	DYV(j),(DYP(j+1)/2.0D+0))
	GINT = flux(MuT(i,j)	,MuT(i,j+1),	DYV(j),(DYP(j+1)/2.0D+0)) 		/Theta_e
	
	GISL = flux(Mu(i,j),    Mu(i,j-1)	,DYV(j-1),(DYP(j-1)/2D+0))
	GIST = flux(MuT(i,j),	MuT(i,j-1)	,DYV(j-1),(DYP(j-1)/2.0D+0)) 		/Theta_e


	GIE = GIEL + GIET
	GIW = GIWL + GIWT
	GIN = GINL + GINT
	GIS = GISL + GIST

	if (i==2) 		GIW = Mu(i-1,j)  + MuT(i-1,j)/Theta_e	! fluid boundary
	if (i==Nx-1) 	GIE = Mu(i+1,j)  + MuT(i+1,j)/Theta_e	! fluid boundary
	if (j==2) 		GIS = Mu(i,j-1)  + MuT(i,j-1)/Theta_e   ! fluid boundary
	if (j==Ny-1) 	GIN = Mu(i,j+1)  + MuT(i,j+1)/Theta_e   !  fluid boundary


!	Calculo do flux difusivo
	
      De=GIE*DYP(j)/( (x(i+1)-x(i)) )
	
      Dw=GIW*DYP(j)/( (x(i)-x(i-1)) )
	
      Dn=GIN*DXP(i)/( (y(j+1)-y(j)) )
	
      Ds=GIS*DXP(i)/( (y(j)-y(j-1)) )


!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


	Pee=Fe/De

    Pew=Fw/Dw

	Pen=Fn/Dn

	Pes=Fs/Ds


!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					High_order_scheme				!
!-----------------------------------------------------------------------------------------------!
	sc_hos = 0.0D+0
	sp_hos = 0.0D+0
	Sdc=0.0D+0

	if ( (i.gt. 2).and.(i.lt. Nx-1).and.(j.gt. 2).and.(j.lt. Ny-1).and.(n_high_order_scheme.ne. 0) &
		.and. high_order_on_scalars ) then

		call high_order_scheme( i,j,Fe,Fw,Fn,Fs,						&
			e(i-2,j),e(i-1,j),e(i,j),e(i+1,j),e(i+2,j),					&
			e(i,j-2),e(i,j-1),e(i,j+1),e(i,j+2),						&
			dxu(i-2),dxu(i-1),dxu(i),dxu(i+1),							&
			dyv(j-2),dyv(j-1),dyv(j),dyv(j+1),Sdc						)
		

		sc_hos = MAX(Sdc,0.0D+0)
		sp_hos = Min(Sdc / (e(i,j)+1.0D-25),0.0D+0)

	end if
!-----------------------------------------------!
!	Source terms				!

	Call rans_source_terms_e(i,j,u,v,Mut,mu,T,rho,e,TK,sc,sp)




	ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

	aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

	an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

	as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

	apo(i,j)=(RHO(i,j)*DXP(i)*DYP(j))/dt
	
     	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo(i,j)  - Sp * DXP(i)*DYP(j) - sp_hos

      	b(i,j)=apo(i,j)*e(i,j) + Sc  * DXP(i)*DYP(j) + sc_hos

		end do
	end do 



!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS															!
!-----------------------------------------------------------------------------------------------!

!					SOUTH					!
	Do i=1,(Nx)
		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=0.0D+00
		as(i,1)=0.0D+00
		b(i,1) =0.0D+00
	end do
   
!					NORTH					!
	Do i=1,(Nx)
		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00
		b(i,Ny) =0.0D+00! rans_boundary_e(i,Ny)
	end do
	
 
!					WEST					!
	Do j=2,Ny-1
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j) = 0.0D+00
	end do


!					EAST					!
	Do j=2,Ny-1
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=1.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j) =0.0D+00
	END DO



!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


	summ=0.0D+00
	Rmax_e=0.0D+00

	do  j=3,Ny-2
		do i=3,Nx-2

		Resi=dabs(ap(i,j)*e(i,j)-ae(i,j)*e(i+1,j)-aw(i,j)*e(i-1,j)-as(i,j)*e(i,j-1)-an(i,j)*e(i,j+1)-b(i,j))
		summ=summ+(Resi*Resi)

		if (Resi .gt. Rmax_e) Rmax_e=Resi

		end do
	end do

	Res_e=dsqrt(summ)

	
	return

end Subroutine coef_e


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_U 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine coef_U(DXP,DYP,DXU,DYV,X,Y,XU,YV,				&
			Ux,Vx,Px,					&
			du,ae,aw,as,an,b,RES_U,Rmax_U,mu,apu,dt,RHO,	&
			TxxL, TyyL, TxyL, Fep, Fnp, MUt, TxxT, TyyT, TxyT,APUNB,U_bg)	     

	implicit none

	integer:: i,j,npas

	real(KIND = DP)  :: Fe,Fw,Fn,Fs,			&
			    De,Dw,Dn,Ds,			&
			    Uo,RES_U,resi,summ,			&
		            Pes,Pen,Pee,Pew,			&
			    Rmax_U,apo,dt,RHOn,RHOs,		&
			    dudx_e ,dudx_w , 			&
			    dudy_n,dudy_s ,			&
			    Source1, Source2,			&
			    MUL_n ,MUL_s,MUT_n ,MUT_s,		&
			    MUL_p,MUT_P,			&
			    MUL_J1,MUL_J_1,MUT_J1,MUT_J_1, 	&
			    GIE,GIW,GIN,GIS,			&
				D1P,D2N,D3S,F1	,MU_n,MU_s,sp_hos,sc_hos,Sdc, &
				U_bg					

	real(KIND = DP), dimension(Nx) :: 		&
		DXP,DXU,X,XU

	real(KIND = DP), dimension(Ny) :: 		&
		DYP,Y,DYV,YV

	real(KIND = DP), dimension(Nx,Ny) :: 		&
		ae,aw,an,as,apu,b,			&
		Ux,Vx,Px,mu,du,				&
		RHO,APUNB,M_SDC
	real(KIND = DP), dimension(Nx,Ny), intent(in) :: TxxL, TyyL, TxyL, Fep, Fnp, MUt, TxxT, TyyT, TxyT
!-----------------------------------------------------------------------------------------------!
!					varaibles reset 		!
!-----------------------------------------------------------------------------------------------!
																	 
																	 
	do j=1,Ny
		do i=1,Nx-1							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			apu(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00

			du(i,j)=0.0D+00	
																	 
		end do
	end do 


	do  j=2,Ny-1
		do  i=2,Nx-2

!-----------------------------------------------------------------------------------------------!
!				Mass flux crossing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

			Fe=DYP(j)*( Fep(i,j) + Fep(i+1,j) )/2.0D+00
	
			Fw=DYP(j)*(  Fep(i,j) + Fep(i-1,j) )/2.0D+00
	
			Fn=DXU(i)*flux( Fnp(i,j), Fnp(i+1,j), DXU(i), X(i+1)-XU(i) )
	
			Fs=DXU(i)*flux( Fnp(i,j-1), Fnp(i+1,j-1), DXU(i), X(i+1)-XU(i) )

!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

	!-----------------------------------------------------------------------------------------------!
	!				Laminar viscosity (MuL)						!
	!-----------------------------------------------------------------------------------------------!
		MuL_j1 = flux(Mu(i,j+1),Mu(i+1,j+1),DXU(i),DXP(i+1)/2.0D+0)
		MuL_p = flux(Mu(i,j),Mu(i+1,j),DXU(i),DXP(i+1)/2.0D+0)

		MuL_n = flux(MuL_p,MuL_j1,DYV(j),DYP(j+1)/2.0D+0)

		MuL_j_1 = flux(Mu(i,j-1),Mu(i+1,j-1),DXU(i),DXP(i+1)/2.0D+0)
		MuL_s = flux(MuL_j_1,MuL_p,DYV(j-1),DYP(j)/2.0D+0)

	!-----------------------------------------------------------------------------------------------!
	!				Turbulent viscosity (MuT)					!
	!-----------------------------------------------------------------------------------------------!
		MuT_j1 = flux(MuT(i,j+1),MuT(i+1,j+1),DXU(i),DXP(i+1)/2.0D+0)
		MuT_p = flux(MuT(i,j),MuT(i+1,j),DXU(i),DXP(i+1)/2.0D+0)

		MuT_n = flux(MuT_p,MuT_j1,DYV(j),DYP(j+1)/2.0D+0)

		MuT_j_1 = flux(MuT(i,j-1),MuT(i+1,j-1),DXU(i),DXP(i+1)/2.0D+0)
		MuT_s = flux(MuT_j_1,MuT_p,DYV(j-1),DYP(j)/2.0D+0)

	!-----------------------------------------------------------------------------------------------!
	!				Total viscosity 						!
	!-----------------------------------------------------------------------------------------------!
		GIE = (MU(i+1,j)+MUt(i+1,j))
		GIW = (MU(i,j)+MUt(i,j))
		GIN = MuL_n + MuT_n
		GIS = MuL_s + MuT_s

		if (j.eq.2) 		GIS= MuL_j_1 + MuT_j_1
		if (j.eq.(Ny-1)) 	GIN= MuL_j1  + MuT_j1 
			

			De=( GIE * DYP(j) ) / ( XU(i+1) - XU(i) )
      
			Dw=( GIW * DYP(j) ) / ( XU(i) - XU(i-1) )

			Dn=( GIN * DXU(i) ) / (Y(j+1)-Y(j)) 
    
      		Ds=( GIS * DXU(i) ) / (Y(j) - Y(j-1))



!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


			Pee=Fe/De

      		Pew=Fw/Dw

			Pen=Fn/Dn

			Pes=Fs/Ds

!-----------------------------------------------------------------------------------------------!
!					High_order_scheme				!
!-----------------------------------------------------------------------------------------------!
	sc_hos = 0.0D+0
	sp_hos = 0.0D+0
	Sdc =0.0D+0
	!if ( (i.gt. 2).and.(i.lt. Nx-2).and.(j.gt. 2).and.(j.lt. Ny-1).and.(n_high_order_scheme.ne. 0) &
	!	.and. high_order_on_velocities ) then
	if ( (n_high_order_scheme.ne.0) .and. high_order_on_velocities_u ) then	

		call high_order_scheme_u( i,j,Fe,Fw,Fn,Fs,						&
			Ux(i-2,j),Ux(i-1,j),Ux(i,j),Ux(i+1,j),Ux(i+2,j),			&
			Ux(i,j-2),Ux(i,j-1),Ux(i,j+1),Ux(i,j+2),					&
			dxp(i-1),dxp(i),dxp(i+1),dxp(i+2),							&
			dyv(j-2),dyv(j-1),dyv(j),dyv(j+1),Sdc						)

	if ( (Sdc / (Ux(i,j)+1.0D-30)) .gt. 0.0D+0) then 

			sc_hos = Sdc
			sp_hos=0.0D+0

		else

			sc_hos=0.0D+0
			sp_hos = Sdc / (Ux(i,j)+1.0D-30)

		end if

	end if

!Sdc_ij_u(i,j)=
!-----------------------------------------------------------------------------------------------!
!					Source terms from turbulent stresses			!
!-----------------------------------------------------------------------------------------------!
			dudx_e = ( Ux(i+1,j) - Ux(i,j)  ) / DXP(i+1)
			dudx_w = ( Ux(i,j)   - Ux(i-1,j)) / DXP(i)	
	
			Source1 = ( TxxL(i+1,j)-TxxL(i,j) + TxxT(i+1,j)-TxxT(i,j) &
				 - GIE*dudx_e +  GIW*dudx_w )*DYP(j)

			dudy_n = ( Ux(i,j+1) - Ux(i,j)  ) / DYV(j)
			dudy_s = ( Ux(i,j)   - Ux(i,j-1)) / DYV(j-1)

			Source2 = ( TxyL(i,j)-TxyL(i,j-1) + TxyT(i,j)-TxyT(i,j-1) &
				- GIN*dudy_n +  GIS*dudy_s )*DXU(i)	

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!
	

			ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

			aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

			an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

			as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

			apo=(flux(RHO(i,j),RHO(i+1,j),DXU(i),DXP(i+1)/2.0D+0)*DXU(i)*DYP(j)) /dt
	
     			apu(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo! - sp_hos

      			b(i,j)=( Px(i,j)-Px(i+1,j) )*DYP(j)+apo*Ux(i,j)!	+  &
				!Source1 + Source2 + sc_hos



		end do
	end do 

!if (iter .eq. 3) then 
!		write(*,*) M_SDC					
!end if 	


!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!
     
!Is it right?
!					SOUTH					!
	Do i=1,(Nx-1)
		ae(i,1)		=	0.0D+00
		aw(i,1)		=	0.0D+00
		apu(i,1)	=	1.0D+00
		an(i,1)		=	0.0D+00
		as(i,1)		=	0.0D+00
		b(i,1)		= 	0.0D+00
	end do
   
!					NORTH					!
	Do i=1,(Nx-1)
		ae(i,Ny)	=	0.0D+00
		aw(i,Ny)	=	0.0D+00
		apu(i,Ny)	=	1.0D+00
		an(i,Ny)	=	0.0D+00
		as(i,Ny)	=	1.0D+00!0.0D+00
		b(i,Ny)		=	0.0D+00!u_coflow
	end do
	

!					WEST					!
	Do j=2,Ny-1
		ae(1,j)		=	0.0D+00
		aw(1,j)		=	0.0D+00
		apu(1,j)	=	1.0D+00
		an(1,j)		=	0.0D+00
		as(1,j)		=	0.0D+00
		b(1,j)		= 	0.0D+00
	end do

	!if (street_canyon) then
	if (.TRUE.) then

	!				WEST					!
		Do j=2,Ny-1
			outside_canyon = Y(j).gt.B1
			if (outside_canyon) then
			ae(1,j)		=	0.0D+00
			aw(1,j)		=	0.0D+00
			apu(1,j)	=	1.0D+00
			an(1,j)		=	0.0D+00
			as(1,j)		=	0.0D+00
			b(1,j)		=       U_bg
			end if  	
		end do


	end if


!					EAST					!
	Do j=2,Ny-1
		ae(Nx-1,j)	=	0.0D+00
		aw(Nx-1,j)	=	1.0D+00
		apu(Nx-1,j)	=	1.0D+00
		an(Nx-1,j)	=	0.0D+00
		as(Nx-1,j)	=	0.0D+00
		b (Nx-1,j)	=	0.0D+00
	END DO




!-----------------------------------------------------------------------------------------------!
!!			d's Matrix determination for SIMPLE algorithm										!
!-----------------------------------------------------------------------------------------------!

	do j=1,Ny
		do i=1,Nx-1

			IF (ALGORITMO.EQ.1) APUNB(I,J)=0.0D+00              	        !ALGORITMO SIMPLE
			IF (ALGORITMO.EQ.2) APUNB(I,J)=AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)  !ALGORITMO SIMPLEC

			!if (apu(i,j) == APUNB(I,J)) then
			!	du = 0.0D+0
			!else
				du(i,j)=DYP(j)/(apu(i,j)-APUNB(I,J))
			!end if

 		end do
	end do

!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


	summ=0.0D+00
	Rmax_U=0.0D+00

if (street_canyon .eqv. .false.) then 

	do  j=2,Ny-1
		do i=2,Nx-2

		Resi=dabs(apu(i,j)*Ux(i,j)-ae(i,j)*Ux(i+1,j)-aw(i,j)*Ux(i-1,j)-as(i,j)*Ux(i,j-1)-an(i,j)*Ux(i,j+1)-b(i,j))
		summ=summ+(Resi*Resi)
		Rij_u(i,j) = Resi
		if (Resi .gt. Rmax_U) Rmax_U=Resi

		end do
	end do

else

	do j=2,Ny-1
	do i=2,Nx-2
		outside_canyon =(X(i).ge.((BD+BW)-1.0D+00)).and.(Y(j).le.B2+1.0D+00).or.&
		(X(i).le.BW+1.0D+00).and.(Y(j).le.B1+1.0D+00)
			if (outside_canyon) then
			else
				Resi=dabs(apu(i,j)*Ux(i,j)-ae(i,j)*Ux(i+1,j)-aw(i,j)*Ux(i-1,j)-as(i,j)*Ux(i,j-1)-an(i,j)*Ux(i,j+1)-b(i,j))
				summ=summ+(Resi*Resi)
				Rij_u(i,j) = Resi
				if (Resi .gt. Rmax_U) Rmax_U=Resi
			end if
	end do
	end do

end if

	Res_U=dsqrt(summ)
	
	return
end Subroutine coef_U




!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_V 						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!


Subroutine COEF_V(DXP,DYP,DXU,DYV,X,Y,XU,YV,				&
		  Ux,Vx,Px,dv,ae,aw,as,an,b,				&
		  RES_V,Rmax_V,mu,apv,dt,Tref,RHO,RHOmed,		&
		  TxxL, TyyL, TxyL, Fep, Fnp, MuT, TxxT, TyyT, TxyT,APVNB)	

	implicit none

	integer:: i,j,npas

	real(KIND = DP)  :: Fe,Fw,Fn,Fs,De,Dw,Dn,Ds,RES_V,resi,summ,Pes,Pen,Pee,Pew,Rmax_V,  &
			    apo,dt,Tref,RHOmed,RHOw,RHOe,	&
			    dvdx_e ,dvdx_w , 			&
			    dvdy_n,dvdy_s ,			&
			    Source1, Source2,			&
			    MUL_e ,MUL_w,MUt_e ,MUt_w,		&
			    MuT_ne,MuL_ne,MuT_se,MuL_se,	&
			    MuT_nw,MuL_nw,MuT_sw,MuL_sw,	&
			    GIE,GIW,GIN,GIS,			&
			    D1P,D2E,D3W,F1,vv_in,sc_hos,sp_hos,Sdc,Pw,Pm,outflow			  

	real(KIND = DP), dimension(Nx) :: 		&
		DXP,DXU,X,XU

	real(KIND = DP), dimension(Ny) :: 		&
		DYP,Y,DYV,YV

	real(KIND = DP), dimension(Nx,Ny) :: 		&
		ae,aw,an,as,apv,b,			&
		Ux,Vx,Px,mu,dv,			&
		RHO,APVNB
	real(KIND = DP), dimension(Nx,Ny), intent(in) :: TxxL, TyyL, TxyL, TxxT, TyyT, TxyT, Fep, Fnp, MuT 


!-----------------------------------------------------------------------------------------------!
!					varaibles reset 					!
!-----------------------------------------------------------------------------------------------!

																	 
																	 
	do j=1,Ny
		do i=1,Nx							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			apv(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00
			
			dv(i,j)=0.0D+00
																	 
		end do
	end do 




	do  j=2,Ny-2
		do  i=2,Nx-1
!-----------------------------------------------------------------------------------------------!
!				Mass flux corssing U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

	
			Fn=DXP(i)*( Fnp(i,j) + Fnp(i,j+1) )/2.0D+00
	
			Fs=DXP(i)*( Fnp(i,j) + Fnp(i,j-1) )/2.0D+00

			Fe=DYV(j) * flux ( Fep(i,j+1),Fep(i,j),DYV(j), YV(j)-Y(j) )
	
			Fw=DYV(j) * flux ( Fep(i-1,j+1),Fep(i-1,j),DYV(j), YV(j)-Y(j) )

!-----------------------------------------------------------------------------------------------!
!				DIFFUSION IN THE U-VC FACES					!
!-----------------------------------------------------------------------------------------------!

	!-----------------------------------------------------------------------------------------------!
	!				Laminar viscosity						!
	!-----------------------------------------------------------------------------------------------!
		MuL_se = flux(Mu(i,j), Mu(i+1,j), DXU(i), DXP(i+1)/2.0D+0)
		MuL_ne = flux(Mu(i,j+1), Mu(i+1,j+1), DXU(i), DXP(i+1)/2.0D+0 )
		MuL_e = flux(MuL_se, MuL_ne, DYV(j), DYP(j+1)/2.0D+0)

		MuL_sw = flux(Mu(i-1,j), Mu(i,j), DXU(i-1), DXP(i)/2.0D+0)
		MuL_nw = flux(Mu(i-1,j+1),Mu(i,j+1),dxu(i-1), DXP(i)/2.0D+0)
		MuL_w = flux(MuL_sw,MuL_nw,DYV(j),DYP(j+1)/2.0D+0)		

	!-----------------------------------------------------------------------------------------------!
	!				Turbulent viscosity						!
	!-----------------------------------------------------------------------------------------------!
		MuT_se = flux(Mut(i,j), Mut(i+1,j), DXU(i), DXP(i+1)/2.0D+0)
		MuT_ne = flux(Mut(i,j+1), Mut(i+1,j+1), DXU(i), DXP(i+1)/2.0D+0 )
		MuT_e = flux(Mut_se, Mut_ne, DYV(j), DYP(j+1)/2.0D+0)

		MuT_sw = flux(Mut(i-1,j), Mut(i,j), DXU(i-1), DXP(i)/2.0D+0)
		MuT_nw = flux(Mut(i-1,j+1),Mut(i,j+1),dxu(i-1), DXP(i)/2.0D+0)
		MuT_w = flux(Mut_sw,Mut_nw,DYV(j),DYP(j+1)/2.0D+0)
	!-----------------------------------------------------------------------------------------------!
	!				Total viscosity							!
	!-----------------------------------------------------------------------------------------------!
		GIE = MuL_e + Mut_e
		GIw = MuL_w + Mut_w
		GIn = (MU(i,j+1)+MUt(i,j+1))
		GIs = (MU(i,j)+MUt(i,j))

		if (i .eq. Nx-1) GIE = flux( Mu(Nx,j)+MuT(Nx,j),Mu(Nx,j+1)+MuT(Nx,j+1),DYV(j),DYP(j+1)/2.0D+0 ) ! fluid boundary
		if (i .eq. 2) GIW = flux( Mu(1,j)+MuT(1,j),Mu(1,j+1)+MuT(1,j+1),DYV(j),DYP(j+1)/2.0D+0 ) 	! fluid boundary

! Nota: Cuando en laminar si impongo GI's = Muo converge mas rapido (200 iteraciones menos)


			De=( GIE * DYV(j) ) /  (X(i+1) - X(i)) 
		      
			Dw=( GIw * DYV(j) ) /  (X(i) - X(i-1)) 

			Dn=( GIn * DXP(i) ) / ( YV(j+1) - YV(j) )
		    
			Ds=( GIs * DXP(i) ) / ( YV(j) - YV(j-1) )



!-----------------------------------------------------------------------------------------------!
!					PECLET NUMBER						!
!-----------------------------------------------------------------------------------------------!


			Pee=Fe/De

			Pew=Fw/Dw

			Pen=Fn/Dn

			Pes=Fs/Ds


!-----------------------------------------------------------------------------------------------!
!					High_order_scheme				!
!-----------------------------------------------------------------------------------------------!
	sc_hos = 0.0D+0
	sp_hos = 0.0D+0
	Sdc =0.0D+0
!	if ( (i.gt. 2).and.(i.lt. Nx-1).and.(j.gt. 2).and.(j.lt. Ny-2).and.(n_high_order_scheme.ne. 0)   &
!		.and. high_order_on_velocities ) then
	if ( (n_high_order_scheme.ne. 0) .and. high_order_on_velocities_v ) then		

		call high_order_scheme_v( i,j,Fe,Fw,Fn,Fs,						&
			Vx(i-2,j),Vx(i-1,j),Vx(i,j),Vx(i+1,j),Vx(i+2,j),			&
			Vx(i,j-2),Vx(i,j-1),Vx(i,j+1),Vx(i,j+2),					&
			dxu(i-2),dxu(i-1),dxu(i),dxu(i+1),							&
			dyp(j-1),dyp(j),DYP(j+1),DYP(j+2),Sdc						)
		

		if ( (Sdc / (Vx(i,j)+1.0D-30)) .gt. 0.0D+0) then 

			sc_hos = Sdc
			sp_hos=0.0D+0

		else

			sc_hos=0.0D+0
			sp_hos = Sdc / (Vx(i,j)+1.0D-30)

		end if



	end if

		Sdc_ij_v(i,j) = Sdc

!-----------------------------------------------------------------------------------------------!
!					Source terms from stresses
!-----------------------------------------------------------------------------------------------!
			dvdx_e = ( Vx(i+1,j) - Vx(i,j)  ) / DXU(i)
			dvdx_w = ( Vx(i,j)   - Vx(i-1,j)) / DXU(i-1)	
	
			Source1 = ( TxyL(i,j)-TxyL(i-1,j)+TxyT(i,j)-TxyT(i-1,j) &
				- GIE*dvdx_e +  GIw*dvdx_w )*DYV(j)

			dvdy_n = ( Vx(i,j+1) - Vx(i,j)  ) / DYP(j+1)

			Source2 = ( TyyL(i,j+1)-TyyL(i,j)+TyyT(i,j+1)-TyyT(i,j) &
				- GIn*dvdy_n +  GIs*dvdy_s )*DXP(i)

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------!


			ae(i,j)=De*A(n,Pee)+DMAX1(-Fe,0.0D+00)

			aw(i,j)=Dw*A(n,Pew)+DMAX1(Fw,0.0D+00)

			an(i,j)=Dn*A(n,Pen)+DMAX1(-Fn,0.0D+00)

			as(i,j)=Ds*A(n,Pes)+DMAX1(Fs,0.0D+00)

			apo=(flux(rho(i,j),rho(i,j+1),DYV(j),DYP(j+1)/2.0D+00)*DXP(i)*DYV(j) ) /dt	

			apv(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+apo! - sp_hos
			

			if (Boussinesq.eq.1) then
								 
				b(i,j)=(Px(i,j)-Px(i,j+1))*DXP(i)+apo*Vx(i,j) + RHOO*g*BBETA* 	&
				(flux(T(i,j),T(i,j+1),DYV(j),(DYP(j+1))/2.0D+00)-Tref )*DXP(i)*DYV(j)
			else		
			
				b(i,j)=(Px(i,j)-Px(i,j+1))*DXP(i)+apo*Vx(i,j)-(RHO(i,j)-RHOo)*g*DXP(i)*DYV(j)
			end if

			b(i,j) = b(i,j)! + Source1 + Source2 + sc_hos



		end do
	end do 



!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!
     

!C	Frontera sur
	Do  i=1,Nx
		ae(i,1)		=	0.0D+00
		aw(i,1)		=	0.0D+00
		apv(i,1)	=	1.0D+00
		an(i,1)		=	0.0D+00
		as(i,1)		=	0.0D+00
		b(i,1)		=	0.0D+00
	end do
 

 !C	Frontera norte
  	Do  i=1,Nx	
   			ae(i,Ny-1)	=	0.0D+00
   			aw(i,Ny-1)	=	0.0D+00
   			apv(i,Ny-1)	=	1.0D+00
   			an(i,Ny-1)	=	0.0D+00
   			as(i,Ny-1)	=	0.0D+00 
   			b(i,Ny-1)	=	0.0D+00

  	end do

	
!					WEST					!
	Do j=2,Ny-2
      
		ae(1,j)		=	0.0D+00
		aw(1,j)		=	0.0D+00
		apv(1,j)	=	1.0D+00
		an(1,j)		=	0.0D+00
		as(1,j)		=	0.0D+00
		b(1,j)		=	0.0D+00
  
	end do


!					EAST					!
	Do  j=2,Ny-2

		ae(Nx,j)	=	0.0D+00
		aw(Nx,j)	=	0.0D+00
		apv(Nx,j)	=	1.0D+00
		an(Nx,j)	=	0.0D+00
		as(Nx,j)	=	0.0D+00
		b(Nx,j)		=	0.0D+00
  
	end do

!-----------------------------------------------------------------------------------------------!
!!			d's Matrix determination for SIMPLE algorithm				!
!-----------------------------------------------------------------------------------------------!   
	do j=1,Ny-1
		do i=1,Nx

      		IF (ALGORITMO.EQ.1) APVNB(I,J)=0.0D+00      			!ALGORITMO SIMPLE
      		IF (ALGORITMO.EQ.2) APVNB(I,J)=AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J) 	!ALGORITMO SIMPLEC

		dv(i,j)=DXP(i)/(apv(i,j)-APVNB(i,j))

		end do
	end do


!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!

	Resi=0.0D+0
	summ=0.0D+00
	Rmax_V=0.0D+00

if (street_canyon .eqv. .false.) then 

	do j=2,Ny-2
		do  i=2,Nx-1

			Resi=apv(i,j)*Vx(i,j)-ae(i,j)*Vx(i+1,j)-aw(i,j)*Vx(i-1,j)-as(i,j)*Vx(i,j-1)-an(i,j)*Vx(i,j+1)-b(i,j)
			summ=summ+(Resi*Resi)
			Rij_v(i,j) = Resi
			if (Resi .gt. Rmax_V) Rmax_V=Resi
 
 
		end do
	end do

else

	do j=2,Ny-2
	do i=2,Nx-1
		outside_canyon =(X(i).ge.((BD+BW)-1.0D+00)).and.(Y(j).le.B2+1.0D+00).or.&
		(X(i).le.BW+1.0D+00).and.(Y(j).le.B1+1.0D+00)
			if (outside_canyon) then
			else
				Resi=apv(i,j)*Vx(i,j)-ae(i,j)*Vx(i+1,j)-aw(i,j)*Vx(i-1,j)-as(i,j)*Vx(i,j-1)-an(i,j)*Vx(i,j+1)-b(i,j)
				summ=summ+(Resi*Resi)
				Rij_v(i,j) = Resi
				if (Resi .gt. Rmax_V) Rmax_V=Resi
			end if
	end do
	end do

end if

	Res_V=dsqrt(summ)
	
	return

end Subroutine coef_V





!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					SUBROTINA COEF P					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!



Subroutine COEF_P(DXP,DYP,Ux,Vx,du,dv,ap,ae,aw,as,an,b,RES_P,Rmax_P,RHO,apu,apv,apunb,apvnb)	


	implicit none

	integer:: i,j,npas

	real(KIND = DP)  :: RES_P,De,Dw,Dn,Ds,resi,summ,Rmax_P,RHOe,RHOw,RHOs,RHOn,P_f

	real(KIND = DP), dimension(Nx) :: 	&
		DXP

	real(KIND = DP), dimension(Ny) :: 	&
		DYP

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		ap,ae,aw,an,as,b,		&
		Ux,Vx,Px,du,dv,			&
		RHO,apu,apv,apunb,apvnb

																	 
																	 
	do j=1,Ny
		do i=1,Nx							  
																	 
			ae(i,j)=0.0D+00
																			 
			aw(i,j)=0.0D+00
																			 
			an(i,j)=0.0D+00
																			 
			as(i,j)=0.0D+00
																			 
			ap(i,j)=0.0D+00
																			 
			b(i,j)=0.0D+00	
																	 
		end do
	end do 



	do  j=2,Ny-1
		do i=2,Nx-1



!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				!
!-----------------------------------------------------------------------------------------------


			RHOe=flux(RHO(i,j),RHO(i+1,j),DXU(i),DXP(i+1)/2.0D+00)
			RHOw=flux(rho(i-1,j),rho(i,j),DXU(i-1),DXP(i)/2.0D+00)
			RHOn=flux(rho(i,j),rho(i,j+1),DYV(j),DYP(j+1)/2.0D+00)
			RHOs=flux(rho(i,j-1),rho(i,j),DYV(j-1),DYP(j)/2.0D+00)


			ae(i,j)=RHOe*du(i,j)*DYP(j)

			aw(i,j)=RHOw*du(i-1,j)*DYP(j)

			an(i,j)=RHOn*dv(i,j)*DXP(i)

			as(i,j)=RHOs*dv(i,j-1)*DXP(i)
	
			IF (I.EQ.NX-1)	AE(I,J) = 0.0D+00	!"C.F.=LADO ESTE"
			IF (I.EQ.2)		AW(I,J) = 0.0D+00       !"C.F.=LADO OESTE"
			IF (J.EQ.NY-1)	AN(I,J) = 0.0D+00   !"C.F.=LADO NORTE"
			IF (J.EQ.2)		AS(I,J) = 0.0D+00       !"C.F.=LADO SUR"

			ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)
	
			b(i,j)=RHOw*DYP(j)*Ux(i-1,j) 		&
			      -RHOe*DYP(j)*Ux(i,j)  		&
			      +RHOs*DXP(i)*Vx(i,j-1)   		&
			      -RHOn*DXP(i)*Vx(i,j)             
			      !+(RHO(i,j)-RHO(i,j) )*DXP(i)*DYP(j)/dt	  


		end do
	end do


!-----------------------------------------------------------------------------------------------!
!		Force a reference point of Pc=0 in the domain					!
!-----------------------------------------------------------------------------------------------! 


!		ae(Nx-1,Ny-1)=0.0D+00
!		aw(Nx-1,Ny-1)=0.0D+00
!		ap(Nx-1,Ny-1)=1.0D+00
!		an(Nx-1,Ny-1)=0.0D+00
!		as(Nx-1,Ny-1)=0.0D+00
!		b (Nx-1,Ny-1)=0.0D+00

! 		ae(2,Ny-1)=0.0D+00
! 		aw(2,Ny-1)=0.0D+00
! 		ap(2,Ny-1)=1.0D+00
! 		an(2,Ny-1)=0.0D+00
! 		as(2,Ny-1)=0.0D+00
! 		b (2,Ny-1)=0.0D+00



! Nota: si activo el punto de referencia (Nx-1,Ny-1) tarda unas 35 iteraciones mas en converger

	!	ae(2,2)=0.0D+00
	!	aw(2,2)=0.0D+00
	!	ap(2,2)=1.0D+00
	!	an(2,2)=0.0D+00
	!	as(2,2)=0.0D+00
	!	b(2,2)=0.0D+00

! Nota: si activo el punto de referencia (2,2) se reduce a la mitad el num. de iteraciones pero da una solucion falsa!
!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS					!
!-----------------------------------------------------------------------------------------------!    

!	Frontera sur
	Do  i=1,Nx

		ae(i,1)=0.0D+00
		aw(i,1)=0.0D+00
		ap(i,1)=1.0D+00
		an(i,1)=0.0D+00
		as(i,1)=0.0D+00
		b(i,1)=0.0D+00
	end do
   
 
!	Frontera norte

	Do i=1,Nx

		ae(i,Ny)=0.0D+00
		aw(i,Ny)=0.0D+00
		ap(i,Ny)=1.0D+00
		an(i,Ny)=0.0D+00
		as(i,Ny)=0.0D+00
		b(i,Ny)=0.0D+00
	end do

 

!     Frontera Oeste
	Do j=2,Ny-1
      
		ae(1,j)=0.0D+00
		aw(1,j)=0.0D+00
		ap(1,j)=1.0D+00
		an(1,j)=0.0D+00
		as(1,j)=0.0D+00
		b(1,j)=0.0D+00
  
	end do


!     Frontera Este
      Do j=2,Ny-1
  
		ae(Nx,j)=0.0D+00
		aw(Nx,j)=0.0D+00
		ap(Nx,j)=1.0D+00
		an(Nx,j)=0.0D+00
		as(Nx,j)=0.0D+00
		b(Nx,j)=0.0D+00
  
	end do

   
      
!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!


	summ=0.0D+00
	Rmax_P=0.0D+00

if (street_canyon .eqv. .false.) then 

	do j=2,Ny-1
		do i=2,Nx-1

			summ=summ+DABS(b(i,j))
			if (DABS(b(i,j)) .gt. Rmax_P) Rmax_P=DABS(b(i,j))  
		end do
	end do

else

		do j=2,Ny-1
			do i=2,Nx-1
			outside_canyon =(X(i).ge.((BD+BW)-1.0D+00)).and.(Y(j).le.B2+1.0D+00).or.&
		(X(i).le.BW+1.0D+00).and.(Y(j).le.B1+1.0D+00)
					if (outside_canyon) then
					else
						summ=summ+DABS(b(i,j))
						if (DABS(b(i,j)) .gt. Rmax_P) Rmax_P=DABS(b(i,j))
					end if
			end do
		end do

end if

	Res_P=summ
	
	return

end Subroutine COEF_P



!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						reset(Pc)					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine reset(Pc)	
	
	real(KIND = DP), dimension(Nx,Ny) :: Pc
	integer :: i,j

     	Do j=1,Ny
		Do i=1,Nx

			Pc(i,j)=0.0D+0
		end do
	end do
	
	return
end Subroutine reset

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						converegence					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function converegence(Res_U,Res_V,Res_P,Res_T,Res_C,Res_e,Res_TK,iter) 


	real(KIND = DP):: Res_U,Res_V,Res_P,Res_T,Res_C,Res_e,Res_TK,Rmax
	logical:: converegence
	integer:: iter	       

	if (iter.le.20) then
	Res_U=1.0D+00
	end if

	Rmax = DMAX1(Res_U,Res_V,Res_P,Res_T,Res_C,Res_e,Res_TK)
	if (Rmax .gt. 1.0D+20) stop 'A friendly message: diverged'
	
	converegence=(( (Res_U.LE.epsilonU).and.(Res_V.LE.epsilonV)).and.(Res_P.LE.epsilonP).and.(Res_T.LE.epsilonT) &
			.and.(Res_C.LE.epsilonC).and.(Res_TK.LE.epsilonTK).and.(Res_e.LE.epsilone))
	
	return     
end 



!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Mass_flux				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine Mass_Flux(U,V,RHO,Fep,Fnp)	

	implicit none
	integer:: i,j
	real(KIND = DP), dimension(Nx,Ny),intent(in) :: U,V,RHO
	real(KIND = DP), dimension(Nx,Ny),intent(out) :: Fep,Fnp
	real(KIND = DP)  :: RHOe,RHOn




	Do j=1,Ny
		Do i=2,(Nx-2)

			RHOe=flux(RHO(i,j),RHO(i+1,j),DXU(i),DXP(i+1)/2.0D+00)	

			Fep(i,j) = U(i,j)*RHOe!( RHO(i,j)+RHO(i+1,j) ) / 2.0D+0
		end do

		Fep(1,j)    = U(1,j)   * RHO(1,j)
		Fep(Nx-1,j) = U(Nx-1,j)* RHO(Nx,j)
	end do

	Do i=1,Nx
		Do j=2,(Ny-2)

			RHOn=flux(rho(i,j),rho(i,j+1),DYV(j),DYP(j+1)/2.0D+00)	

			Fnp(i,j) = V(i,j)*( RHO(i,j)+RHO(i,j+1) ) / 2.0D+0

		end do

		Fnp(i,1)    = V(i,1)   * RHO(i,1)
		!Fnp(i,1)    = 0.0D+0 ! Solid boundary
		Fnp(i,Ny-1) = V(i,Ny-1)* RHO(i,Ny)
		!Fnp(i,Ny-1) = 0.0D+0 ! Solid boundary

		!Fep(i,1)    =  U(i,1)   * RHO(i,1) 
		!Fep(i,Ny)    = U(i,Ny)  * RHO(i,Ny)
	end do 


    return
end subroutine Mass_Flux


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						Avarege						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

	function flux(Tp,Te,dx,dxmas)

	!	Calculo del valor de una propiedad en la frontera asumiendo variación linial
	! 	dxmas: distancia de la frontera al punto nodal Te
	!	dx: distancia entre puntos nodales

		real(KIND = DP):: Tp,Te,dx,dxmas,fe,flux	       

		fe=dxmas/dx
		flux=( (1.0D+00-fe)*Te )+( fe*Tp )
	
		return     
      end 


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					Harmonic mean						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

	function harmo(CP,CE,dx,dxmas )

	!	Calculo de la media harmonica una propiedad en la frontera
	!	dxmas: distancia de la frontera al punto nodal Ce
	!	dx: distancia entre puntos nodales

	       
		real(KIND = DP):: CP,CE,dx,dxmas,fe,harmo

		fe=dxmas/dx
		harmo=1.0D+00/( ((1.0D+00-fe)/CP) + (fe/CE) )
	
		return     
      
	end 






!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						Properites					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine PROPERTIES(Mu,CONTER,CP,CONTERo,T,Muo,CPo,RHO,Po,Pmed,Hx,Hy,DXP,DYP,Tref,RHOo,MuT,e,TK)	

	real(KIND = DP):: MA0,MA1,MA2,MA3,MA4, KA0,KA1,KA2,KA3,KA4,KA5, CA0,CA1,CA2,CA3,CA4,CA5, &
			  Muo,contero,cpo,Tp, RHOmed, summ, Pmed ,Po,Hx,Hy ,Tref,RHOo, f_mu_p, Re_t

	real(KIND = DP), dimension(Nx) :: DXP
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: Mu,CONTER,CP,T,RHO,MuT,e,TK
	integer :: i,j,k


!-----------------------------------------------------------------------------------------------!
!						Viscocity					!
!-----------------------------------------------------------------------------------------------! 
!	MA0=-9.8601D-1
!	MA1=9.080125D-2
!	MA2=-1.17635575D-4
!	MA3=1.2349703D-7
!	MA4=-5.7971299D-11

	Do j=1,Ny
		Do i=1,Nx
			if (CttProperties.eq.1) then			
				Mu(i,j)=Muo 			!Constant value
			else
				Tp=T(i,j)		
				!Mu(i,j)= (-4.6763d-12  * Tp * Tp   &   !Water vapor
                !         +4.8259d-8   * Tp        &
                !         -5.9340d-6)

				Mu(i,j)= 1.68D-5*( (Tp/273.0D+0)**(3.0D+0/2.0D+0) )*( 383.5D+0/(Tp+110.5) ) !SUTHERLAND's Law	
			end if

			if (street_canyon) then
				outside_canyon =(X(i).ge.(BD+BW)).and.(Y(j).le.B2).or.&
		(X(i).le.BW).and.(Y(j).le.B1)
				if (outside_canyon) Mu(i,j)=1.0D+10
				if (outside_canyon) CONTER(i,j)=1.0D+10
			end if


		end do
	end do

!-----------------------------------------------------------------------------------------------!
!				Turbulent Viscocity (MuT)					!
!-----------------------------------------------------------------------------------------------! 
	Do j=1,Ny
		Do i=1,Nx		
			MuT(i,j)= MUT_model(i,j,Tk,e,rho,mu)
			!write(*,*) MuT(i,j)
			if ( (MuT(i,j).gt.1.0D+4*MUo).and.(iter.ne.1) ) MuT(i,j) = 1.0D+4*MUo
		end do
	end do
	
!	BC Turbulent Viscocity (MuT)					!
!	Do i = 1,Nx
!		MuT(i,1)  = 0.0D+0 ! South boundary is a Solid 
!		MuT(i,Ny) = 0.0D+0 ! North boundary is a Solid boundary
!	end do
!	Do j = 1,Ny
		!MuT(1,j)  = 0.0D+0 ! inlet Fluid Boundary
		!MuT(Nx,j) = 0.0D+0 ! outlet Fluid Boundary
!	end do

!-----------------------------------------------------------------------------------------------!
!					Thermal conductivity					!
!-----------------------------------------------------------------------------------------------! 
!	KA0=-2.276501D-3
!	KA1=1.2598485D-4
!	KA2=-1.4815235D-7
!	KA3=1.73550646D-10
!	KA4=-1.066657D-13
!	KA5=2.47663035D-17

	Do j=1,Ny
		Do i=1,Nx
			if (CttProperties.eq.1) then	
				CONTER(i,j)=CONTERo  			!Constant value
			else
				Tp=T(i,j)
				CONTER(i,j)=     (+5.3644d-8 * Tp * Tp   & 
                           		 +3.9243d-5 * Tp        &
                            	 +5.0855d-3)
                !CONTER(i,j)=CONTERo             	 
!				CONTER(i,j)=Mu(i,j)*401.8D+0/2.84D-1
			end if
		end do
	end do
!-----------------------------------------------------------------------------------------------!
!					Heat capacity						!
!-----------------------------------------------------------------------------------------------! 
	CA0= 0.103409D+01
	CA1=-0.284887D-3
	CA2= 0.7816818D-6
	CA3=-0.4970786D-9
	CA4= 0.1077024D-12

	Do j=1,Ny
		Do i=1,Nx
			if (CttProperties.eq.1) then	
				Cp(i,j)=CPo
			else
				Tp=T(i,j)
				Cp(i,j)= (CA0 + CA1*Tp + CA2*Tp*Tp + CA3*Tp*Tp*Tp + CA4*Tp*Tp*Tp*Tp)*1.0D+3
			end if
		end do
	end do

!-----------------------------------------------------------------------------------------------!
!						Density						!
!-----------------------------------------------------------------------------------------------! 

!write(*,*) 
!stop

	if (Boussinesq.eq.2) then	
		summ=0.0D+00
		Do j=2,Ny-1
		Do i=2,Nx-1
			summ=summ+( DXP(i)*DYP(j)/T(i,j) )
		end do
		end do
	
		Pmed=Po*(Hx*Hx/Tref)/summ
	end if
	Pmed = Po

	Do j=1,Ny
	Do i=1,Nx
		if (Boussinesq.eq.1) then
			RHO(i,j)=RHOo
		else 		
			RHO(i,j)=Pmed/(R*T(i,j))
		end if

	end do
	end do


	return
end Subroutine PROPERTIES

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						correct						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine correct (U,Ux,V,Vx,Pc,P,Px,du,dv,Fep,Fnp,rho)	
	

	integer :: i,j	
	real(KIND = DP), dimension(Nx,Ny),intent(in) :: rho,Ux,Vx,Pc,Px,dv,du
	real(KIND = DP), dimension(Nx,Ny),intent(out) :: Fep,Fnp,U,V,P

    Do j=2,Ny-1
		Do i=2,Nx-2
			U(i,j)=Ux(i,j)+du(i,j)*( Pc(i,j)-Pc(i+1,j) )
		end do
	end do

	Do j=2,(Ny-2)
		Do i=2,(Nx-1)
			V(i,j)=Vx(i,j)+dv(i,j)*( Pc(i,j)-Pc(i,j+1) )
		end do
	end do

	Do j=2,Ny-1
		Do i=2,Nx-1
			P(i,j)=Px(i,j)+relaxP*Pc(i,j)
		end do
	end do 

	call Mass_Flux(U,V,RHO,Fep,Fnp)	

	return
end Subroutine correct


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						 A(n,Pe)					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function A(n,Pe) 
	
	real(KIND = DP)  :: Pe,A
	integer :: n


	if (n.EQ.1) A=1.0D+00-0.5D+00*DABS(Pe)
	if (n.EQ.2) A=1.0D+00
	if (n.EQ.3) A=DMAX1( 0.0D+00,1.0D+00-( 5.0D-01*DABS(Pe)) ) 
	if (n.EQ.4) A=DMAX1(0.0D+00,(1.0D+00-( 1.0D-01*DABS(Pe) ))**5.0D+0)
	if (n.EQ.5) A=DABS(Pe)/(DEXP(ABS(Pe))-1.0D+00)

	return
end


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					high_order_scheme					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine high_order_scheme( i,j,Fe,Fw,Fn,Fs,					&
			Tww,Tw,Tp,Te,Tee,					&
			Tss,Ts,Tn,Tnn,						&
			dw_ww,dp_w,de_p,dee_e,					&
			ds_ss,dp_s,dn_p,dnn_n,Sdc				)	

	implicit none

	real(KIND = DP):: Fe,Fw,Fn,Fs,							&
			Tww,Tw,Tp,Te,Tee,						&
			Tss,Ts,Tn,Tnn,							&
			dw_ww,dp_w,de_p,dee_e,						&
			ds_ss,dp_s,dn_p,dnn_n,Sdc,					&
			alpha_e,alpha_w,alpha_s,alpha_n,				&
			re_p,re_n,rw_p,rw_n,rs_p,rs_n,rn_p,rn_n,			&
			be, bw, bn, bs

	integer :: i,j

	alpha_e=1.0D+0
	alpha_w=1.0D+0
	alpha_s=1.0D+0
	alpha_n=1.0D+0

	if (Fe .lt. 0.0D+0) alpha_e=0.0D+0
	if (Fw .lt. 0.0D+0) alpha_w=0.0D+0
	if (Fn .lt. 0.0D+0) alpha_n=0.0D+0	
	if (Fs .lt. 0.0D+0) alpha_s=0.0D+0


	re_p = ( (Tp-Tw) )   / ( (Te-Tp)+1.0D-30 ) 
	re_n = ( (Tee-Te) ) / ( (Te-Tp)+1.0D-30 )	

	rw_p = ( (Tw-Tww) ) / ( (Tp-Tw)+1.0D-30 )
	rw_n = ( (Te-Tp))   / ( (Tp-Tw)+1.0D-30 )

	rn_p = ( (Tp-Ts) )   / ( (Tn-Tp)+1.0D-30 )
	rn_n = ( (Tnn-Tn) ) / ( (Tn-Tp)+1.0D-30 )

	rs_p = ( (Ts-Tss) ) / ( (Tp-Ts)+1.0D-30 )
	rs_n = ( (Tn-Tp) )   / ( (Tp-Ts)+1.0D-30 )

!	re_p = ( (Tp-Tw)/dp_w )   / ( (Te-Tp)/de_p ) 
!	re_n = ( (Tee-Te)/dee_e ) / ( (Te-Tp)/de_p )	

!	rw_p = ( (Tw-Tww)/dw_ww ) / ( (Tp-Tw)/dp_w )
!	rw_n = ( (Te-Tp)/de_p )   / ( (Tp-Tw)/dp_w )

!	rn_p = ( (Tp-Ts)/dp_s )   / ( (Tn-Tp)/dn_p )
!	rn_n = ( (Tnn-Tn)/dnn_n ) / ( (Tn-Tp)/dn_p )

!	rs_p = ( (Ts-Tss)/ds_ss ) / ( (Tp-Ts)/dp_s )
!	rs_n = ( (Tn-Tp)/dn_p )   / ( (Tp-Ts)/dp_s )

	be=1.0D+0
	bw=1.0D+0
	bn=1.0D+0
	bs=1.0D+0

	if (i.eq. 3) 	bw = 0.0D+0
	if (i.eq. nx-2) be = 0.0D+0
	if (j.eq. 3) 	bs = 0.0D+0
	if (j.eq. ny-2) bn = 0.0D+0 

	Sdc =	0.5D+0 * be  * Fe * ( (1.0D+0-alpha_e) * phsi(re_n) - alpha_e          * phsi(re_p) ) * (Te-Tp)	&
	    +	0.5D+0 * bw  * Fw * (          alpha_w * phsi(rw_p) - (1.0D+0-alpha_w) * phsi(rw_n) ) * (Tp-Tw)	&
	    +	0.5D+0 * bn  * Fn * ( (1.0D+0-alpha_n) * phsi(rn_n) - alpha_n          * phsi(rn_p) ) * (Tn-Tp)	&
	    +	0.5D+0 * bs  * Fs * (          alpha_s * phsi(rs_p) - (1.0D+0-alpha_s) * phsi(rs_n) ) * (Tp-Ts)	



	return
end Subroutine  high_order_scheme


Subroutine high_order_scheme_u( i,j,Fe,Fw,Fn,Fs,					&
			Tww,Tw,Tp,Te,Tee,					&
			Tss,Ts,Tn,Tnn,						&
			dw_ww,dp_w,de_p,dee_e,					&
			ds_ss,dp_s,dn_p,dnn_n,Sdc				)	

	implicit none

	real(KIND = DP):: Fe,Fw,Fn,Fs,							&
			Tww,Tw,Tp,Te,Tee,						&
			Tss,Ts,Tn,Tnn,							&
			dw_ww,dp_w,de_p,dee_e,						&
			ds_ss,dp_s,dn_p,dnn_n,Sdc,					&
			alpha_e,alpha_w,alpha_s,alpha_n,				&
			re_p,re_n,rw_p,rw_n,rs_p,rs_n,rn_p,rn_n,			&
			be, bw, bn, bs

	integer :: i,j

	alpha_e=1.0D+0
	alpha_w=1.0D+0
	alpha_s=1.0D+0
	alpha_n=1.0D+0

	if (Fe .lt. 0.0D+0) alpha_e=0.0D+0
	if (Fw .lt. 0.0D+0) alpha_w=0.0D+0
	if (Fn .lt. 0.0D+0) alpha_n=0.0D+0	
	if (Fs .lt. 0.0D+0) alpha_s=0.0D+0


	re_p = ( (Tp-Tw) )   / ( (Te-Tp)+1.0D-30 ) 
	re_n = ( (Tee-Te) ) / ( (Te-Tp)+1.0D-30 )	

	rw_p = ( (Tw-Tww) ) / ( (Tp-Tw)+1.0D-30 )
	rw_n = ( (Te-Tp))   / ( (Tp-Tw)+1.0D-30 )

	rn_p = ( (Tp-Ts) )   / ( (Tn-Tp)+1.0D-30 )
	rn_n = ( (Tnn-Tn) ) / ( (Tn-Tp)+1.0D-30 )

	rs_p = ( (Ts-Tss) ) / ( (Tp-Ts)+1.0D-30 )
	rs_n = ( (Tn-Tp) )   / ( (Tp-Ts)+1.0D-30 )

! 	re_p = ( (Tp-Tw)/dp_w )   / ( (Te-Tp)/de_p ) 
! 	re_n = ( (Tee-Te)/dee_e ) / ( (Te-Tp)/de_p )	

! 	rw_p = ( (Tw-Tww)/dw_ww ) / ( (Tp-Tw)/dp_w )
! 	rw_n = ( (Te-Tp)/de_p )   / ( (Tp-Tw)/dp_w )

! 	rn_p = ( (Tp-Ts)/dp_s )   / ( (Tn-Tp)/dn_p )
! 	rn_n = ( (Tnn-Tn)/dnn_n ) / ( (Tn-Tp)/dn_p )

! 	rs_p = ( (Ts-Tss)/ds_ss ) / ( (Tp-Ts)/dp_s )
! 	rs_n = ( (Tn-Tp)/dn_p )   / ( (Tp-Ts)/dp_s )

	be=1.0D+0
	bw=1.0D+0
	bn=1.0D+0
	bs=1.0D+0

	if (i.eq. 3) 	bw = 0.0D+0
	if (i.eq. nx-3) be = 0.0D+0
	if (j.eq. 3) 	bs = 0.0D+0
	if (j.eq. ny-2) bn = 0.0D+0 

	Sdc =	0.5D+0 * be  * Fe * ( (1.0D+0-alpha_e) * phsi(re_n) - alpha_e          * phsi(re_p) ) * (Te-Tp)	&
	    +	0.5D+0 * bw  * Fw * (          alpha_w * phsi(rw_p) - (1.0D+0-alpha_w) * phsi(rw_n) ) * (Tp-Tw)	&
	    +	0.5D+0 * bn  * Fn * ( (1.0D+0-alpha_n) * phsi(rn_n) - alpha_n          * phsi(rn_p) ) * (Tn-Tp)	&
	    +	0.5D+0 * bs  * Fs * (          alpha_s * phsi(rs_p) - (1.0D+0-alpha_s) * phsi(rs_n) ) * (Tp-Ts)	

	if (i.le. 2) 	Sdc = 0.0D+0
	if (i.ge. nx-2) Sdc = 0.0D+0
	if (j.le. 2) 	Sdc = 0.0D+0
	if (j.ge. ny-1) Sdc = 0.0D+0


	return
end Subroutine  high_order_scheme_u

Subroutine high_order_scheme_v( i,j,Fe,Fw,Fn,Fs,					&
			Tww,Tw,Tp,Te,Tee,					&
			Tss,Ts,Tn,Tnn,						&
			dw_ww,dp_w,de_p,dee_e,					&
			ds_ss,dp_s,dn_p,dnn_n,Sdc				)	

	implicit none

	real(KIND = DP):: Fe,Fw,Fn,Fs,							&
			Tww,Tw,Tp,Te,Tee,						&
			Tss,Ts,Tn,Tnn,							&
			dw_ww,dp_w,de_p,dee_e,						&
			ds_ss,dp_s,dn_p,dnn_n,Sdc,					&
			alpha_e,alpha_w,alpha_s,alpha_n,				&
			re_p,re_n,rw_p,rw_n,rs_p,rs_n,rn_p,rn_n,			&
			be, bw, bn, bs

	integer :: i,j

	alpha_e=1.0D+0
	alpha_w=1.0D+0
	alpha_s=1.0D+0
	alpha_n=1.0D+0

	if (Fe .le. 0.0D+0) alpha_e=0.0D+0
	if (Fw .le. 0.0D+0) alpha_w=0.0D+0
	if (Fn .le. 0.0D+0) alpha_n=0.0D+0	
	if (Fs .le. 0.0D+0) alpha_s=0.0D+0


	re_p = ( (Tp-Tw) )   / ( (Te-Tp)+1.0D-25 ) 
	re_n = ( (Tee-Te) ) / ( (Te-Tp)+1.0D-25 )	

	rw_p = ( (Tw-Tww) ) / ( (Tp-Tw)+1.0D-25 )
	rw_n = ( (Te-Tp))   / ( (Tp-Tw)+1.0D-25 )

	rn_p = ( (Tp-Ts) )   / ( (Tn-Tp)+1.0D-25 )
	rn_n = ( (Tnn-Tn) ) / ( (Tn-Tp)+1.0D-25 )

	rs_p = ( (Ts-Tss) ) / ( (Tp-Ts)+1.0D-25 )
	rs_n = ( (Tn-Tp) )   / ( (Tp-Ts)+1.0D-25 )

! 	re_p = ( (Tp-Tw)/dp_w )   / ( (Te-Tp)/de_p ) 
! 	re_n = ( (Tee-Te)/dee_e ) / ( (Te-Tp)/de_p )	

! 	rw_p = ( (Tw-Tww)/dw_ww ) / ( (Tp-Tw)/dp_w )
! 	rw_n = ( (Te-Tp)/de_p )   / ( (Tp-Tw)/dp_w )

! 	rn_p = ( (Tp-Ts)/dp_s )   / ( (Tn-Tp)/dn_p )
! 	rn_n = ( (Tnn-Tn)/dnn_n ) / ( (Tn-Tp)/dn_p )

! 	rs_p = ( (Ts-Tss)/ds_ss ) / ( (Tp-Ts)/dp_s )
! 	rs_n = ( (Tn-Tp)/dn_p )   / ( (Tp-Ts)/dp_s )

	be=1.0D+0
	bw=1.0D+0
	bn=1.0D+0
	bs=1.0D+0

	if (i.eq. 3) 	bw = 0.0D+0
	if (i.eq. nx-2) be = 0.0D+0
	if (j.eq. 3) 	bs = 0.0D+0
	if (j.eq. ny-3) bn = 0.0D+0 

	Sdc =	0.5D+0 * be  * Fe * ( (1.0D+0-alpha_e) * phsi(re_n) - alpha_e          * phsi(re_p) ) * (Te-Tp)	&
	    +	0.5D+0 * bw  * Fw * (          alpha_w * phsi(rw_p) - (1.0D+0-alpha_w) * phsi(rw_n) ) * (Tp-Tw)	&
	    +	0.5D+0 * bn  * Fn * ( (1.0D+0-alpha_n) * phsi(rn_n) - alpha_n          * phsi(rn_p) ) * (Tn-Tp)	&
	    +	0.5D+0 * bs  * Fs * (          alpha_s * phsi(rs_p) - (1.0D+0-alpha_s) * phsi(rs_n) ) * (Tp-Ts)	

	if (i.le. 2) 	Sdc = 0.0D+0
	if (i.ge. nx-1) Sdc = 0.0D+0
	if (j.le. 2) 	Sdc = 0.0D+0
	if (j.ge. ny-2) Sdc = 0.0D+0 

	return
end Subroutine  high_order_scheme_v
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!											 phsi(r)											!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function phsi(r) 
	
	real(KIND = DP)  :: phsi,r

! 2 not converging
! 4 not converging
! 5 not converging

	if (n_high_order_scheme.EQ.1) phsi = (r + abs(r)) / (1.0D+0 + r) ! Van Leer

	if (n_high_order_scheme.EQ.2) phsi = (r + r*r) / (1.0D+0 + r*r ) ! Van Albada

	if (n_high_order_scheme.EQ.3) then				! Min-Mod
		if (r .gt. 0.0D+0 ) phsi = min(r,1.0D+0) 
		if (r .le. 0.0D+0 ) phsi = 0.0D+0
	end if

	if (n_high_order_scheme.EQ.4) phsi = max( 0.0D+0, min(2.0D+0*r,1.0D+0), min(r,2.0D+0) ) ! SUPERBEE

	if (n_high_order_scheme.EQ.5) phsi = max( 0.0D+0, min( 2.0D+0*r , (3.0D+0+r)/4.0D+0 ,2.0D+0) ) ! QUICK

	if (n_high_order_scheme.EQ.6) phsi = max( 0.0D+0, min( 2.0D+0*r , (1.0D+0 + 3.0D+0*r)/4.0D+0, &
							(3.0D+0+r)/4.0D+0 ,2.0D+0) ) !UMIST

	return
end

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			SUBROUTINE LBL X(I0,J0,I1,J1,AP,AE,AW,AN,AS,B,FHI,NPAS,RELAX)		!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!


SUBROUTINE LBL_X (I0,J0,I1,J1,AP,AE,AW,AN,AS,B,FHI,NPAS,RELAX)                 
   
    
	implicit none

	integer:: i2,i1,i0,j2,j1,j0,iCUENTA,i,j,npas

	real(KIND = DP)  :: beta,denom,relax

	real(KIND = DP), dimension(Nx,Ny) :: FHI,FHIOLD,ae,aw,an,as,ap,b
	real(KIND = DP), dimension(Nx) :: pt,qt


	I2=I1-1
	ICUENTA=0

	321 CONTINUE

	DO  I=I0,I1              
		DO  J=J0,J1 
                      FHIOLD(I,J)=FHI(I,J)
	  	end do
	end do

	DO J=J0,J1             

		BETA=B(I0,J)+AS(I0,J)*FHIOLD(I0,J-1)+AN(I0,J)*FHIOLD(I0,J+1)

		IF(J.EQ.J0)   BETA=B(I0,J)+AN(I0,J)*FHIOLD(I0,J+1)   
		IF(J.EQ.J1)   BETA=B(I0,J)+AS(I0,J)*FHIOLD(I0,J-1) 
	      
		DENOM=AP(I0,J)
		PT(I0)=AE(I0,J)/DENOM
		QT(I0)=BETA/DENOM

		DO I=I0+1,I2
			BETA=B(I,J)+AS(I,J)*FHIOLD(I,J-1)+AN(I,J)*FHIOLD(I,J+1)

	      		IF(J.EQ.J0)      BETA=B(I,J)+AN(I,J)*FHIOLD(I,J+1)    
	      		IF(J.EQ.J1)      BETA=B(I,J)+AS(I,J)*FHIOLD(I,J-1)   

	      		DENOM=AP(I,J)-PT(I-1)*AW(I,J)
	      		PT(I)=AE(I,J)/DENOM
	     		QT(I)=(BETA+QT(I-1)*AW(I,J))/DENOM
	  
		end do

		BETA=B(I1,J)+AS(I1,J)*FHIOLD(I1,J-1)+AN(I1,J)*FHIOLD(I1,J+1)

		IF(J.EQ.J0)   BETA=B(I,J1)+AN(I1,J)*FHIOLD(I1,J+1)   
		IF(J.EQ.J1)   BETA=B(I,J1)+AS(I1,J)*FHIOLD(I1,J-1)  
		IF(J.EQ.J0)   BETA=B(I1,J)+AN(I1,J)*FHIOLD(I1,J+1)   
		IF(J.EQ.J1)   BETA=B(I1,J)+AS(I1,J)*FHIOLD(I1,J-1)  

		DENOM=AP(I1,J)-PT(I2)*AW(I1,J)
		QT(I1)=(BETA+QT(I2)*AW(I1,J))/DENOM
     
		FHI(I1,J)=QT(I1)

		DO I=I2,I0,-1
			FHI(I,J)=PT(I)*FHI(I+1,J)+QT(I)
		end do

	end do   

	DO  I=I0,I1              
		DO J=J0,J1 
			FHI(I,J)=RELAX*FHI(I,J)+(1.0D+00-RELAX)*FHIOLD(I,J)
		end do
	end do

	ICUENTA=ICUENTA+1
	
	IF (ICUENTA.LT.NPAS) THEN            
                                  GOTO 321
                              ELSE
	END IF                                            

      RETURN  
END

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			SUBROUTINE LBL ADI(I0,J0,I1,J1,AP,AE,AW,AN,AS,B,FHI,NPAS,RELAX)		!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

!**************************************************************************
!**********         RESUELVE EL SIST. DE ECS. ORDINARIAS         **********        
!*****      (METHOD: LINE BY LINE ALTERNATING DIRECTION IMPLICIT)    ******        
!**************************************************************************
        
!***************   SOLUCIÓN DE LA MATRIZ TRIDIAGONAL   ********************
!**********  LA MATRIZ ES RESUELTA USANDO EL ALGORITMO DE THOMAS   ********

SUBROUTINE LBL_ADI(I0,J0,I1,J1,AP,AE,AW,AN,AS,B,FHI,NPAS,RELAX)                 
!C        

implicit none

integer :: i2,i1,i0,j2,j1,j0,iCUENTA,i,j,npas 
real(KIND = DP)  :: beta,denom,relax
real(KIND = DP), dimension(Nx,Ny) :: FHI,AE,AW,AN,AS,AP,B,FHIOLD,FHIOLDX
real(KIND = DP), dimension(max(Nx,Ny)) :: PT,QT

!CC    I0 = 1er. NODO COMPUTACIONAL EN X
!CC    J0 = 1er. NODO COMPUTACIONAL EN J
!CC    I1 = ULTIMO NODO COMPUTACIONAL EN X
!CC    J1 = ULTIMO NODO COMPUTACIONAL EN Y
!CC    FHI = VARIABLE A DETERMINAR
!CC    PT = RELACIÓN DE RECURRENCIA
!CC    QT = RELACIÓN DE RECURRENCIA

I2=I1-1
J2=J1-1

i=0
j=0
beta=0.0d+0
denom=0.0d+0
FHIOLD(:,:) = 0.0d+0	
FHIOLDX(:,:) = 0.0d+0
PT(:) = 0.0d+0
QT(:) =0.0d+0

Do ICUENTA=1,NPAS

!**************************************************************************
!**************     RENOMBRAMIENTO DE LA VARIABLE FHI      ****************
!**************************************************************************
	DO I=I0,I1              
		DO J=J0,J1 
			FHIOLD(I,J)=FHI(I,J)
		end do
	end do
!**************************************************************************
!************   BARRIDO DE PUNTOS NODALES EN DIRECCIÓN "X"   **************
!**************************************************************************
	DO J=J0,J1             
!CC    ****  1er.Punto en "X"  ****
!C
		IF(J.EQ.J0)   THEN
			BETA=B(I0,J)+AN(I0,J)*FHIOLD(I0,J+1) 
  
		ELSE IF(J.EQ.J1)  THEN
			BETA=B(I0,J)+AS(I0,J)*FHIOLD(I0,J-1) 

		ELSE               
			BETA=B(I0,J)+AS(I0,J)*FHIOLD(I0,J-1)+ AN(I0,J)*FHIOLD(I0,J+1)
		END IF
      
		DENOM=AP(I0,J)
		PT(I0)=AE(I0,J)/DENOM
		QT(I0)=BETA/DENOM
!C      
!CC    ****  Puntos Centales en "X"  ****
!C
		DO I=I0+1,I2
			IF(J.EQ.J0) THEN
				BETA=B(I,J)+AN(I,J)*FHIOLD(I,J+1)    
			ELSE IF(J.EQ.J1) THEN
				BETA=B(I,J)+AS(I,J)*FHIOLD(I,J-1)  
			ELSE
				BETA=B(I,J)+AS(I,J)*FHIOLD(I,J-1)+AN(I,J)*FHIOLD(I,J+1)
			END IF

			DENOM=AP(I,J)-PT(I-1)*AW(I,J)
			PT(I)=AE(I,J)/DENOM
			QT(I)=(BETA+QT(I-1)*AW(I,J))/DENOM
		end do
!C      
!CC    ****  Ultimo Punto en "X"  ****
!C
		IF(J.EQ.J0) THEN
			BETA=B(I1,J)+AN(I1,J)*FHIOLD(I1,J+1)   
		ELSE IF(J.EQ.J1) THEN
			BETA=B(I1,J)+AS(I1,J)*FHIOLD(I1,J-1)  
		ELSE
			BETA=B(I1,J)+AS(I1,J)*FHIOLD(I1,J-1) + AN(I1,J)*FHIOLD(I1,J+1)
		END IF

		DENOM=AP(I1,J)-PT(I2)*AW(I1,J)
		QT(I1)=(BETA+QT(I2)*AW(I1,J))/DENOM
!C
!CC    Solución de la Variable FHI(i,j)
!C      
		FHI(I1,J)=QT(I1)
!C      
		DO I=I2,I0,-1
			FHI(I,J)=PT(I)*FHI(I+1,J)+QT(I)
		end do
	end do
 
!C
!CC    Finaliza el Barrido en dirección "X"
!C                                                              
!**************************************************************************
!**************     RENOMBRAMIENTO DE LA VARIABLE FHI      ****************
!***********    SE UTILIZA LA SOLUCIÓN DEL BARRIDO EN "X"     *************
!**************************************************************************
!C
	DO I=I0,I1              
		DO J=J0,J1 
			FHIOLDX(I,J)=FHI(I,J)
		end do
	end do
!C
!**************************************************************************
!************   BARRIDO DE PUNTOS NODALES EN DIRECCIÓN "Y"   **************
!**************************************************************************
!
	DO I=I0,I1             
!C      
!CC    ****  1er.Punto en "Y"  ****
!C
		IF(I.EQ.I0) THEN
			BETA=B(I,J0)+AE(I,J0)*FHIOLDX(I+1,J0)   
		ELSE IF(I.EQ.I1) THEN
			BETA=B(I,J0)+AW(I,J0)*FHIOLDX(I-1,J0)  
		ELSE  
			BETA=B(I,J0)+AW(I,J0)*FHIOLDX(I-1,J0)+AE(I,J0)*FHIOLDX(I+1,J0)
		END IF

		DENOM=AP(I,J0)
		PT(J0)=AN(I,J0)/DENOM
		QT(J0)=BETA/DENOM
!C      
!CC    ****  Puntos Centales en "Y"  ****
!C
		DO J=J0+1,J2	
			IF(I.EQ.I0) THEN
				BETA=B(I,J)+AE(I,J)*FHIOLDX(I+1,J)    
			ELSE IF(I.EQ.I1) THEN
				BETA=B(I,J)+AW(I,J)*FHIOLDX(I-1,J)   
			ELSE
				BETA=B(I,J)+AW(I,J)*FHIOLDX(I-1,J)+AE(I,J)*FHIOLDX(I+1,J)  
			END IF
			DENOM=AP(I,J)-PT(J-1)*AS(I,J)
			PT(J)=AN(I,J)/DENOM
			QT(J)=(BETA+QT(J-1)*AS(I,J))/DENOM  
		end do
!C      
!CC    ****  Ultimo Punto en "Y"  ****
!C
		IF(I.EQ.I0)   THEN
			BETA=B(I,J1)+AE(I,J1)*FHIOLDX(I+1,J1)   
		ELSE IF(I.EQ.I1) THEN
			BETA=B(I,J1)+AW(I,J1)*FHIOLDX(I-1,J1)    
		ELSE
			BETA=B(I,J1)+AW(I,J1)*FHIOLDX(I-1,J1)+AE(I,J1)*FHIOLDX(I+1,J1)  
		END IF      

		DENOM=AP(I,J1)-PT(J2)*AS(I,J1)
		QT(J1)=(BETA+QT(J2)*AS(I,J1))/DENOM
!C
!CC    Solución de la Variable FHI(i,j)
!C      
		FHI(I,J1)=QT(J1)
!C      
		DO J=J2,J0,-1
			FHI(I,J)=PT(J)*FHI(I,J+1)+QT(J)
		end do
	end do
!C
!CC    Finaliza el Barrido en dirección "Y"
!C
!**************************************************************************
!*************     BAJA-RELAJACIÒN DE LA VARIABLE FHI      ****************
!**************************************************************************
!C
	DO I=I0,I1              
		DO J=J0,J1 
			FHI(I,J)=RELAX*FHI(I,J)+(1.0D+00-RELAX)*FHIOLD(I,J)
		end do
	end do

end do                                         

RETURN  
END

!**************************************************************************
!**********                   SUBRUTINA MASSFLUXES               **********         
!***********************        (FLUJOS MASICOS)      *********************        
!**************************************************************************
!C        
      SUBROUTINE MASSFLUXES(RHO,U,V,FEP,FNP,k)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	integer :: k
!C  VARIABLES DE FLUJO.
      DIMENSION U(NX,Ny),V(NX,Ny)
!C  FLUJOS MASICOS
      DIMENSION FEP(NX,Ny),FNP(NX,Ny)
!C
!*+++++++++++++++

!*+++++++++++++++
!C
        DO  I=1,NX        
        DO  J=1,NY
                      FEP(I,J)=0.0D+00                                               
                      FNP(I,J)=0.0D+00                                               
 		end do
 		end do
!C       
!CC    CALCULO DE LOS FLUJOS MASICOS   
!C
      Do J=2,NY-1 
      Do I=2,NX-2
!C
!CC    CALCULO DEL FLUJO MASICO "FEP" EN LA DERECHA DE LA MALLA PRINCIPAL "P,T"
!CC    (DOMINIO INTERNO (NODOS CENTRALES) DE LA MALLA DESPLAZADA EN "X")    
!CC    EN EL DOMINIO INTERNO (NODOS CENTRALES)   
!C
                      FEP(I,J)=RHO*U(I,J)*DYP(J)
 		end do
 		end do
!C
!CC    CALCULO DEL FLUJO MASICO "FNP" ARRIBA DE LA MALLA PRINCIPAL "P,T"
!CC    (DOMINIO INTERNO (NODOS CENTRALES) DE LA MALLA DESPLAZADA EN "Y")    
!CC    EN EL DOMINIO INTERNO (NODOS CENTRALES)   
!C
      Do J=2,NY-2 
      Do I=2,NX-1
                      FNP(I,J)=RHO*V(I,J)*DXP(I)
 	  end do
 	  end do
!C
!CC    CALCULO DEL FLUJO MASICO "FEP" EN LA DERECHA DE LA MALLA PRINCIPAL "P,T"
!CC    (FRONTERAS (NODOS FRONTERAS) DE LA MALLA DESPLAZADA EN "X")    
!CC    EN LAS FRONTERAS   
!C
!CC     LADO OESTE
!C
!***      I=1
      Do J=2,NY-1 
                      FEP(1,J)=0.0D+00                                        ! "FRONTERA SOLIDA"
!**                      FEP(1,J)=RHO*U(1,J)*DYP(J)                               "FRONTERA FLUIDA"
 	  end do
!C
!CC     LADO ESTE
!C
!***      I=NX-1
      Do J=2,NY-1 
                      FEP(NX-1,J)=0.0D+00                                   ! "FRONTERA SOLIDA"
!**                      FEP(NX-1,J)=RHO*U(NX-1,J)*DYP(J)                         "FRONTERA FLUIDA"
 	  end do

!C
!CC    CALCULO DEL FLUJO MASICO "FNP" ARRIBA DE LA MALLA PRINCIPAL "P,T"
!CC    (FRONTERAS (NODOS FRONTERAS) DE LA MALLA DESPLAZADA EN "Y")    
!CC    EN LAS FRONTERAS   
!C
!CC     LADO SUR
!C
!***      J=1
      Do I=2,NX-1
                      FNP(I,1)=RHO*V(I,1)*DXP(I) 
!****                      FNP(I,1)=0.0D+00                                         "FRONTERA SOLIDA"
             IF (K.NE.1)  FNP(I,1)=0.0D+00  
!**                      FNP(I,1)=RHO*V(I,1)*DXP(I)                               "FRONTERA FLUIDA"
 	  end do

!C
!CC     LADO NORTE
!C
!***      J=NY-1
      Do I=2,NX-1
                      FNP(I,NY-1)=RHO*V(I,NY-1)*DXP(I) 
!****                      FNP(I,NY-1)=0.0D+00                                      "FRONTERA SOLIDA"
             IF (K.NE.1)  FNP(I,NY-1)=0.0D+00  
!**                      FNP(I,NY-1)=RHO*V(I,NY-1)*DXP(I)                         "FRONTERA FLUIDA"
 		end do
   

      RETURN  
      END


!**************************************************************************
!*****   FUNCIÒN PROPIEDAD (RESUELVE LA PROPIEDAD EN LA INTERFACE)   ******        
!**************************************************************************
!C        
      FUNCTION PROPIEDAD(X1,Y1,F1)
!C        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!C
!CC    CALCULO DE LA PROPIEDAD EFECTIVA EN LA INTERFACE
!C
      PROPIEDAD=1.0D+00/(((1.0D+00-F1)/X1)+(F1/Y1))
!C
      RETURN            
      END 





end module cfdlib
