module rans_models
	use user
	implicit none

	
contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						f_mu						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function MUT_model(i,j,Tk,e,rho,mu)

!			MuT(i,j)=Cmu*RHO(i,j)*TK(i,j)*TK(i,j) / (e(i,j)+1D-25)
	implicit none
	integer :: i,j
	real(KIND = DP)  ::  Re_t, f_mu, arg, arg2, MUT_model

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		mu,rho,e,TK


	SELECT CASE (rans_model)
		CASE (0) 
			MUT_model = 0.0D+0

		CASE (1) ! k-epsilon HH

			Re_t = ABS( RHO(i,j)*TK(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			f_mu = 1.D+0

			MUT_model = Cmu*abs(f_mu)*RHO(i,j)*TK(i,j)*TK(i,j) / (e(i,j)+1D-25)

		CASE (2) ! k-epsilon JL

			Re_t = ABS( RHO(i,j)*TK(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			arg = 2.5D+0 / (1.0D+0 + Re_t/50.0D+0 )  

			if (arg .lt. 25.0D+0 ) then
				f_mu = exp(-arg)
			else
				f_mu=0.0D+0
				!write(*,*) " warning : limiting turbulent vicosity "
			end if

			MUT_model = Cmu*abs(f_mu)*RHO(i,j)*TK(i,j)*TK(i,j) / (e(i,j)+1D-25)		
			

		CASE (3) ! k-epsilon LS

			Re_t = ABS( RHO(i,j)*TK(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			arg = 3.4D+0 / ( (1.0D+0 + Re_t/50.0D+0 )*(1.0D+0 + Re_t/50.0D+0 ) )

			if (arg .lt. 25.0D+0 ) then
				f_mu = exp(-arg)
			else
				f_mu=0.0D+0
				!write(*,*) " warning : limiting turbulent vicosity "
			end if


			MUT_model = Cmu*abs(f_mu)*RHO(i,j)*TK(i,j)*TK(i,j) / (e(i,j)+1D-25)

		CASE (4) ! k-w PHD

			Re_t = ABS( RHO(i,j)* TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )


                     ARG=DABS((RE_T/10.0D+00)**0.75D+00)
                     ARG2=DABS((RE_T/200.0D+00)**2.0D+00)
     			F_MU=0.025D+00 + ( ( 1.0D+00-DEXP(-ARG)) * &
	    		( 0.975D+00+((0.001/(RE_T+1.0D-30))*DEXP(-ARG2))) )  !"MODELO PENG-DAVIDSON-HOLMBERG (PDH)"

			MUT_model = Cmu*abs(f_mu)*RHO(i,j)*TK(i,j) / (e(i,j)+1D-25)




		CASE (5) ! k-w WX 1994

			Re_t = ABS( RHO(i,j)* TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			f_mu = ( 0.025D+00 + Re_t/6.0D+00 ) / ( 1.0D+00 + Re_t/6.0D+00 )

			MUT_model = Cmu*abs(f_mu)*RHO(i,j)*TK(i,j) / (e(i,j)+1D-25)
			
		CASE DEFAULT
			stop "Selected RANS model unknown"
	END SELECT

	return

end function MUT_model

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				init_rans_constants					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine init_rans_constants()

	implicit none

	SELECT CASE (rans_model)
		CASE (0)
			Cmu=1.0D+0
			C1e=1.0D+0
			C2e=1.0D+0
			Theta_k=1.0D+0
			Theta_e=1.0D+0
			Theta_t=1.0D+0

		CASE (1) ! k-epsilon HH

			Cmu=0.09D+0
			C1e=1.44D+0
			C2e=1.92D+0
			Theta_k=1.0D+0
			Theta_e=1.3D+0
			Theta_t=0.9D+0

		CASE (2) ! k-epsilon JL

			Cmu=0.09D+0
			C1e=1.55D+0
			C2e=2.0D+0
			C3e = 0.0D+0
			Theta_k=1.0D+0
			Theta_e=1.3D+0
			Theta_t=0.9D+0	
			

		CASE (3) ! k-epsilon LS

		!	Cmu=0.09D+0
			Cmu=0.115D+0		!aumento valor per veure si el jet s'obre més	
			C1e=1.44D+0
			C2e=1.92D+0
			C3e = 0.0D+0

			Theta_k=1.0D+0

			Theta_e=1.3D+0
!			Theta_t=0.9D+0
			Theta_t=0.5D+0
		CASE (4) ! k-w PHD

			Cmu = 1.0D+0
			c_k = 0.09D+0
			cw1 = 0.42D+0
			cw2 = 0.075D+0
			crw = 0.75D+0

			Theta_k=0.8D+0
			Theta_e=1.35D+0
			Theta_t=0.9D+0

		CASE (5) ! k-w WX 1994

			Cmu = 1.0D+0
			c_k = 0.09D+0
			cw1 = 0.56D+0
			cw2 = 0.075D+0

			Theta_k=2.0D+0
			Theta_e=2.3D+0
			Theta_t=0.9D+0
			
		CASE DEFAULT
			stop "Selected RANS model unknown"
	END SELECT

	return

end subroutine init_rans_constants



!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			 		Rans_source_terms_TK					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine rans_source_terms_TK(i,j,u,v,Mut,mu,T,rho,e,TK,sc,sp)

	implicit none
		   
	integer:: i,j

	real(KIND = DP)  :: 				&
		dtdy,Txyt_P, Un, Up, Us, Vw, Vp, Ve, Re_t, Re_te3, Re_te4,	&
		dudx, dudy, dvdx, dvdy, sp, sc, dk12_dx ,dk12_dy, D_source, E_source, f_g , f_k, &
		factor1, factor2

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		U,V,Tk,mu,MuT,e,rho,		&
		Fep,Fnp,T


		Un = ( U(i-1,j+1) + U(i,j+1) ) / 2.0D+0
		Up = ( U(i-1,j)   + U(i,j) )   / 2.0D+0
		Us = ( U(i-1,j-1) + U(i,j-1) ) / 2.0D+0

		Vw = ( V(i-1,j) + V(i-1,j-1) ) / 2.0D+0
		Vp = ( V(i,j)   + V(i,j-1) ) / 2.0D+0
		Ve = ( V(i+1,j) + V(i+1,j-1) ) / 2.0D+0

		dudx = ( U(i,j) - U(i-1,j) ) / DXP(i)
		dudy = Deri(Us, Up, Un, dyv(j-1), dyv(j))
		dvdx = Deri(Vw, Vp, Ve, dxu(i-1), dxu(i))
		dvdy = (V(i,j) - V(i,j-1)) / DYP(j)
		dtdy = Deri( T(i,j-1),T(i,j),T(i,j+1),DYV(j-1),DYV(j)) 

		Txyt_P = muT(i,j) * (dudy + dvdx)

		Pk(i,j) = Txxt(i,j)*dudx + Txyt_P*(dudy + dvdx) + Tyyt(i,j)*dvdy
		Gk(i,j) = (-Mut(i,j) / theta_t) * g * bbeta * dtdy



	SELECT CASE (rans_model)
		CASE (0) 
			Sc = 0.0D+0; Sp = 0.0D+0
		CASE (1) ! k-epsilon HH

			Sc =  ( Pk(i,j) + DMAX1(Gk(i,j),0.0D+0) )
			Sp = -( rho(i,j)*e(i,j) - DMIN1(Gk(i,j),0.0D+0) ) / (Tk(i,j)+1.0D-25)

		CASE (2:3) ! k-epsilon JL and LS
	
			dk12_dx = Deri( sqrt(abs(TK(i-1,j))) , sqrt(abs(TK(i,j))), sqrt(abs(TK(i+1,j))), dxu(i-1), dxu(i))
			dk12_dy = Deri( sqrt(abs(TK(i,j-1))) , sqrt(abs(TK(i,j))), sqrt(abs(TK(i,j+1))),DYV(j-1),DYV(j)) 

			D_source = 2.0D+0 * mu(i,j) * (dk12_dx*dk12_dx + dk12_dy*dk12_dy)

			Sc =  ( Pk(i,j) + DMAX1(Gk(i,j),0.0D+0) - DMIN1(D_source,0.0D+0) )
			Sp = -( rho(i,j)*e(i,j) + DMAX1(D_source,0.0D+0) - DMIN1(Gk(i,j),0.0D+0) ) / (Tk(i,j)+1.0D-25)		
		

		CASE (4) ! k-w PHD

			Re_t = abs( RHO(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )
			
			Re_te3 = (Re_t/12.0D+0)*(Re_t/12.0D+0)*(Re_t/12.0D+0)
			Re_te4 = (Re_t/10.0D+0)*(Re_t/10.0D+0)*(Re_t/10.0D+0)*(Re_t/10.0D+0)

			factor1 =  1.0D+0 - exp(-(Re_te3))
			factor2 = 1.0D+00 + (10.0D+00 / ( (Re_t**3.25D+00)+1.0D-25) ) 
			f_g = factor1 * factor2
			f_k =  1.0D+0 - 0.722D+0 * exp(- Re_te4)		

			D_source = C_k * f_k * rho(i,j) * TK(i,j) * e(i,j)

			Gk(i,j) = f_g * Gk(i,j)

			Sc =  ( Pk(i,j) + DMAX1(Gk(i,j),0.0D+0) - DMIN1(D_source,0.0D+0) )
			Sp = -( DMAX1(D_source,0.0D+0) - DMIN1(Gk(i,j),0.0D+0) ) / (Tk(i,j)+1.0D-25)		

		CASE (5) ! k-w WX 1994

			Re_t = abs( RHO(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )
			
			Re_te4 = (Re_t/8.0D+0)*(Re_t/8.0D+0)*(Re_t/8.0D+0)*(Re_t/8.0D+0)

			f_g = 0.0D+0 
			f_k = ( 5.0D+0/18.0D+0 + Re_te4 ) / ( 1.0D+0 + Re_te4 )		

			D_source = C_k * f_k * rho(i,j) * TK(i,j) * e(i,j)

			Gk(i,j) = f_g * Gk(i,j)

			Sc =  ( Pk(i,j) + DMAX1(Gk(i,j),0.0D+0) - DMIN1(D_source,0.0D+0) )
			Sp = -( DMAX1(D_source,0.0D+0) - DMIN1(Gk(i,j),0.0D+0) ) / (Tk(i,j)+1.0D-25)
			
		CASE DEFAULT
			stop "Selected RANS model unknown"
	END SELECT

end subroutine rans_source_terms_tk


!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			 		Rans_source_terms_e					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine rans_source_terms_e(i,j,u,v,Mut,mu,T,rho,e,TK,sc,sp)

	implicit none
		   
	integer:: i,j

	real(KIND = DP)  :: 								&
		sc,sp,Dis,up,vp,E_source,f1,f2,A_source,D_source,f_mu,Re_t,Re_te2,	&
		d2u_dy2,d2v_dx2, dkdx,dkdy,dwdx,dwdy,Un,Us,Vw,Ve

	real(KIND = DP), dimension(Nx,Ny) :: 	&
		U,V,Mut,mu,T,rho,e,TK



	SELECT CASE (rans_model)
		CASE (0) 
			Sc = 0.0D+0; Sp = 0.0D+0
		CASE (1) ! k-epsilon HH

			Up = ( U(i-1,j) + U(i,j) ) / 2.0D+0
			Vp = ( V(i,j) + V(i,j-1) ) / 2.0D+0

			C3e = Dtanh( dabs(Vp/(Up+1.0D-25)) )	

			Pk(i,j) = C1e*Pk(i,j)*e(i,j)     / TK(i,j)
			Gk(i,j) = C1e*C3e*Gk(i,j)*e(i,j) / Tk(i,j)
			Dis = C2e*RHO(i,j)*e(i,j)*e(i,j) / TK(i,j)	

			Sc = ( Pk(i,j) + DMAX1(Gk(i,j),0.0D+0) )
			Sp = -( Dis - DMIN1(Gk(i,j),0.0D+0) ) / (e(i,j)+1.0D-25)

		CASE (2:3) ! k-epsilon JL and LS

			Re_t = ABS( RHO(i,j)*TK(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			f1 = 1.0D+0

			Re_te2 = Re_t * Re_t

			if (Re_te2 .Lt. 25.0D+0) then
				f2 = 1.0D+0 - 0.3D+0 * exp( - Re_te2)
			else 	
				f2=1.0D+0
			end if
	
			Pk(i,j) = C1e*abs(f1)*Pk(i,j)*e(i,j)     / (TK(i,j)+1.0D-25)
			Gk(i,j) = C1e*C3e*Gk(i,j) *e(i,j)        / (Tk(i,j)+1.0D-25)
			Dis = C2e*abs(f2)*RHO(i,j)*e(i,j)*e(i,j) / (TK(i,j)+1.0D-25)


			Un = ( U(i-1,j+1) + U(i,j+1) ) / 2.0D+0
			Up = ( U(i-1,j)   + U(i,j) )   / 2.0D+0
			Us = ( U(i-1,j-1) + U(i,j-1) ) / 2.0D+0

			Vw = ( V(i-1,j) + V(i-1,j-1) ) / 2.0D+0
			Vp = ( V(i,j)   + V(i,j-1) ) / 2.0D+0
			Ve = ( V(i+1,j) + V(i+1,j-1) ) / 2.0D+0


			d2u_dy2 = Deri2(Us, Up, Un, dyv(j-1), dyv(j))
			d2v_dx2 = Deri2(Vw, Vp, Ve, dxu(i-1), dxu(i))

			E_source = 2.0D+0* mu(i,j)* MuT(i,j)* (d2u_dy2*d2u_dy2 + d2v_dx2*d2v_dx2) / RHO(i,j)	

			Sc = ( Pk(i,j) + DMAX1(Gk(i,j),0.0D+0)  + DMAX1(E_source,0.0D+0) )
			Sp = -( Dis - DMIN1(Gk(i,j),0.0D+0)  - DMIN1(E_source,0.0D+0) ) / (e(i,j)+1.0D-25)	
		

			!if ( (j.eq.120) .and. ((i.ge.1).and.(i.le.10)) ) then

			!	write(*,*) 'i,j,Pk(i,j)'
			!	write(*,*)  i,',',j,',',Pk(i,j)
			!	write(*,*)   
			!	write(*,*) 'i,j,TK(i,j)'
			!	write(*,*)  i,',',j,',',TK(i,j)
			!	write(*,*)                        
			!end if 


		CASE (4) ! k-w PHD

			Re_t = abs( RHO(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			f1 = 1.0D+0 + 4.3D+0 * exp (- sqrt(Re_t/1.5D+0))
			f2 = 1.0D+0

			dkdx = Deri(TK(i-1,j), TK(i,j), TK(i+1,j), dxu(i-1), dxu(i)) 
			dkdy = Deri(TK(i,j-1), TK(i,j), TK(i,j+1) , dyv(j-1), dyv(j))
			dwdx = Deri(e(i-1,j), e(i,j), e(i+1,j), dxu(i-1), dxu(i))
			dwdy = Deri(e(i,j-1), e(i,j), e(i,j+1) , dyv(j-1), dyv(j)) 
			
			Pk(i,j)  = Cw1 * f1 * Pk(i,j) * e(i,j) / TK(i,j)
			A_source = Cw2 * f2 * rho(i,j)* e(i,j) * e(i,j) 				
			D_source = Crw * muT(i,j) * (dkdx*dwdx + dkdy*dwdy) / TK(i,j)

			Sc = ( DMAX1(Pk(i,j),0.0D+0) + DMAX1(D_source,0.0D+0)  )
			Sp = -( -DMIN1(Pk(i,j),0.0D+0) + A_source  - DMIN1(D_source,0.0D+0) ) / (e(i,j)+1.0D-25)

		CASE (5) ! k-w WX 1994

			Re_t = abs( RHO(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )

			f_mu = ( 0.025D+00 + Re_t/6.0D+00 ) / ( 1.0D+00 + Re_t/6.0D+00 )

			f1 = ( (0.1D+0 + (Re_t/2.7D+0))  /  (1.0D+0 + (Re_t/2.7D+0))  ) / (abs(f_mu) + 1.0D-30)
			f2 = 1.0D+0
			
			Pk(i,j)  = Cw1 * abs(f1) * Pk(i,j) * e(i,j) / TK(i,j)
			A_source = Cw2 * abs(f2) * rho(i,j)* e(i,j) * e(i,j) 				

			Sc = ( DMAX1(Pk(i,j),0.0D+0) )
			Sp = -( -DMIN1(Pk(i,j),0.0D+0) + A_source ) / (e(i,j)+1.0D-25)
			
		CASE DEFAULT
			stop "Selected RANS model unknown"
	END SELECT

end subroutine rans_source_terms_e

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					rans_boundary						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function rans_boundary_e(i,j)


	implicit none
	integer :: i,j
	real(KIND = DP)  :: rans_boundary_e

	SELECT CASE (rans_model)
		CASE (0) 
			rans_boundary_e = 0.0D+0
		CASE (1) ! k-epsilon HH

			rans_boundary_e = 100.0D+0

		CASE (2:3) ! k-epsilon JL and k-epsilon LS

			rans_boundary_e = 0.0D+0
			

		CASE (4:5) ! k-w PHD and k-w WX 1994

			if (i==1) rans_boundary_e =  6.0D+0*Mu(i,j) / (Cw2 * RHO(i,j) * DXU(i)* DXU(i))
			if (i==Nx) rans_boundary_e = 6.0D+0*Mu(i,j) / (Cw2 * RHO(i,j) * DXU(Nx-1)* DXU(Nx-1))

			if (j==1) rans_boundary_e =  6.0D+0*Mu(i,j) / (Cw2 * RHO(i,j) * DYV(j)* DYV(j))
			if (j==Ny) rans_boundary_e = 6.0D+0*Mu(i,j) / (Cw2 * RHO(i,j) * DYV(Ny-1)* DYV(Ny-1))

		CASE DEFAULT
			stop "Selected RANS model unknown"
	END SELECT

	return

end function rans_boundary_e
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				First order derivative						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function deri(y1,y0,y2,x01,x20) 

	implicit none
	real(KIND = DP):: y1,y0,y2,x01,x20,deri
       

	deri = ( (y2-y0)*(x01)*(x01) - (y1-y0)*(x20)*(x20) ) / (x01*x20*(x20+x01))
	
	return     
end 

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Second order derivative						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

function deri2(y1,y0,y2,x01,x20) 

	implicit none
	real(KIND = DP):: y1,y0,y2,x01,x20,deri2
       

	deri2 = 2.0D+0*( y1*x20 - y0*(x20+x01)+(y2*x01) ) / (x01*x20*(x20+x01))
	
	return     
end 

end module rans_models
