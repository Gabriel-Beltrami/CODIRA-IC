module visual
	use user
	implicit none
	
contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!	Calculo del Nusselt,temperaturas y velocidades adimensionales				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine literatura_compare()

	implicit none
	real(KIND = DP), dimension(Ny) :: aux

	if (n_high_order_scheme.EQ.0) &
			open (2,FILE='./outputdata/diag.dat',STATUS='unknown',access='sequential')
	if (n_high_order_scheme.EQ.1) &
			open (2,FILE='./outputdata/diag_1.dat',STATUS='unknown',access='sequential')
	if (n_high_order_scheme.EQ.2) &
			open (2,FILE='./outputdata/diag_2.dat',STATUS='unknown',access='sequential')
	if (n_high_order_scheme.EQ.3) &
			open (2,FILE='./outputdata/diag_3.dat',STATUS='unknown',access='sequential')
	if (n_high_order_scheme.EQ.4) &
			open (2,FILE='./outputdata/diag_4.dat',STATUS='unknown',access='sequential')
	if (n_high_order_scheme.EQ.5) &
			open (2,FILE='./outputdata/diag_5.dat',STATUS='unknown',access='sequential')
	if (n_high_order_scheme.EQ.6) &
			open (2,FILE='./outputdata/diag_6.dat',STATUS='unknown',access='sequential')


	do i = 1, Nx
 		aux(i) = T(i,Ny-i+1)
	end do

	do i=1,Nx
		write(2,*) aux(i)
	end do

	close(2)
	

end subroutine literatura_compare

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing 					!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine write_vtk(U,V,P,T,C,b,iter_t,TK,e,MuT)

	implicit none

	character(len = 27) :: str, tmp
	integer :: i,j,iter_t
	real(KIND = DP), dimension(Nx+1) :: x_print
	real(KIND = DP), dimension(Ny+1) :: Y_print
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,C,b,PSIP,TK,e,inter_U,inter_V,MuT,	&
										 Pk,epsi,div_T,DTK_Dt,vorticity_z

	x_print(1)=XU(1)
	do i = 1,Nx-1
		x_print(i+1)=XU(i)
	end do
	x_print(Nx+1)=XU(Nx-1)+1.0D-5 	

	y_print(1)=yv(1)
	do j = 1,Ny-1
		y_print(j+1)=yv(j)
	end do
	y_print(Ny+1)=yv(Ny-1)+1.0D-5
	x_print(2)=XU(1)+1.0D-5
	y_print(2)=yv(1)+1.0D-5


	write (tmp,'(I5.5)') case_jet
	str = 'outputdata/output_'//trim(tmp)//'.vtk'
	open(unit=1, file = str)

	write ( 1,'(a)' )'# vtk DataFile Version 3.0'
	write ( 1,'(a)' ) 'Non-uniform Rectilinear - Rectilinear Grid'
	write ( 1,'(a)' ) 'ASCII'
	write ( 1,'(a)' ) 'DATASET RECTILINEAR_GRID'
	write ( 1,'(a,3i4)' )'DIMENSIONS', nx+1,ny+1,1
	write ( 1,'(a13,1i4,a6)' ) 'X_COORDINATES',nx+1,' float'
	write ( 1,'(f22.14)' ) (x_print(i),i=1,nx+1)
	write ( 1,'(a13,1i4,a6)' ) 'Y_COORDINATES',ny+1,' float'
	write ( 1,'(f22.14)' ) (y_print(j),j=1,ny+1)
	write ( 1,'(a13,1i4,a6)' ) 'Z_COORDINATES',1,' float'
	write ( 1,'(f22.14)' ) 1.0

	write ( 1, '(a, I16)' ) 'CELL_DATA ', (nx)*(ny)*1

	write ( 1, '(a)' ) 'SCALARS P float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((P(i,j),i=1,(nx)),j=1,(ny))
	
	write ( 1, '(a)' ) 'SCALARS T float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((T(i,j),i=1,(nx)),j=1,(ny))
	
	write ( 1, '(a)' ) 'SCALARS C float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((C(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS e float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((e(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS Tk float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((TK(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS MuT float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((MuT(i,j),i=1,(nx)),j=1,(ny))

    CALL STREAMLINES(DYP,U,PSIP,RHO)

	write ( 1, '(a)' ) 'SCALARS Str float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((PSIP(i,j),i=1,(nx)),j=1,(ny))

	CALL energy_budget(Pk,epsi,div_T,DTK_Dt,vorticity_z)

	write ( 1, '(a)' ) 'SCALARS Production float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((Pk(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS -Dissipation float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((epsi(i,j),i=1,(nx)),j=1,(ny))	

	write ( 1, '(a)' ) 'SCALARS Tubulent_transport float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((div_T(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS -Advection float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((DTK_Dt(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS vorticity_z float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((vorticity_z(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS Rij_T float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((Rij_T(i,j),i=1,(nx)),j=1,(ny))
	
	write ( 1, '(a)' ) 'SCALARS Rij_C float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((Rij_C(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS Rij_u float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((Rij_u(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS Rij_V float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((Rij_v(i,j),i=1,(nx)),j=1,(ny))

! 	write ( 1, '(a)' ) 'SCALARS Sdc_V float'
! 	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
! 	write ( 1, '(f22.14)' ) ((Sdc_ij_v(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS RHO(Po,T) float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((rho(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS P_rad float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.8)' ) ((P_rad(i,j),i=1,(nx)),j=1,(ny))

	CALL interpol_UV(U,V,inter_U,inter_V)

	write ( 1, '(a)' ) 'SCALARS u float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.8)' ) ((inter_U(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'SCALARS v float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.8)' ) ((inter_V(i,j),i=1,(nx)),j=1,(ny))

	write ( 1, '(a)' ) 'VECTORS velocity float'
	do j=1,ny
		do i=1,nx
			write ( 1, "(f22.14,f22.14,f22.14)" ) inter_U(i,j), inter_V(i,j), 0.0D+00
		end do
	end do


	close(1)

!	TIMECODE=TIMEF()

end subroutine write_vtk

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					print_Residual						!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine print_Residual(iter,iter_t,Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Res_C,Rmax_T,Rmax_C,	&
					Res_Tk,Rmax_Tk,Res_e,Rmax_e) 

	real(KIND = DP) :: Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Res_C,Rmax_T,Rmax_C,Res_Tk,Rmax_Tk,Res_e,Rmax_e
	integer :: iter,iter_t

!	write(*,*) 'Iteracao no tempo:',iter_t
	write(*,*) "______________________________________________"
	write(*,*) 'Iteracao no SIMPLE:',iter	
	write(*,*)
	write(*,*) 'Global error of U:',Res_U,"  |  ",'MAX local error of U:',Rmax_U
	write(*,*) "-----------------------------------------------"	
	write(*,*) 'Global error of V:',Res_V,"  |  ",'MAX local error of V:',Rmax_V
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of P:',Res_P,"  |  ",'MAX local error of P:',Rmax_P
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of T:',Res_T,"  |  ",'MAX local error of T:',Rmax_T
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of C:',Res_C,"  |  ",'MAX local error of C:',Rmax_C
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of TK:',Res_TK,"  |  ",'MAX local error of TK:',Rmax_TK
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Global error of e:',Res_e,"  |  ",'MAX local error of e:',Rmax_e	
	write(*,*) "-----------------------------------------------"
	write(*,*) 'Maximum Turvulent viscosity:',maxval(maxval(Mut,2))
	return

end Subroutine print_Residual

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!!							Interpoling stagared values to the main mesh			!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine interpol_UV(U,V,inter_U,inter_V)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Nx,Ny) :: U,V,inter_U,inter_V


	do j=1,Ny
		do i=2,Nx-1
			inter_U(i,j) = (U(i,j)+U(i-1,j)) /2.0D+0	
		end do	 
	end do
	do j=1,Ny
		inter_U(1,j) = U(1,j)
		inter_U(Nx,j) = U(Nx-1,j)	
	end do

	do j=2,Ny-1
		do i=1,Nx
			inter_V(i,j) = (V(i,j)+V(i,j-1)) /2.0D+0	
		end do	 
	end do

	do i=1,nx
		inter_V(i,1) = V(i,1)
		inter_V(i,Ny) = V(i,Ny-1)	
	end do

end subroutine interpol_UV




subroutine interpol_txyt(inter_txyt)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Nx,Ny) :: inter_txyt


	do j=2,Ny-1
		do i=2,Nx-1

			inter_txyt(i,j) = (txyt(i,j) + txyt(i-1,j) + txyt(i,j-1) + txyt(i-1,j-1)) /4.0D+0	

		end do	 
	end do

! ! ! ! ! East and west Boundary! ! ! ! !! ! ! ! !! ! ! ! !
	do j=2,Ny-1

		inter_txyt(1,j) =  (txyt(1,j)  +  txyt(1,j-1))  / 2.0D+0
		inter_txyt(Nx,j) = (txyt(Nx-1,j) +  txyt(Nx-1,j-1)) / 2.0D+0

	end do
! ! ! ! !! ! ! ! !! ! ! ! !! ! ! ! !! ! ! ! !! ! ! ! !


! ! ! ! ! Nort and south Boundary ! ! ! ! !! ! ! ! !! ! ! ! !
	do i=2,Nx-1

		inter_txyt(i,1) =  (txyt(i,1)  +  txyt(i-1,1))  / 2.0D+0
		inter_txyt(i,Ny) = (txyt(i,Ny-1) +  txyt(i-1,Ny-1)) / 2.0D+0

	end do
! ! ! ! !! ! ! ! !! ! ! ! !! ! ! ! !! ! ! ! !! ! ! ! !	


	inter_txyt(1,1) = txyt(1,1)
	inter_txyt(1,Ny) = txyt(1,Ny-1)
	inter_txyt(Nx,1) = txyt(Nx-1,1)
	inter_txyt(Nx,Ny) = txyt(Nx-1,Ny-1)


end subroutine interpol_txyt
!-----------------------------------------------------------------------------------------------!




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
!!								  Energy_transport(ETT)												!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine Energy_transport(ETT)

	implicit none

	integer :: i,j

	real(KIND = DP) ::	dT_dy,GIT		

	real(KIND = DP), dimension(Nx,Ny) :: ETT,ETT_tempo

	do j=2,Ny-1
		do i=2,Nx-1

			GIT=MUT(i,j) / Theta_t

			dT_dy = Deri( T(i,j-1) , T(i,j), T(i,j+1), DYV(j-1), DYV(j)) 

!			ETT_tempo(i,j) = rho(i,j) * Cpo * GIT * dT_dy
			ETT_tempo(i,j) = - Cpo * GIT * dT_dy
		end do	 
	end do
	ETT_tempo(:,Ny)  = ETT_tempo(:,Ny-1) ; ETT_tempo(:,1)  = ETT_tempo(:,2)
	ETT_tempo(1,:)   = ETT_tempo(2,:)	 ; ETT_tempo(Nx,:) = ETT_tempo(Nx-1,:) 
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
	do j=2,Ny-1
		do i=2,Nx-1

			ETT(i,j) =	Deri( ETT_tempo(i,j-1) , ETT_tempo(i,j), ETT_tempo(i,j+1), DYV(j-1), DYV(j)) 
		end do	 
	end do
	ETT(:,Ny)  = ETT(:,Ny-1) ; ETT(:,1)  = ETT(:,2)
	ETT(1,:)   = ETT(2,:)	 ; ETT(Nx,:) = ETT(Nx-1,:) 	

end subroutine Energy_transport

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!!									energy_budget												!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine energy_budget(Pk,epsi,div_T,DTK_Dt,vorticity_z)

	implicit none

	integer :: i,j

	real(KIND = DP) ::	Txyt_P, Un, Up, Us, Vw, Vp, Ve, 							&
						dudx, dudy, dvdx, dvdy,										&
						dk12_dx, dk12_dy,D_source, f_k, Re_te4, Re_t,				&
						dk_dx_p,dk_dy_p,dk_dy_s,dk_dy_n,dk_dx_w,dk_dx_e, 			&
						Mu_total_w,Mu_total_p,Mu_total_e,Mu_total_s,Mu_total_n,		&
						uk_e, uk_p, uk_w, vk_s, vk_p, vk_n, duk_dx, dvk_dy			

	real(KIND = DP), dimension(Nx,Ny) :: Pk,epsi,div_T,DTK_Dt,inter_U,inter_V,vorticity_z





	do j=2,Ny-1
		do i=2,Nx-1

!-----------------------------------------------------------------------------------------------!
!!										Pk														!
!-----------------------------------------------------------------------------------------------!
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

		vorticity_z(i,j) = dvdx - dvdy

		Txyt_P = muT(i,j) * (dudy + dvdx)

		Pk(i,j) = Txxt(i,j)*dudx + Txyt_P*(dudy + dvdx) + Tyyt(i,j)*dvdy
		Pk(i,j) = Pk(i,j) / rho(i,j)		
!-----------------------------------------------------------------------------------------------!



!-----------------------------------------------------------------------------------------------!
!!										epsilon													!
!-----------------------------------------------------------------------------------------------!
	SELECT CASE (rans_model)
		CASE (0) 
			epsi(i,j) = 0.0D+0
		CASE (1) ! k-epsilon HH

			epsi(i,j) = e(i,j)

		CASE (2:3) ! k-epsilon JL and LS
	
			dk12_dx = Deri( sqrt(abs(TK(i-1,j))) , sqrt(abs(TK(i,j))), sqrt(abs(TK(i+1,j))), dxu(i-1), dxu(i))
			dk12_dy = Deri( sqrt(abs(TK(i,j-1))) , sqrt(abs(TK(i,j))), sqrt(abs(TK(i,j+1))),DYV(j-1),DYV(j)) 

			D_source = 2.0D+0 * mu(i,j) * (dk12_dx*dk12_dx + dk12_dy*dk12_dy)

			epsi(i,j) = e(i,j) + D_source/rho(i,j)	

		!	write(*,*) ' D/e = ', 	(D_source/rho(i,j)) / e(i,j)	

		CASE (4) ! k-w PHD

			Re_t = abs( RHO(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )			
			Re_te4 = (Re_t/10.0D+0)*(Re_t/10.0D+0)*(Re_t/10.0D+0)*(Re_t/10.0D+0)

			f_k =  1.0D+0 - 0.722D+0 * exp(- Re_te4)		

			epsi(i,j) = - C_k * f_k * TK(i,j) * e(i,j)
	

		CASE (5) ! k-w WX 1994

			Re_t = abs( RHO(i,j)*TK(i,j) / (Mu(i,j)*e(i,j) +1.0D-25) )	
			Re_te4 = (Re_t/8.0D+0)*(Re_t/8.0D+0)*(Re_t/8.0D+0)*(Re_t/8.0D+0)

			f_k = ( 5.0D+0/18.0D+0 + Re_te4 ) / ( 1.0D+0 + Re_te4 )		

			epsi(i,j) = - C_k * f_k * TK(i,j) * e(i,j)
			
		CASE DEFAULT
			stop "Selected RANS model unknown"
	END SELECT
!-----------------------------------------------------------------------------------------------!


!-----------------------------------------------------------------------------------------------!
!!										Div·T'													!
!-----------------------------------------------------------------------------------------------!


			dk_dx_p = Deri( TK(i-1,j) , TK(i,j), TK(i+1,j), dxu(i-1), dxu(i))
			dk_dy_p = Deri( TK(i,j-1) , TK(i,j), TK(i,j+1), DYV(j-1), DYV(j)) 

			if (j.eq.2)  then
					dk_dy_s = ( TK(i,2)-TK(i,1) ) / DYV(1)
				else
					dk_dy_s = Deri( TK(i,j-2) , TK(i,j-1), TK(i,j), DYV(j-2), DYV(j-1))
			end if			
			if (j.eq.Ny-1)  then
					dk_dy_n = ( TK(i,Ny)-TK(i,Ny-1) ) / DYV(Ny-1)
				else
					dk_dy_n = Deri( TK(i,j) , TK(i,j+1), TK(i,j+2), DYV(j), DYV(j+1))
			end if

			if (i.eq.2)  then
					dk_dx_w = ( TK(2,j)-TK(1,j) ) / dxu(1)
				else
					dk_dx_w = Deri( TK(i-2,j) , TK(i-1,j), TK(i,j), dxu(i-2), dxu(i-1))
			end if			
			if (i.eq.Nx-1)  then
					dk_dx_e = ( TK(Nx,j)-TK(Nx-1,j) ) / Dxu(Nx-1)
				else
					dk_dx_e = Deri( TK(i,j) , TK(i+1,j), TK(i+2,j), dxu(i), dxu(i+1))
			end if


			Mu_total_w = Mu(i-1,j) 	+ 	MuT(i-1,j)/theta_k 

			Mu_total_p = Mu(i  ,j) 	+ 	MuT(i,j)/theta_k 

			Mu_total_e = Mu(i+1,j) 	+ 	MuT(i+1,j)/theta_k 

			Mu_total_s = Mu(i,j-1) 	+ 	MuT(i,j-1)/theta_k 

			Mu_total_n = Mu(i,j+1) 	+ 	MuT(i,j+1)/theta_k

			div_t(i,j) = Deri( Mu_total_w *dk_dx_w , Mu_total_p*dk_dx_p, Mu_total_e*dk_dx_e, dxu(i-1), dxu(i)) &
				  + Deri( Mu_total_s *dk_dy_s , Mu_total_p*dk_dy_p, Mu_total_e*dk_dy_n, DYV(j-1), DYV(j))

			div_t(i,j) = div_t(i,j) / rho(i,j)		
!-----------------------------------------------------------------------------------------------!


!-----------------------------------------------------------------------------------------------!
!!										DTK_Dt													!
!-----------------------------------------------------------------------------------------------!

			call interpol_UV(U,V,inter_U,inter_V)

			uk_e = inter_U(i-1,j)*TK(i-1,j) 
			uk_p = inter_U(i,j)  *TK(i,j)
			uk_w = inter_U(i+1,j)*TK(i+1,j)

			vk_s = inter_V(i,j-1)*TK(i,j-1)
			vk_p = inter_V(i,j)  *TK(i,j)			
			vk_n = inter_V(i,j+1)*TK(i,j+1)
			
			duk_dx = Deri(uk_e,uk_p,uk_w, dxu(i-1), dxu(i))

			dvk_dy = Deri(vk_s,vk_p,vk_n, DYV(j-1), DYV(j)) 

			DTK_dt(i,j) = duk_dx + dvk_dy

!-----------------------------------------------------------------------------------------------!
		end do	 
	end do
			DTK_dt(:,Ny) = DTK_dt(:,Ny-1); DTK_dt(:,1) = DTK_dt(:,2)
			DTK_dt(1,:)  = DTK_dt(2,:)	; DTK_dt(Nx,:) = DTK_dt(Nx-1,:) 

			Pk(:,Ny) = Pk(:,Ny-1); Pk(:,1) = Pk(:,2)
			Pk(1,:)  = Pk(2,:)	; Pk(Nx,:) = Pk(Nx-1,:) 

			epsi(:,Ny) = epsi(:,Ny-1); epsi(:,1) = epsi(:,2)
			epsi(1,:)  = epsi(2,:)	; epsi(Nx,:)= epsi(Nx-1,:) 

			div_T(:,Ny) = div_T(:,Ny-1); div_T(:,1) = div_T(:,2)
			div_T(1,:)  = div_T(2,:)	; div_T(Nx,:) = div_T(Nx-1,:) 			

			vorticity_z(:,Ny) = vorticity_z(:,Ny-1); vorticity_z(:,1) = vorticity_z(:,2)
			vorticity_z(1,:)  = vorticity_z(2,:)	; vorticity_z(Nx,:) = vorticity_z(Nx-1,:)			
end subroutine energy_budget



!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing for Matlab				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine write_matlab()

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Nx,Ny) :: PSIP,inter_U,inter_V,inter_txyt,Pk,epsi,div_T,DTK_Dt,vorticity_z,ETT
	character(len=8) :: case_jet_string


	case_jet_string="diff_cavity"


	open(unit=1, file ='./outputdata/matlab/x' //trim(case_jet_string)// '.dat')
	open(unit=2, file ='./outputdata/matlab/y' //trim(case_jet_string)// '.dat')

	open(unit=33, file ='./outputdata/matlab/T' //trim(case_jet_string)// '.dat')
	open(unit=34, file ='./outputdata/matlab/C' //trim(case_jet_string)// '.dat')

	open(unit=4, file ='./outputdata/matlab/stream.dat')

	open(unit=55, file ='./outputdata/matlab/u' //trim(case_jet_string)// '.dat')
	open(unit=56, file ='./outputdata/matlab/u_re' //trim(case_jet_string)// '.dat')	
	open(unit=66, file ='./outputdata/matlab/v'//trim(case_jet_string)// '.dat')
	open(unit=67, file ='./outputdata/matlab/v_re' //trim(case_jet_string)// '.dat')	

	open(unit=7, file ='./outputdata/matlab/TK.dat')
	open(unit=8, file ='./outputdata/matlab/e.dat')

	open(unit=10, file ='./outputdata/matlab/Ruu' //trim(case_jet_string)// '.dat')
	open(unit=11, file ='./outputdata/matlab/Rvv' //trim(case_jet_string)// '.dat')
	open(unit=12, file ='./outputdata/matlab/Ruv' //trim(case_jet_string)// '.dat')

	open(unit=13, file ='./outputdata/matlab/Pk' //trim(case_jet_string)// '.dat')
	open(unit=14, file ='./outputdata/matlab/epsi' //trim(case_jet_string)// '.dat')
	open(unit=15, file ='./outputdata/matlab/div_T' //trim(case_jet_string)// '.dat')
	open(unit=16, file ='./outputdata/matlab/DTK_Dt' //trim(case_jet_string)// '.dat')

!	open(unit=17, file ='./outputdata/matlab/vorticity_z.dat')

	open(unit=18, file ='./outputdata/matlab/rho' //trim(case_jet_string)// '.dat')
	open(unit=19, file ='./outputdata/matlab/Energy_TT' //trim(case_jet_string)// '.dat')
	open(unit=88, file ='./outputdata/matlab/pressure' //trim(case_jet_string)// '.dat')

	open(unit=99, file ='./outputdata/matlab/ite_coupling' //trim(case_jet_string)// '.dat',STATUS='REPLACE')
	write(99,*) ite_coupling
	close(99)

    CALL STREAMLINES(DYP,U,PSIP,RHO)
	CALL interpol_UV(U,V,inter_U,inter_V)
	CALL interpol_txyt(inter_txyt)
	CALL energy_budget(Pk,epsi,div_T,DTK_Dt,vorticity_z)
	CALL Energy_transport(ETT)	


	do i=1,Nx
		write(1,*) X(i)
	end do

	do j=1,Ny
		write(2,*) y(j)
	end do

	do i=1,Nx
		write(33,*) T(i,:)
	end do

	do i=1,Nx
		write(34,*) C(i,:)
	end do

	do i=1,Nx
		write(4,*) PSIP(i,:)
	end do

	do i=1,Nx
		write(55,*) inter_U(i,:)
	end do

	do i=1,Nx
		write(56,*) U(i,:)
	end do

	do i=1,Nx
		write(66,*) inter_V(i,:)
	end do

	do i=1,Nx
		write(67,*) V(i,:)
	end do


	do i=1,Nx
		write(7,*) TK(i,:)
	end do

	do i=1,Nx
		write(8,*) e(i,:)
	end do

	do i=1,Nx
		write(10,*) Txxt(i,:)/RHO(i,:)
	end do

	do i=1,Nx
		write(11,*) Tyyt(i,:)/RHO(i,:)
	end do

	do i=1,Nx
		write(12,*) inter_txyt(i,:)/RHO(i,:)
	end do

	do i=1,Nx
		write(13,*) Pk(i,:)
	end do

	do i=1,Nx
		write(14,*) epsi(i,:)
	end do

	do i=1,Nx
		write(15,*) div_T(i,:)
	end do

	do i=1,Nx
		write(16,*) DTK_Dt(i,:)
	end do

	do i=1,Nx
		write(17,*) vorticity_z(i,:)
	end do

	do i=1,Nx
		write(18,*) rho(i,:)
	end do

	do i=1,Nx
		write(88,*) P(i,:)
	end do

	do i=1,Nx
		write(19,*) ETT(i,:)
	end do

	close(1)
    close(2)
	close(33)
	close(34)
    close(4)
    close(55)
    close(56)    
    close(66)
    close(67)
    close(7)
	close(8)
    close(10)
	close(11)
    close(12)
    close(13)
    close(14)
	close(15)
    close(16) 
    close(17)   
    close(18)
    close(88) 
    close(19)
end subroutine write_matlab

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!				Varaible Fields Data Writing for TECPLOT				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine write_tecplot(U,V,P,T,C,RHO,MuT,TK,e)

	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,C,PSIP,RHO,MuT,TK,e,inter_U,inter_V

	open(unit=100, file ='./outputdata/tecplot.dat')


      	CALL STREAMLINES(DYP,U,PSIP,RHO)
	CALL interpol_UV(U,V,inter_U,inter_V)

	write(100,*)"Title= RANS"
 	write(100,*)"Variables=X,Y,u,v,PSIP,T,C,P,TK,Mut"	
	write(100,*)"Zone I=",Nx, "J=" ,Ny, "F=point"
	write(100,*)
	
	do  j=1,Ny
		do i=1,Nx
	
		Write(100,"(200F25.15)") X(i),Y(j),inter_U(i,j),inter_V(i,j),PSIP(i,j),T(i,j),C(i,j),	&
		P(i,j),TK(i,j),e(i,j),Mut(i,j)

		end do
	end do


	close(100)

end subroutine write_tecplot



subroutine write_vtk_old(U,V,P,T,b,iter_t,DYP)

	implicit none

	character(len = 27) :: str, tmp
	integer :: i,j,iter_t
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,V,P,T,b,PSIP


	write (tmp,'(I5.5)') iter_t
	str = 'outputdata/output_'//trim(tmp)//'.vtk'
	open(unit=1, file = str)

	write ( 1,'(a)' )'# vtk DataFile Version 3.0'
	write ( 1,'(a)' ) 'Non-uniform Rectilinear - Rectilinear Grid'
	write ( 1,'(a)' ) 'ASCII'
	write ( 1,'(a)' ) 'DATASET RECTILINEAR_GRID'
	write ( 1,'(a,3i4)' )'DIMENSIONS', nx,ny,1
	write ( 1,'(a13,1i4,a6)' ) 'X_COORDINATES',nx,' float'
	write ( 1,'(f14.6)' ) (x(i),i=1,nx)
	write ( 1,'(a13,1i4,a6)' ) 'Y_COORDINATES',ny,' float'
	write ( 1,'(f14.6)' ) (y(j),j=1,ny)
	write ( 1,'(a13,1i4,a6)' ) 'Z_COORDINATES',1,' float'
	write ( 1,'(f14.6)' ) 1.0

	write ( 1, '(a, I16)' ) 'CELL_DATA ', (nx-1)*(ny-1)*1

	write ( 1, '(a)' ) 'SCALARS P float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((P(i,j),i=1,(nx-1)),j=1,(ny-1))
	
	write ( 1, '(a)' ) 'SCALARS T float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((T(i,j),i=1,(nx-1)),j=1,(ny-1))

      	CALL STREAMLINES(DYP,U,PSIP,RHO)

	write ( 1, '(a)' ) 'SCALARS Str float'
	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
	write ( 1, '(f22.14)' ) ((PSIP(i,j),i=1,(nx-1)),j=1,(ny-1))

	write ( 1, '(a)' ) 'VECTORS velocity float'
	do j=1,ny-1
		do i=1,nx-1
			write ( 1, "(f22.14,f22.14,f22.14)" ) U(i,j), V(i,j), 0.0D+00
		end do
	end do

!	write ( 1, '(a)' ) 'SCALARS u float'
!	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
!	write ( 1, '(f22.14)' ) ((U(i,j),i=1,(nx-1)),j=1,(ny-1))

!	write ( 1, '(a)' ) 'SCALARS v float'
!	write ( 1, '(a)' ) 'LOOKUP_TABLE default'
!	write ( 1, '(f22.14)' ) ((V(i,j),i=1,(nx-1)),j=1,(ny-1))   



	close(1)

!	TIMECODE=TIMEF()

end subroutine write_vtk_old
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			STREAM LINES (PSIP = LINEAS DE CORRIENTE) 				!
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

SUBROUTINE STREAMLINES(DYP,U,PSIP,RHO) 
                
	implicit none

	integer :: i,j
	real(KIND = DP), dimension(Ny) :: DYP
	real(KIND = DP), dimension(Nx,Ny) :: U,RHO,PSI,PSIP
	
	DO I=1,Nx
	DO J=1,Ny
		PSI(I,J)=0.0D+00
		PSIP(I,J)=0.0D+00
	end do
	end do
	
!  NOTA: CALCULO DE PSI SOBRE LA MALLA DE LA VELOCIDAD "V"

        DO I=2,NX-1
        DO J=2,NY-2
 
              PSI(I,J)=(((U(I,J)+U(I-1,J))/2.0D+00)*DYP(J) * RHO(i,j)*1000 )+PSI(I,J-1)

	END DO
	END DO

!  NOTA: CALCULO DE PSI SOBRE LA MALLA PRINCIPAL (T,P)

        DO I=2,NX-1
        DO J=3,NY-2 

                  PSIP(I,J)=(PSI(I,J)+PSI(I,J-1))/2.0D+00

	END DO
	END DO

!     LADO SUR
        DO I=2,NX-1
               J=2 
                  PSIP(I,J)=(PSI(I,J)+PSI(I,J-1))/2.0D+00
	END DO

!     LADO NORTE
        DO I=2,NX-1
               J=NY-1 
                  PSIP(I,J)=(PSI(I,J)+PSI(I,J-1))/2.0D+00
	END DO
  
      RETURN            
END SUBROUTINE STREAMLINES

 
end module visual
