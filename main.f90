!    Copyright (C) 2022 Jan Mateu Armengol 
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.

program CONDIRA

	use user
	use grids
	use cfdlib
	use visual
	use mpi 

	implicit none

	real(kind = dp) :: Sch_ma, Th_ma, U_ma, emiss_ma, C_ma
	integer :: rank, ierr, np, zed




	real, dimension(:), allocatable :: Sch_par, Th_par, U_par, emiss_par, C_par
	integer :: len_Sch, len_Th, len_U, len_emiss, len_C

	call CPU_TIME(Time_start)

	open(unit = 99, file = "SamplesSch.txt", status = 'old', action = 'read')
	read(99,*) len_Sch
	allocate(Sch_par(len_Sch))
	do i = 1, len_Sch
		read(99,*) Sch_par(i)
		end do

	open(unit = 98, file = "SamplesSkinTemp.txt", status = 'old', action = 'read')
	read(98,*) len_Th
	allocate(Th_par(len_Th))
	do i = 1, len_Th
		read(98,*) Th_par(i)
	end do

		open(unit = 97, file = "SamplesWdSpd.txt", status = 'old', action = 'read')
	read(97,*) len_U
	allocate(U_par(len_U))
	do i = 1, len_U
		read(97,*) U_par(i)
	end do

		open(unit = 96, file = "SamplesEmiss.txt", status = 'old', action = 'read')
	read(96,*) len_emiss
	allocate(emiss_par(len_emiss))
	do i = 1, len_emiss
		read(96,*) emiss_par(i)
	end do
	
		open(unit = 95, file = "SamplesBackCon.txt", status = 'old', action = 'read')
	read(95,*) len_C
	allocate(C_par(len_C))
	do i = 1, len_C
		read(95,*) C_par(i)
	end do

	call mpi_init(ierr)
	call mpi_comm_size(MPI_COMM_WORLD, np, ierr)
	call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

	do zed = rank+1, size(Sch_par), np
	Sch_ma = Sch_par(zed)
	Th_ma = Th_par(zed)
	U_ma = U_par(zed)
	emiss_ma = emiss_par(zed)
	C_ma = C_par(zed)

	write(*,*) Sch_ma, Th_ma, U_ma, emiss_ma, C_ma, zed, rank



	call dependent_var(DeltaT,Hx,Hy,dt,Tref,CONTERo,Po,Muo,U_bg,Sch_ma,Th_ma,U_ma,emiss_ma,C_ma) ! Compute dependent variables
	!call GRID_half_jet(DXP,DYP,DXU,DYV,X,Y,XU,YV)					! Grid generation
	call GRID2d (DXP,DYP,DXU,DYV,X,Y,XU,YV,Hx,Hy)
	call start(U,V,P,T,C,Tref,rho,Fep,Fnp,TKx,e)						! Initaial values for U,V and P
	call write_vtk(U,V,P,T,C,b,iter_t,zed,TK,e,MuT,Sch_ma,Th_ma,U_ma,emiss_ma,C_ma)

	!--------------------------------------- SIMPLE/SIMPLEC -------------------------------------------------------!
	do iter=1,itermax

		call update_U (Ux,U)						! Replace Ux=U
		call update_V (Vx,V)						! Replace Vx=V		
		call update_P (Px,P)						! Replace Px=P
		call update_P (Tk,Tkx)						! Replace Tkx=Tk		

		   ! Detemine the local thermal properties
        CALL PROPERTIES(Mu,CONTER,CP,CONTERo,T,Muo,CPo,RHO,Po,Pmed,Hx,Hy,DXP,DYP,Tref,RHOo,MuT,e,TK)


		if (rans_model .eq. 0) then
			CALL stressLT(Mu,U,V,DXP,DYP,DYV,DXU,TxxL,TyyL,TxyL, TxxT, TyyT, TxyT, Tk, MuT, RHO)
			TxxT(:,:) = 0.0D+0 ; TxyT(:,:) = 0.0D+0 ; TyyT(:,:) = 0.0D+0 ; MuT(:,:) = 0.0D+0 ;
		else
			CALL stressLT(Mu,U,V,DXP,DYP,DYV,DXU,TxxL,TyyL,TxyL, TxxT, TyyT, TxyT, Tk, MuT, RHO)
		end if


		   ! Compute U 	
	    call coef_U(DXP,DYP,DXU,DYV,X,Y,XU,YV,				&
			Ux,Vx,Px,							&
			du,ae,aw,as,an,b,RES_U,Rmax_U,mu,apu,dt,RHO,			&
			TxxL, TyyL, TxyL, Fep, Fnp, MUt, TxxT, TyyT, TxyT,APUNB,U_bg)

		call LBL_adi(1,1,(Nx-1),Ny,apu,ae,aw,an,as,b,Ux,npas_U,RELAXU)
!		call LBL_x(1,1,(Nx-1),Ny,apu,ae,aw,an,as,b,Ux,npas_U,RELAXU)

		   ! Compute V 	
	    call COEF_V(DXP,DYP,DXU,DYV,X,Y,XU,YV,				&
		 	 U,Vx,Px,dv,ae,aw,as,an,b,					&
		 	 RES_V,Rmax_V,mu,apv,dt,Tref,RHO,RHOmed,			&
		  	TxxL, TyyL, TxyL, Fep, Fnp, MuT, TxxT, TyyT, TxyT,APVNB)

	    call LBL_adi(1,1,Nx,(Ny-1),apv,AE,AW,AN,AS,B,Vx,npas_V,RELAXV)
			!	call LGS_adi(1,1,Nx,(Ny-1),apv,AE,AW,AN,AS,B,Vx,npas_V,RELAXV)
!	    call LBL_x(1,1,Nx,(Ny-1),apv,AE,AW,AN,AS,B,Vx,npas_V,RELAXV)

		   ! Compute P 

		call COEF_P(DXP,DYP,Ux,Vx,du,dv,ap,ae,aw,as,an,b,RES_P,Rmax_P,RHO,apu,apv,apunb,apvnb)	
	    call reset(Pc)							
		call LBL_ADI(1,1,Nx,Ny,AP,AE,AW,AN,AS,B,Pc,npas_P,RelaxP)

		   ! Update corrected values
		call update_U (U,Ux)						
		call update_V (V,Vx)
		call correct (U,Ux,V,Vx,Pc,P,Px,du,dv,Fep,Fnp,rho)


	!---------------------------------------END SIMPLE/SIMPLEC-----------------------------------------------------!

	!--------------------------------------- Compute scalars  -----------------------------------------------------!
		   ! Compute T 	
	      	call COEF_T(DXP,DYP,DXU,DYV,X,Y,XU,YV,							&
			CONTER,Cp,RHO,U,V,T,ap,ae,aw,as,an,b,Res_T,Rmax_T,dt,Fep,Fnp,MUt)

!	      	call LBL_ADI(1,1,Nx,Ny,ap,AE,AW,AN,AS,B,T,npas_T,RELAXT)
	      	call LBL_x(1,1,Nx,Ny,ap,AE,AW,AN,AS,B,T,npas_T,RELAXT)
	      	
	      	   ! Compute C 	
	      	call COEF_C(DXP,DYP,DXU,DYV,X,Y,XU,YV,							&
			CONTER,Cp,RHO,U,V,C,ap,ae,aw,as,an,b,Res_C,Rmax_C,dt,Fep,Fnp,MUt)

!	      	call LBL_ADI(1,1,Nx,Ny,ap,AE,AW,AN,AS,B,C,npas_C,RELAXC)
	      	call LBL_x(1,1,Nx,Ny,ap,AE,AW,AN,AS,B,C,npas_C,RELAXC)

		if (rans_model .ne. 0) then

		   ! Compute TK	
			call COEF_Tk(DXP,DYP,DXU,DYV,X,Y,XU,YV,							&
				Mu,Mut,RHO,U,V,TKx,e,ap,ae,aw,as,an,b,Res_Tk,Rmax_Tk,dt,Fep,Fnp,T,Txxt,Tyyt,Txyt,Pk,Gk)
			call LBL_ADI(1,1,Nx,Ny,ap,AE,AW,AN,AS,B,Tkx,npas_T,RELAXTk)

		   ! Compute e	
			call COEF_e(DXP,DYP,DXU,DYV,X,Y,XU,YV,			&
				Mu,MuT,RHO,U,V,Tk,e,ap,ae,aw,as,an,b,Res_e,Rmax_e,dt,Fep,Fnp,T,Txxt,Tyyt,Txyt,Pk,Gk)	
			call LBL_ADI(1,1,Nx,Ny,ap,AE,AW,AN,AS,B,e,npas_T,RELAXe)

		end if

	!--------------------------------------- check convergence  -----------------------------------------------------!

		   ! Monitoring residual values during execution
		if (mod(iter,1000)==0) then
			call print_Residual (iter,iter_t,Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Res_C,Rmax_T,Rmax_C,	&
				Res_Tk,Rmax_Tk,Res_e,Rmax_e)
		end if

		if ( converegence(Res_U,Res_V,Res_P,Res_T,Res_C,Res_e,Res_TK,iter)  ) then
			exit
		end if

	end do
!--------------------------------------------- 		 END Loop 	-----------------------------------------------------!
	call print_Residual (iter,iter_t,Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Res_C,Rmax_T,Rmax_C,	&
				Res_Tk,Rmax_Tk,Res_e,Rmax_e)
	call write_vtk(U,V,P,T,C,b,iter_t,zed,TK,e,MuT,Sch_ma,Th_ma,U_ma,emiss_ma,C_ma)
	!call write_matlab()
	call write_matlab(zed,iter,iter_t,Res_U,Res_V,Res_P,Rmax_U,Rmax_V,Rmax_P,Res_T,Res_C,Rmax_T,Rmax_C,	&
					Res_Tk,Rmax_Tk,Res_e,Rmax_e)

end do

! Deallocate arrays
  deallocate(Sch_par)
  deallocate(Th_par)
  deallocate(U_par)
  deallocate(emiss_par)
  deallocate(C_par)

call mpi_finalize(ierr)

	call CPU_TIME(Time_finish)
    write(*,*) "computational time",(time_finish-time_start)

end program CONDIRA


