program CONDIRA

	use user
	use grids
	use cfdlib
	use visual

	implicit none
	call GRID2d (DXP,DYP,DXU,DYV,X,Y,XU,YV,Hx,Hy)
	call start(U,V,T,T_ref, Told)						! Initaial values for U,V and P

	!--------------------------------------- SIMPLE/SIMPLEC -------------------------------------------------------!
	do iter=1,itermax

		call update_U (Ux,U)						! Replace Ux=U
		call update_V (Vx,V)						! Replace Vx=V

				DO I=1,NX-1
					   U(I,1)=U_const                                         !      "FRONTERA SUR"
					   U(I,NY)=U_const                                        !      "FRONTERA NORTE"
				END DO
				DO J=2,NY-1
					   U(1,J)=U_const                                         !     "FRONTERA OESTE"
					   U(NX-1,J)=U_const                                       !     "FRONTERA ESTE"
				END DO

				DO I=1,NX
					   V(I,1)=V_const                                           !    "FRONTERA SUR"
					   V(I,NY-1)=V_const
				END DO
				DO J=2,NY-2
					   V(1,J)=V_const                                           !    "FRONTERA OESTE"
					   V(NX,J)=V_const                                          !    "FRONTERA ESTE"
        			END DO


	!---------------------------------------END SIMPLE/SIMPLEC-----------------------------------------------------!

	!--------------------------------------- Compute scalars  -----------------------------------------------------!
		   ! Compute T
	      	call COEF_T(DXP,DYP,DXU,DYV,X,Y,XU,YV,			&
							CONTERo,Cpo,RHOo,U,V,T,ap,ae,aw,as,an,b,Res_T,Rmax_T, dt, T_s, T_n, T_e, T_w)


				do i = 1, Nx
					do j = 1, ny
						T(i,j) = (aw(i,j)*Told(i-1,j) + ae(i,j)*Told(i+1,j) + an(i,j)*Told(i,j+1) + as(i,j)*Told(i,j-1) + b(i,j))/ap(i,j)
						end do
						end do
				do i = 1, nx
					do j = 1, ny
						Told(i,j) = relaxT*(T(i,j)+Told(i,j))
						end do
						end do










	!--------------------------------------- check convergence  -----------------------------------------------------!

		   ! Monitoring residual values during execution
		if (mod(iter,10)==0) then
			call print_Residual(iter,Res_T,Rmax_T)
		end if

		if ( converegence(Res_T,iter) .and. iter .gt. 2 ) then
			write(*,*) ' Convergence - Iter', iter
			!write(*,*) ' Res = ', Res_T
			!write(*,*) T
			call write_vtk(U,V,T)
			exit
		   end if

	end do
!--------------------------------------------- 		 END Loop 	-----------------------------------------------------!
!write(*,*) ' Res = ', Res_T
!write(*,*) T

	!call write_vtk(U,V,P,T,b,iter_t,TK,e,MuT)

	!close(2)
	!close(3)

end program CONDIRA
