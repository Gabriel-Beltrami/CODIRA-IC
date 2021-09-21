program CONDIRA

	use user
	use grids
	use cfdlib
	use visual 

	implicit none

	call CPU_TIME(Time_start)

	call GRID1d (DXP,DX,X,Hx)
	call start(T,Ti,Tref,DXP)									! Initial values for T
      	call COEF_T(DXP,DX,X,							&
			ap,ae,aw,b)

	
	!--------------------------------------- Jacobi's Method -------------------------------------------------------!
	do iter=2,itermax	
	
	!--------------------------------------- Compute scalar  -----------------------------------------------------!
			DO I=2,NX-1
					   T(I)=(aw(I)*Ti(I-1)+ae(I)*Ti(I+1)+b(I))/ap(I)
					   
			END DO         		
			
			T(Nx)=(aw(Nx)*Ti(Nx-1)+b(Nx))/ap(Nx)
	
	!--------------------------------------- check convergence  -----------------------------------------------------!

	call residual(T,ap,ae,aw,b,Res_T,Rmax_T)

		! Monitoring residual values during execution
		if (mod(iter,1)==0) then
			call print_Residual(iter,Res_T,Rmax_T)
		end if
		
		if ( converegence(Res_T,iter)  ) then
			exit
                end if
                
                DO I=2,NX
				Ti(I)=T(I)
					   
		END DO

	end do
!--------------------------------------------- 		 END Loop 	-----------------------------------------------------!
	
	call write_vtk(T,b)

	call CPU_TIME(Time_finish)
        write(*,*) (time_finish-time_start)

	close(2)
	close(3)
	
end program CONDIRA
