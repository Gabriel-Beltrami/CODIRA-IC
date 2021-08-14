program CONDIRA

	use user
	use grids
	use cfdlib
	use visual 

	implicit none

	call CPU_TIME(Time_start)

	call GRID1d (DXP,DX,X,Hx)
	call start(T,Tref)									! Initial values for T
      	call COEF_T(DXP,DX,X,							&
			ap,ae,aw,b)


	
	!--------------------------------------- Jacobi's Method -------------------------------------------------------!
	do iter=2,itermax	
	
	!--------------------------------------- Compute scalar  -----------------------------------------------------!
			DO I=2,NX-1
					   T(I,iter)=(aw(I)*T(I-1,iter-1)+ae(I)*T(I+1,iter-1)+b(I))/ap(I)
					   
			END DO         		
		
	!--------------------------------------- check convergence  -----------------------------------------------------!

	call residual(T,ap,ae,aw,b,Res_T,Rmax_T)

		! Monitoring residual values during execution
		if (mod(iter,1)==0) then
			call print_Residual(iter,Res_T,Rmax_T)
		end if
		
		if ( converegence(Res_T,iter)  ) then
			exit
                end if

	end do
!--------------------------------------------- 		 END Loop 	-----------------------------------------------------!
	
	call write_vtk(T,b)

	call CPU_TIME(Time_finish)
        write(*,*) (time_finish-time_start)

	close(2)
	close(3)
	
end program CONDIRA
