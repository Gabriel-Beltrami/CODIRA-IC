program CONDIRA

	use user
	use grids
	use cfdlib
	use visual 

	implicit none

	call CPU_TIME(Time_start)

	iter_t=1
	
	call GRID1d (DXP,DX,X,Hx)								! Grid 1-dimensional
	call start(T,X)									! Initial values for T
      	call COEF_T(DXP,DX,X,							&		! Coefficients of T
			ap,apo,ae,aw,b,T,iter_t)
	
	!--------------------------------------- Implicit Method -------------------------------------------------------!

	do iter_t=2,itermax_t		
			
		if (iter>2) then
			
			call COEF_T(DXP,DX,X,							&
					ap,apo,ae,aw,b,T,iter_t)
		end if
			
		do iter=2,itermax

	!--------------------------------------- Compute scalar  -----------------------------------------------------!

			DO I=2,NX-1
					   T(I,iter,iter_t)=(aw(I)*T(I-1,iter-1,iter_t)+ae(I)*T(I+1,iter-1,iter_t)+b(I))/ap(I)
					   
			END DO         		
		
	!--------------------------------------- check convergence  -----------------------------------------------------!

			call residual(T,ap,apo,ae,aw,b,Res_T,Rmax_T)

			 !Monitoring residual values during execution
			if (mod(iter,10)==0) then
				call print_Residual(iter,iter_t,Res_T,Rmax_T)
			end if
			
			if ( converegence(Res_T,iter)  ) then	
				exit
		       end if
			
		end do		
		
		write(*,*) 'T1:',T(1,iter,iter_t)
		write(*,*) 'T2:',T(2,iter,iter_t)
		write(*,*) 'T3:',T(3,iter,iter_t)
		write(*,*) 'T4:',T(4,iter,iter_t)
		write(*,*) 'T5:',T(5,iter,iter_t)
		write(*,*) 'T6:',T(6,iter,iter_t)
		write(*,*) 'T7:',T(7,iter,iter_t)		
		
		if (iter_t<(itermax_t-1)) then
			DO I=2,NX-1
						   T(I,1,iter_t+1)=T(I,iter,iter_t)
			END DO
		
		end if
	
	end do
	
!--------------------------------------------- 		 END Loop 	-----------------------------------------------------!
	
	call write_vtk(T,b)

	call CPU_TIME(Time_finish)
        write(*,*) (time_finish-time_start)

	close(2)
	close(3)
	
end program CONDIRA
