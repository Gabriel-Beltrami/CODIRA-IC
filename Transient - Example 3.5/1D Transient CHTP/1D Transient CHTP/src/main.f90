program CONDIRA

	use user
	use grids
	use cfdlib
	use visual 

	implicit none

	call CPU_TIME(Time_start)

	iter_t=1
	
	call GRID1d (DXP,DX,X,Hx)								! Grid 1-dimensional
	call start(T,Ti,Tt,X)									! Initial values for T
      	call COEF_T(DXP,DX,X,							&		! Coefficients of T
			ap,apo,ae,aw,b,T,Ti,Tt,iter_t)
	
	!--------------------------------------- Implicit Method -------------------------------------------------------!

	do iter_t=2,itermax_t		
			
		if (iter>2) then
			
			call COEF_T(DXP,DX,X,							&
					ap,apo,ae,aw,b,T,Ti,Tt,iter_t)
		end if
			
		do iter=2,itermax

	!--------------------------------------- Compute scalar  -----------------------------------------------------!

			DO I=2,NX-1
					   T(I)=(aw(I)*Ti(I-1)+ae(I)*Ti(I+1)+b(I))/ap(I)
					   
			END DO         		
		
	!--------------------------------------- check convergence  -----------------------------------------------------!

			call residual(T,Ti,Tt,ap,apo,ae,aw,b,Res_T,Rmax_T)

			 !Monitoring residual values during execution
			if (mod(iter,10)==0) then
				call print_Residual(iter,iter_t,Res_T,Rmax_T)
			end if
			
			if ( converegence(Res_T,iter)  ) then	
				exit
		       end if
			
			
			DO I=2,NX-1
				Ti(I)=T(I)
					   
			END DO
		
		end do		
		
		write(*,*) 'T1:',T(1)
		write(*,*) 'T2:',T(2)
		write(*,*) 'T3:',T(3)
		write(*,*) 'T4:',T(4)
		write(*,*) 'T5:',T(5)
		write(*,*) 'T6:',T(6)
		write(*,*) 'T7:',T(7)		
		
		if (iter_t<(itermax_t-1)) then
			DO I=2,NX-1
						   Tt(I)=T(I)
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
