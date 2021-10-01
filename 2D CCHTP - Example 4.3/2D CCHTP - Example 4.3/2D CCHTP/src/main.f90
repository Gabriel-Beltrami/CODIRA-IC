program CONDIRA

	use user
	use grids
	use cfdlib
	use visual 

	implicit none

	call CPU_TIME(Time_start)

	call GRID2d (DXP,DYP,DX,DY,X,Y,Hx,Hy)
	call start(T,Ti,Tref)									! Initial values for T
      	call Mass_Flux(u,v,Rho,Fep,Fnp)
      	call COEF_T(DXP,DYP,DX,DY,X,Y,							&
			ap,ae,aw,an,as,b,Fep,Fnp)
		
	!--------------------------------------- Jacobi's Method -------------------------------------------------------!
	
	do iter=2,itermax	
	
	!--------------------------------------- Compute scalar  -----------------------------------------------------!
		
		DO J=2,NY-1
			DO I=2,NX-1
					   T(I,J)=(aw(I,J)*Ti(I-1,J)+ae(I,J)*Ti(I+1,J)+an(I,J)*Ti(I,J+1)+	&
					   as(I,J)*Ti(I,J-1)+b(I,J))/ap(I,J)
					   
			END DO
		END DO         		

		!DO J=2,NY-1
		!	DO I=2,NX-1
		!			T(I,J)=RELAXT*(T(I,J)+Ti(I,J))
						   
		!	END DO
		!END DO
		
	!--------------------------------------- check convergence  -----------------------------------------------------!

	call residual(T,ap,ae,aw,an,as,b,Res_T,Rmax_T)

		! Monitoring residual values during execution
		if (mod(iter,1)==0) then
			call print_Residual(iter,Res_T,Rmax_T)
		end if
		
		if ( converegence(Res_T,iter)  ) then
			exit
                end if
	
		DO J=2,NY-1
			DO I=2,NX-1
					Ti(I,J)=T(I,J)
						   
			END DO
		END DO

	end do
!--------------------------------------------- 		 END Loop 	-----------------------------------------------------!
	
	call write_vtk(T,b)

	call CPU_TIME(Time_finish)
        write(*,*) (time_finish-time_start)

	close(2)
	close(3)
	
end program CONDIRA
