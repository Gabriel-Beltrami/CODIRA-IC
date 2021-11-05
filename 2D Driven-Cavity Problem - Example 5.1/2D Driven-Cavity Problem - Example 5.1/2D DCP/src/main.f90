program CONDIRA

	use user
	use grids
	use cfdlib
	use visual 

	implicit none

	call CPU_TIME (Time_start)

	call GRID2d (x,y,h,dom_length)
	call start (nu,Re,u_final,v_final,p_final,u,u_star,d_e,u_new,v, &
			v_star,d_n,v_new,p,p_star,pc,b,p_new,error,iterations)
		
	!--------------------------------------- SIMPLE -------------------------------------------------------!
	
	do while (error > error_req)
	   
	   	! Compute U
		call coef_U(u_E,u_W,v_N,v_S,a_E1,a_E2,a_E3,	&
	        		a_W,a_N2,a_S,nu,h,u,u_star,d_e,	&
					v,p)
	    	! Compute V 
		call coef_V(u_E,u_W,v_N,v_S,a_E2,	&
	         a_W,a_N1,a_N2,a_N3,a_S,nu,h,u,d_n,v,v_star,p)
	    	
	    	! Compute Pc
		call coef_Pc(a_E2,a_W,h,a_N2,a_S,a_P,	&
				u,u_star,d_e,v,v_star,d_n, &
					p,b,pc)
	    
		call correct (alpha,alpha_p,u_star,d_e,u_new,	&
				v_star,d_n,v_new,p,pc,p_new)
			
		!---------------------------------------END SIMPLE/SIMPLEC-----------------------------------------------------!

		!--------------------------------------- check convergence  -----------------------------------------------------!
		
		! Monitoring residual values during execution    
		call residual(error,b)
	      	if (mod(iterations,1000)==0) then
			call print_Residual(iterations,error)
		end if
	   
	   	! Update corrected values
		call update (iterations,u,u_new,v,v_new,p,p_new)	    
	    
	end do

!--------------------------------------------- END Loop -----------------------------------------------------!

!--------------------------------------- Compute final vectors -----------------------------------------------------!

	! Staggered variables to collocated variables
	do i = 1,n_points
	    do j = 1,n_points
		u_final(i,j) = 0.5*(u(i,j) + u(i+1,j))
		v_final(i,j) = 0.5*(v(i,j) + v(i,j+1))
		p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1))
	    end do
	end do

!--------------------------------------- vtk file generation -----------------------------------------------------!

	call write_vtk(u_final,v_final)

	call CPU_TIME(Time_finish)
        write(*,*) (time_finish-time_start)

	close(2)
	close(3)
	
end program CONDIRA
