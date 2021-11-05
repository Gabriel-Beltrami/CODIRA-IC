module cfdlib
	use user
	use grids
	implicit none


contains

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!			           start(U,V,P)						    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine start (nu,Re,u_final,v_final,p_final,u,u_star,d_e,u_new,v,v_star,d_n,v_new,p,p_star,pc,b,p_new,error,iterations)

	implicit none
	integer :: iterations
	real(KIND = DP)  :: nu,Re,error
	real(KIND = DP), dimension(n_points,n_points) :: u_final,v_final,p_final
	real(KIND = DP), dimension(n_points+1,n_points) :: u,u_star,d_e,u_new
	real(KIND = DP), dimension(n_points,n_points+1) :: v,v_star,d_n,v_new
	real(KIND = DP), dimension(n_points+1,n_points+1) :: p,p_star,pc,b,p_new

	! Dependent variable 
	nu=1.0D+00/Re
	
	! Initializing the variables
	
	! Final collocated variables
	u_final(n_points,n_points)=0.0D+00
	v_final(n_points,n_points)=0.0D+00
	p_final(n_points,n_points)=1.0D+00
	u_final(1,:) = 1.0D+00
	
	!Staggered variables
	u(n_points+1,n_points)=0.0D+00
	u_star(n_points+1,n_points)=0.0D+00
	d_e(n_points+1,n_points)=0.0D+00
	v(n_points,n_points+1)=0.0D+00
	v_star(n_points,n_points+1)=0.0D+00
	d_n(n_points,n_points+1)=0.0D+00
	p(n_points+1,n_points+1)=1.0D+00
	p_star(n_points+1,n_points+1)=1.0D+00
	pc(n_points+1,n_points+1)=0.0D+00
	b(n_points+1,n_points+1)=0.0D+00
	u(1,:)=2.0D+00

	u_new(n_points+1,n_points)=0.0D+00
	v_new(n_points,n_points+1)=0.0D+00
	p_new(n_points+1,n_points+1)=1.0D+00
	u_new(1,:)=2.0D+00
	
	error = 1.0D+00
	iterations = 0

    return
end subroutine start

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_U 					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine coef_U(u_E,u_W,v_N,v_S,a_E1,a_E2,a_E3,	&
	         a_W,a_N2,a_S,nu,h,u,u_star,d_e,	&
		 v,p)

	implicit none

	integer:: i,j

	real(KIND = DP)  :: 			&
		u_E,u_W,v_N,v_S,		&
		a_E1,a_E2,a_E3,		& ! A_e,a_E,a_e
		a_W,a_N2,a_S,nu,h		  ! a_N

	real(KIND = DP), dimension(n_points+1,n_points) :: u,u_star,d_e
	real(KIND = DP), dimension(n_points,n_points+1) :: v
	real(KIND = DP), dimension(n_points+1,n_points+1) :: p

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!	

	    do i = 2,n_points
		do j = 2,n_points-1
		    u_E = 0.5*(u(i,j) + u(i,j+1))
		    u_W = 0.5*(u(i,j) + u(i,j-1))
		    v_N = 0.5*(v(i-1,j) + v(i-1,j+1))
		    v_S = 0.5*(v(i,j) + v(i,j+1))
		    
		    a_E2 = -0.5*u_E*h + nu
		    a_W = 0.5*u_W*h + nu
		    a_N2 = -0.5*v_N*h + nu
		    a_S = 0.5*v_S*h + nu
		    
		    a_E3 = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu
		    
		    a_E1 = -h
		    d_e(i,j) = a_E1/a_E3
		    
		    u_star(i,j) = (a_E2*u(i,j+1) + a_W*u(i,j-1) + a_N2*u(i-1,j) + a_S*u(i+1,j))/a_E3 + d_e(i,j)*(p(i,j+1) - p(i,j))
		end do
	    end do
 
!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS                                       !
!-----------------------------------------------------------------------------------------------!

	    u_star(1,:) = 2.0D+00 - u_star(2,:)
	    u_star(n_points+1,:) = -u_star(n_points,:)
	    u_star(2:n_points,1) = 0.0D+00
	    u_star(2:n_points,n_points) = 0.0D+00
	
	return
end Subroutine coef_U

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_V 					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine coef_V(u_E,u_W,v_N,v_S,a_E2,	&
	         a_W,a_N1,a_N2,a_N3,a_S,nu,h,u,d_n,v,v_star,p)

	implicit none

	integer:: i,j

	real(KIND = DP)  :: 			&
		u_E,u_W,v_N,v_S,		&
		a_E2,a_W,nu,h,			& ! a_E
		a_N1,a_N2,a_N3,a_S		  ! A_n,a_N,a_n

	real(KIND = DP), dimension(n_points+1,n_points) :: u
	real(KIND = DP), dimension(n_points,n_points+1) :: v,v_star,d_n
	real(KIND = DP), dimension(n_points+1,n_points+1) :: p

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!	

	    ! y-momentum eq. - Interior
	    do i = 2,n_points-1
		do j = 2,n_points
		    u_E = 0.5*(u(i,j) + u(i+1,j))
		    u_W = 0.5*(u(i,j-1) + u(i+1,j-1))
		    v_N = 0.5*(v(i-1,j) + v(i,j))
		    v_S = 0.5*(v(i,j) + v(i+1,j))
		    
		    a_E2 = -0.5*u_E*h + nu
		    a_W = 0.5*u_W*h + nu
		    a_N2 = -0.5*v_N*h + nu
		    a_S = 0.5*v_S*h + nu
		    
		    a_N3 = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu
		    
		    a_N1 = -h
		    d_n(i,j) = a_N1/a_N3
		    
		    v_star(i,j) = (a_E2*v(i,j+1) + a_W*v(i,j-1) + a_N2*v(i-1,j) + a_S*v(i+1,j))/a_N3 + d_n(i,j)*(p(i,j) - p(i+1,j))
		end do
	    end do
 
!-----------------------------------------------------------------------------------------------!
!					BOUNDARY CONDITIONS                                       !
!-----------------------------------------------------------------------------------------------!

	    v_star(:,1) = -v_star(:,2)
	    v_star(:,n_points+1) = -v_star(:,n_points)
	    v_star(1,2:n_points) = 0.0D+00
	    v_star(n_points,2:n_points) = 0.0D+00
	
	return
end Subroutine coef_V

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						coef_Pc 					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine coef_Pc(a_E2,a_W,h,a_N2,a_S,a_P,	&
		    u,u_star,d_e,v,v_star,d_n, &
		    p,b,pc)

	implicit none

	integer:: i,j

	real(KIND = DP)  :: 			&
		a_E2,a_W,h,			& ! a_E
		a_N2,a_S,a_P		  	  ! a_N

	real(KIND = DP), dimension(n_points+1,n_points) :: u,u_star,d_e
	real(KIND = DP), dimension(n_points,n_points+1) :: v,v_star,d_n
	real(KIND = DP), dimension(n_points+1,n_points+1) :: p,b,pc

!-----------------------------------------------------------------------------------------------!
!						reset(Pc)					    !
!-----------------------------------------------------------------------------------------------!

		pc(1:n_points+1,1:n_points+1)=0.0D+00

!-----------------------------------------------------------------------------------------------!
!					COEFICIENT DETERMINATION				    !
!-----------------------------------------------------------------------------------------------!	
	    
	    do i = 2,n_points
		do j = 2,n_points
		    a_E2 = -d_e(i,j)*h
		    a_W = -d_e(i,j-1)*h
		    a_N2 = -d_n(i-1,j)*h
		    a_S = -d_n(i,j)*h
		    a_P = a_E2 + a_W + a_N2 + a_S
		    b(i,j) = -(u_star(i,j) - u_star(i,j-1))*h + (v_star(i,j) - v_star(i-1,j))*h
		    
		    pc(i,j) = (a_E2*pc(i,j+1) + a_W*pc(i,j-1) + a_N2*pc(i-1,j) + a_S*pc(i+1,j) + b(i,j))/a_P
		end do
	    end do
	
	return
end Subroutine coef_Pc

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!						correct					    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

Subroutine correct (alpha,alpha_p,u_star,d_e,u_new,	&
			v_star,d_n,v_new,p,pc,p_new)

	integer :: i,j
	real(KIND = DP)  :: alpha,alpha_p
	real(KIND = DP), dimension(n_points+1,n_points) :: u_star,d_e,u_new
	real(KIND = DP), dimension(n_points,n_points+1) :: v_star,d_n,v_new
	real(KIND = DP), dimension(n_points+1,n_points+1) :: p,pc,p_new
		
	   ! Correcting the pressure field
	    do i = 2,n_points
		do j = 2,n_points
		    p_new(i,j) = p(i,j) + alpha_p*pc(i,j)
		end do
	    end do
	    
	    ! Continuity eq. - Boundary
	    p_new(1,:) = p_new(2,:)
	    p_new(n_points+1,:) = p_new(n_points,:)
	    p_new(:,1) = p_new(:,2)
	    p_new(:,n_points+1) = p_new(:,n_points)
	    
	    ! Correcting the u-velocities
	    do i = 2,n_points
		do j = 2,n_points-1
		    u_new(i,j) = u_star(i,j) + alpha*d_e(i,j)*(pc(i,j+1) - pc(i,j))
		end do
	    end do
	    
	    ! x-momentum eq. - Boundary
	    u_new(1,:) = 2.0D+00 - u_new(2,:)
	    u_new(n_points+1,:) = -u_new(n_points,:)
	    u_new(2:n_points,1) = 0.0D+00
	    u_new(2:n_points,n_points) = 0.0D+00
	    
	    ! Correcting the v-velocities
	    do i = 2,n_points-1
		do j = 2,n_points
		    v_new(i,j) = v_star(i,j) + alpha*d_n(i,j)*(pc(i,j) - pc(i+1,j))
		end do
	    end do
	    
	    ! y-momentum eq. - Boundary
	    v_new(:,1) = -v_new(:,2)
	    v_new(:,n_points+1) = -v_new(:,n_points)
	    v_new(1,2:n_points) = 0.0D+00
	    v_new(n_points,2:n_points) = 0.0D+00
     	
	return
end Subroutine correct

!-----------------------------------------------------------------------------------------------!
!					RESIDUAL						!
!-----------------------------------------------------------------------------------------------!

subroutine residual(error,b)

	integer :: i,j
	real(KIND = DP) :: error
	real(KIND = DP), dimension(n_points+1,n_points+1) :: b

	! Continuity residual as error measure
	    error = 0.0D+00
	    do i = 2,n_points
		do j = 2,n_points
		    error = error + abs(b(i,j))
		end do
	    end do

	return	
end subroutine residual

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
!					update							    !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

subroutine update (iterations,u,u_new,v,v_new,p,p_new)

	implicit none
	integer:: iterations
	real(KIND = DP), dimension(n_points+1,n_points) :: u,u_new
	real(KIND = DP), dimension(n_points,n_points+1) :: v,v_new
	real(KIND = DP), dimension(n_points+1,n_points+1) :: p,p_new	
	
	    ! Finishing the iteration
	    u = u_new
	    v = v_new
	    p = p_new
	    iterations = iterations + 1

    return
end subroutine update

end module cfdlib
