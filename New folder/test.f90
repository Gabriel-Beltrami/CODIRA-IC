program for_all_example

	implicit none
	
	real, allocatable, dimension(:,:) :: A
	integer :: n, i, j
	
	n = 10
	
	allocate ( A(n,n))
	print *, a
	do j = 1,n 
		forall(i=1:n) A(i,j) = i+j
	end do
	!forall(i=1:n) A(i,i) = i
	!forall(i=1:n-1) A(i, i+1) = 3
	!forall(i=2:n) A(i,i-1) = 5 
	
	print *, A
	
	deallocate (A)
	pause
	

	end program for_all_example	