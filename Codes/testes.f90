program testes

implicit none

real :: a, b, c
integer :: i, j, k
real, dimension (2,2) :: matrix
matrix = 3.
do i = 1, 2
  write(*,*) matrix(i:1, i:1)
  end do

end program testes
