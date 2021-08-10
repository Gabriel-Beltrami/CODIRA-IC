module solver

implicit none

contains

subroutine jacobi1d(coef, vector, b, nx, residual, itermax)

integer :: size, i, j, itermax, nx
real, allocatable :: coef(:,:), vector(:), residual(:,:), b(:)
real :: current

size = nx+2

do i = 1,itermax
  do j = 1,size
    current = (b(j) + coef(j, 3)*vector(j-1) + coef(j,4)*vector(j+1))/coef(j,1)
    residual(i,j) = abs(current - coef(j,2))
    coef(j,2) = current
    end do
  vector(:) = coef(:,2)
  residual(i,size+1) = sum(residual(i,1:size))
  write(*,*) vector
  write(*,*) residual(i, size+1)
  end do
end subroutine jacobi1d
end module solver
