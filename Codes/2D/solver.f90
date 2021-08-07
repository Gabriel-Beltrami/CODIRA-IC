module solver

implicit none

contains

subroutine jacobi2d(coef, vector, b, nx, ny, residual, itermax)

integer :: size, i, j, itermax, nx, ny, k, id
real, allocatable :: coef(:,:), vector(:), residual(:,:), b(:)
real :: current

size = (nx+2)*(ny+2)

do i = 1,itermax
  do j = 1,nx+2
    do k = 1, ny+2
    ID =  j + (Ny+2)*(k-1)
    current = b(id) + coef(id, 3) * vector(id-Ny - 2) + coef(id,4) * vector(id+ny + 2)
    current = (current + coef(id,5) * vector(id - 1) + coef(id, 6) * vector(id + 1))/coef(id,1)
    residual(id,j) = abs(current - coef(id,2))
    coef(id,2) = current
    end do
    end do
  vector(:) = coef(:,2)
  residual(i,size+1) = sum(residual(i,1:size))
  write(*,*) vector
  write(*,*) residual(i, size+1)
  end do
end subroutine jacobi2d
end module solver
