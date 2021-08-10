module solver

implicit none

contains

subroutine jacobi3d(coef, vector, b, nx, ny, nz, residual, itermax)

integer :: size, i, j, itermax, nx, ny, nz, k, id, l
real, allocatable :: coef(:,:), vector(:), residual(:,:), b(:)
real :: current

size = (nx+2)*(ny+2)+(nz+2)

do i = 1,itermax
  residual(i, size+1) = 0
  do j = 1,nx+2
    do k = 1, ny+2
      do l = 1, nz+2
    ID =  l + (k-1)*(Nz+2) + (Nz+2)*(Ny+2)*(j-1)
    current = b(id) + coef(id, 3) * vector(id- (Ny + 2)*(nz+2)) + coef(id,4) * vector(id+(Ny + 2)*(nz+2))
    current = current + coef(id,5) * vector(id - (nz+2)) + coef(id, 6) * vector(id + (nz+2))
    current = (current + coef(id,7) * vector(id - 1) + coef(id, 8) * vector(id + 1))/coef(id,1)
    residual(i,id) = abs(current - coef(id,2))
    residual(i, size + 1) = residual(i, size + 1) + residual(i, Id)
    coef(id,2) = current
    end do
    end do
    end do
  vector(:) = coef(:,2)
  !write(*,*) vector
  write(*,"(F20.15)") residual(i, size+1)
  end do
  !write(*,*) vector
  write(*,"(F20.15)") residual(itermax, size+1)
end subroutine jacobi3d
end module solver
