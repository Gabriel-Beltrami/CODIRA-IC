module solver

implicit none

contains

subroutine jacobi3d(coef, vector, b, nx, ny, nz, residual, itermax, criteria)

integer :: size, i, j, itermax, nx, ny, nz, k, id, l
real, allocatable :: coef(:,:), vector(:), residual(:,:), b(:)
real :: current, c1, c2, c3, c4, c5, c6, criteria

size = (nx+2)*(ny+2)+(nz+2)

do i = 1,itermax
  residual(i, size+1) = 0
  do j = 1,nx+2
    do k = 1, ny+2
      do l = 1, nz+2
    ID =  l + (k-1)*(Nz+2) + (Nz+2)*(Ny+2)*(j-1)

    c1 = coef(id, 3) * vector(id- (Ny + 2)*(nz+2))
    c2 = coef(id, 4) * vector(id+ (Ny + 2)*(nz+2))
    c3 = coef(id, 5) * vector(id - (nz+2))
    c4 = coef(id, 6) * vector(id + (nz+2))
    c5 = coef(id, 7) * vector(id - 1)
    c6 = coef(id, 8) * vector(id + 1)

    current = (b(id) + c1 + c2 + c3 + c4 + c5 + c6)/coef(id,1)
    residual(i,id) = abs(current - coef(id,2))
    residual(i, size + 1) = residual(i, size + 1) + residual(i, Id)
    coef(id,2) = current
    end do
    end do
    end do
  vector(:) = coef(:,2)
  !write(*,*) vector
  !write(*,"(F20.15)") residual(i, size+1)

  if (criteria .gt. residual(i, size+1)) then
  write(*,*) 'Convergency Achieved'
  exit
  end if
  end do
  write(*,*) vector
  write(*,"(F20.15)") residual(itermax, size+1)
end subroutine jacobi3d






subroutine solver_transiente(coef, vector, b, nx, ny, nz, residual, ts, itermax, criteria, f)

integer :: size, i, j, itermax, nx, ny, nz, k, id, l, ts
real, allocatable :: coef(:,:,:), vector(:,:), residual(:,:,:), b(:,:)
real :: it, criteria, f, c1, c2, c3, c4, c5, c6

size = (nx+2)*(ny+2)+(nz+2)

do i = 1,itermax
  do j = 1,nx+2
    do k = 1, ny+2
      do l = 1, nz+2
    ID =  l + (k-1)*(Nz+2) + (Nz+2)*(Ny+2)*(j-1)

    c1 = coef(ts,id,  3) * (vector(ts, id - (Ny + 2)*(nz+2))*f + vector(ts-1, id - (Ny + 2)*(nz+2))*(1-f))
    c2 = coef(ts, id, 4) * (vector(ts, id + (Ny + 2)*(nz+2))*f + vector(ts-1, id + (Ny + 2)*(nz+2))*(1-f))
    c3 = coef(ts, id, 5) * (vector(ts, id - (nz+2))*f + vector(ts-1, id - (nz+2))*(1-f))
    c4 = coef(ts, id, 6) * (vector(ts, id + (nz+2))*f + vector(ts-1, id + (nz+2))*(1-f))
    c5 = coef(ts, id, 7) * (vector(ts, id - 1)*f + vector(ts-1, id -1)*(1-f))
    c6 = coef(ts, id, 8) * (vector(ts, id + 1)*f + vector(ts-1, id +1)*(1-f))


    it = (b(ts, id) + c1 + c2 + c3 + c4 + c5 + c6)/coef(ts, id,1)
    residual(i,ts, id) = abs(it - coef(ts, id,2))
    residual(i, ts, size + 1) = residual(i, ts, size + 1) + residual(i,ts,Id)
    coef(ts, id,2) = it
    end do
    end do
    end do
  if (criteria .gt. residual(i, ts, size+1)) then
    write(*,*) 'Convergency Achieved'
    exit
  end if

  vector(ts, :) = coef(ts, :,2)
  !write(*,*) vector
  !write(*,"(F20.15)") residual(i, size+1)
  end do
  write(*,*) vector(ts,:)
  !write(*,"(F20.15)") residual(itermax, size+1)
end subroutine solver_transiente



end module solver
