module solver

implicit none

contains

subroutine jacobi2d(coef, vector, b, nx, ny, residual, itermax, criteria)

integer :: size, i, j, itermax, nx, ny, k, id
real, allocatable :: coef(:,:), vector(:), residual(:,:), b(:)
real :: current, criteria

size = (nx+2)*(ny+2)

do i = 1,itermax
  do j = 1,nx+2
    do k = 1, ny+2
    ID =  j + (Ny+2)*(k-1)
    current = b(id) + coef(id, 3) * vector(id-Ny - 2) + coef(id,4) * vector(id+ny + 2)
    current = (current + coef(id,5) * vector(id - 1) + coef(id, 6) * vector(id + 1))/coef(id,1)
    residual(i,id) = abs(current - coef(id,2))
    coef(id,2) = current



    end do
    end do


  vector(:) = coef(:,2)
  residual(i,size+1) = sum(residual(i,1:size))
  !write(*,*) vector
  !write(*,"(F20.15)") residual(i, size+1)
    if (criteria .gt. residual(i, size+1)) then
    write(*,*) 'Convergency Achieved, Iteration', i
    exit
    end if
  end do
  write(*,*) vector
  !write(*,"(F20.15)") residual(i, size+1)
end subroutine jacobi2d



subroutine solver_transiente(coef, vector, b, nx, ny, residual, ts, itermax, criteria, f)

integer :: size, i, j, itermax, nx, ny, ts, k, id
real, allocatable :: coef(:,:,:), vector(:,:), residual(:,:,:), b(:,:)
real :: it, criteria, f
real :: c1, c2, c3, c4

size = (nx+2)*(ny+2)

do i = 1,itermax
  do j = 1,nx+2
    do k = 1, ny+2
    ID =  j + (Ny+2)*(k-1)
    c1 = coef(ts,id, 3) * (vector(ts, id-Ny - 2)*f + vector(ts-1, id-Ny-2)*(1-f))
    c2 = coef(ts, id,4) * (vector(ts, id+ny + 2)*f + vector(ts-1, id+ny+2)*(1-f))
    c3 = coef(ts, id,5) * (vector(ts, id-1)*f + vector(ts-1, id-1)*(1-f))
    c4 = coef(ts, id, 6) * (vector(ts, id+1)*f + vector(ts-1, id+1)*(1-f))
    it = (c1 + c2 + c3 + c4 + b(ts, id))/coef(ts, id, 1)
    residual(i,ts, id) = abs(it - coef(ts, id,2))
    coef(ts, id,2) = it
    end do
    end do

  vector(ts, :) = coef(ts, :,2)
  residual(i,ts,size+1) = sum(residual(i,ts, 1:size))

  if (criteria .gt. residual(i, ts, size+1)) then
  !write(*,*) 'Convergency Achieved, Iteration = ', i
  exit

  end if


  end do

  !write(*,*) vector(ts,:)
  !write(*,"(F20.15)") residual(itermax, size+1)
end subroutine solver_transiente




end module solver
