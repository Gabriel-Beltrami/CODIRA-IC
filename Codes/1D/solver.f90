module solver

implicit none

contains

subroutine jacobi1d(coef, vector, b, nx, residual, itermax, criteria)

integer :: size, i, j, itermax, nx
real, allocatable :: coef(:,:), vector(:), residual(:,:), b(:)
real :: current, criteria

size = nx+2

do i = 1,itermax
  do j = 1,size
    current = (b(j) + coef(j, 3)*vector(j-1) + coef(j,4)*vector(j+1))/coef(j,1)
    residual(i,j) = abs(current - coef(j,2))
    coef(j,2) = current
    end do
  vector(:) = coef(:,2)
  residual(i,size+1) = sum(residual(i,1:size))
  if (criteria .gt. residual(i, size+1)) then
  write(*,*) 'Convergency Achieved, Iteration ', i
  residual(:,:) = residual(1:i, :)
  exit
  end if
end do

end subroutine jacobi1d


subroutine solver_transiente(coef, vector, b, nx, residual, ts, itermax, criteria, f)

integer :: size, i, j, itermax, nx, k, ts
real, allocatable :: coef(:,: ,:), vector(:,:), residual(:, :,:), b(:, :)
real :: it, criteria, f

size = nx+2


do i = 1,itermax
  do j = 1,size
    it = b(ts, j) + coef(ts, j, 3)*(vector(ts, j-1)*f + (1-f)*vector(ts-1, j-1))
    it = (it + coef(ts, j,4)*(f*vector(ts, j+1) + (1-f)*vector(ts-1, j+1)))/coef(ts, j,1)
    residual(i,ts,j) = abs(it - coef(ts, j,2))
    coef(ts,j,2) = it
    end do
  vector(ts, :) = coef(ts, :,2)
  residual(i,ts, size+1) = sum(residual(i,ts, 1:size))

  if (criteria .gt. residual(i, ts, size+1)) then
  write(*,*) 'Convergency Achieved'
  exit

  end if
end do

end subroutine solver_transiente


end module solver
