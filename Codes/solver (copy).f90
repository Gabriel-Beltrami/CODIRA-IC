module solver

implicit none

contains

subroutine jacobi(matrix, b, value, residual, itermax)

integer :: i, j, itermax, ss
real, allocatable :: matrix (:,:), b(:), residual(:,:), value(:), iteration(:)
real :: current, res


ss = size(matrix(1,:))

allocate(iteration(ss))

do i = 1, itermax
  do j = 1, ss
    if (j == 1) then
      current = sum(matrix(2:ss, j) * value(2:ss))
    else if (j == ss) then
      current = sum(matrix(1:ss-1, j)*value(1:ss-1))
    else
      current = sum(matrix(1:j-1, j) * value(1:j-1)) + sum(matrix(j+1:ss, j) * value(j+1:ss))
    end if
    iteration(j) = (b(j) - current)/matrix(j,j)
    residual(i, j) = abs(iteration(j) - value(j))
  end do
  value = iteration
  residual(i, ss+1) = sum(residual(i, 1:ss))
  !write(*,*) value
  !write(*,*) residual(i, ss+1)
  end do
write(*,*) value

return
end subroutine jacobi
end module solver
