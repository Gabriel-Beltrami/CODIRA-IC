program main

use mesher
use solver

implicit none


integer :: i, j, k
character(len = 30) :: FMT


!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx
integer :: Nx,vector_dim
real :: stepx
real, allocatable :: vectors(:,:,:), nodes(:,:)

!------------------------------------------------------------------------------!
!----------------------------Temperature Coefficient---------------------------!
!------------------------------------------------------------------------------!

real::ae, as, aw, an, ap
real, allocatable :: Temperature(:), eq_temp(:,:), b(:), residual(:,:)
real:: alpha, cp, Gamma, q_gen, criteria
integer :: itermax
real, dimension(2,2) :: T_cond

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!

Lx = 1.0        !Mesh's X-size

Nx = 5         !Mesh's X-elements

stepx = Lx/Nx

allocate (nodes(Nx+2, 1))
vector_dim = Nx+1
allocate (vectors(vector_dim, 1, 1))

call mesher1d(Lx, Nx, vectors, nodes)




!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 35
cp = 130
itermax = 500
q_gen = 0
criteria = 1e-10

allocate(eq_temp(nx+2, 4), temperature(nx+2), b(nx+2))
allocate(residual(itermax,nx+3))


temperature = 50. !Initial Temperatures
T_cond(1,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Left side
T_cond(2,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Right side

T_cond(1,2) = 0 !Heat transfer or Fixed Temperature on Left side
T_cond(2,2) = 0 !Heat transfer or Fixed Temperature on Right side




!Left Border Condition

eq_temp(1,1) = 1
eq_temp(1,3) = 0

if (T_cond(1,1) == 1) then
  b(1) = T_cond(1,2)
  eq_temp(1,4) = 0
  eq_temp(1,2) = T_cond(1,2)
  temperature(1) = T_cond(1,2)
else if (T_cond(1,1) == 2) then
  eq_temp(1,4) = 1
  b(1) = T_cond(1,2)*(nodes(i+1,1)-nodes(i,1))/alpha
  eq_temp(1,2) = temperature(1)
end if

!Right Border Condition

eq_temp(nx+2,1) = 1
eq_temp(nx+2,4) = 0


if (T_cond(2,1) == 1) then
  b(nx+2) = T_cond(2,2)
  eq_temp(nx+2,3) = 0
  temperature(nx+2) = T_cond(2,2)
  eq_temp(nx+2,2) = T_cond(2,2)
else if (T_cond(2,1) == 2) then
  eq_temp(nx+2,3) = 1
  b(nx+2) = T_cond(2,2)*(nodes(nx+2,1)-nodes(nx+1,1))/alpha
  eq_temp(nx+2,2) = temperature(nx+2)
end if



do i = 2, nx+1
  Gamma = alpha/cp
  aw = gamma/(nodes(i,1)-nodes(i-1,1))
  ae = gamma/(nodes(i+1,1)-nodes(i,1))
  ap = aw + ae
  eq_temp(i, 1) = ap
  eq_temp(i, 2) = temperature(i)
  eq_temp(i, 3) = aw
  eq_temp(i, 4) = ae
  b(i) = q_gen*stepx
end do

do i = 1, 7
write(*,*) eq_temp(i,:)
end do

call jacobi1d(eq_temp, temperature, b, nx, residual, itermax, criteria)
write(*,*) temperature
FMT = "(A, F8.6, A, F12.6)"

open(11, file = 'solution-steady.txt', status = 'replace')
write(11, *) '1-D Temperature Solution'
write(11, *) 'Lx = ', Lx
write(11, *) 'Nx = ', Nx
write(11, *) 'Conditions'
write(11, *) T_cond(1,:)
write(11, *) T_cond(2,:)
write(11, *) 'Q_gen = ', Q_gen
write(11, *) '*---------------------------------*'
do i = 1, nx+2
  write(11, FMT) 'x = ', nodes(i,1), ' ;T = ', temperature(i)
end do

end program main
