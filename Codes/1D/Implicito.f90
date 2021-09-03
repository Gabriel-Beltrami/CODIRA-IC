program main

use mesher
use solver

implicit none


integer :: i, j, k
character(len = 30) :: FMT

!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Lt
integer :: Nx, vector_dim, Nt
real :: stepx, stept
real, allocatable :: vectors(:,:,:), nodes(:,:)


!------------------------------------------------------------------------------!
!----------------------------Temperature Coefficient---------------------------!
!------------------------------------------------------------------------------!

real::ae, aw, ap, ap0
real, allocatable :: Temperature(:,:), eq_temp(:,:,:), b(:,:), residual(:,:,:)
real:: alpha, cp, Gamma, rho, q_gen, criteria, time_factor, T0
real, dimension(2,2) :: T_cond
integer :: itermax

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!

Lx = 1.0        !Mesh's X-size

Nx = 5         !Mesh's X-elements

stept = 0.0001

Nt = 1000

Lt = (nt)*stept

stepx = Lx/Nx

allocate (nodes(Nx+2, 1))
vector_dim = Nx+1
allocate (vectors(vector_dim, 1, 1))

call mesher1d(Lx, Nx, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 2    !W/mK
cp = 1      !J/kgK
rho = 1       !kg/m
q_gen = 0   !W
T0 = 20
time_factor = 1     !1 -> implicit, 0 -> explicit, 0.5 ->Crank-Nielson

itermax = 2000
criteria = 1e-8

allocate(eq_temp(nt+1, nx+2, 4), temperature(nt + 1, nx+2), b(nt+1, nx+2))
allocate(residual(itermax, nt + 1, nx+3))

!Initial Conditions
temperature(1,:) = T0 !Initial Temperature
do i = 2, nx+1
  temperature(1,i) = T0*sin(nodes(i,1)*2.D0*DASIN(1.D0))
end do

T_cond(1,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Left side
T_cond(2,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Right side

T_cond(1,2) = 0 !Heat transfer or Fixed Temperature on Left side
T_cond(2,2) = 0 !Heat transfer or Fixed Temperature on Right side


temperature(:, 1) = T_cond(1,2)
eq_temp(:, 1,3) = 0
eq_temp(:, 1,1) = 1
eq_temp(:, 1,2) = T_cond(1,2)

if (T_cond(1,1) == 1) then
  eq_temp(:, 1,4) = 0
  b(:, 1) = T_cond(1,2)
else
  eq_temp(:,1,4) = 1
  b(:,1) = T_cond(1,2)*(nodes(i+1,1)-nodes(i,1))/alpha
end if


eq_temp(:, nx+2,1) = 1
eq_temp(:, nx+2,2) = T_cond(2,2)
temperature(:, nx+2) = T_cond(2,2)
eq_temp(:, nx+2,4) = 0

if (T_cond(2,1) == 1) then
  eq_temp(:, nx+2,3) = 0
  b(:, nx+2) = T_cond(2,2)
else
  eq_temp(:,nx+2,3) = 1
  b(:,nx+2) = T_cond(2,2)*(nodes(nx+2,1)-nodes(nx+1,1))/alpha
end if



do i = 2, nt+1

do j = 2, nx+1
  Gamma = alpha/cp
  aw = gamma/(nodes(j,1)-nodes(j-1,1))
  ae = gamma/(nodes(j+1,1)-nodes(j,1))
  ap = rho*stepx/stept + time_factor*(ae + aw)
  ap0 = rho*stepx/stept

  b(i,j) = (ap0 - (1-time_factor)*(ae+aw))*temperature(i-1, j) + q_gen*stepx
  eq_temp(i, j, 1) = ap
  eq_temp(i, j, 2) = temperature(i, j)
  eq_temp(i, j, 3) = aw
  eq_temp(i, j, 4) = ae

end do

call solver_transiente(eq_temp, temperature, b, nx, residual, i, itermax, criteria, time_factor)

end do

open(11, file = 'solution-transient.txt', status = 'replace')
write(11, *) '1-D Temperature Solution'
write(11, *) 'Lx = ', Lx, ' ;Lt = ', nt*stept
write(11, *) 'Nx = ', Nx, ' ;Nt = ', nt
write(11, *) 'Conditions'
write(11, *) T_cond(1,:)
write(11, *) T_cond(2,:)
write(11, *) T0
write(11, *) 'Q_gen = ', Q_gen
write(11, *) '*---------------------------------*'


do j = 1, nt+1
FMT = "(A, F8.6)"
write(11, FMT) 'Time = ', (j-1)*stept
FMT = "(A, F11.6, A, F11.6)"
  do i = 1, nx+2
    write(11, FMT) 'x = ', nodes(i,1), ' ;T = ', temperature(j, i)
    end do
  end do



end program main
