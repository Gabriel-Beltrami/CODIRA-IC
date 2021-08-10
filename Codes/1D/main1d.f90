program main

use mesher
use solver

implicit none


integer :: i, j, k

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
real:: alpha, cp, Gamma
real:: Tw, Te, Tn, Ts, Tb, Tt  !Initial Conditions
integer :: itermax, tracker

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
itermax = 20

allocate(eq_temp(nx+2, 4), temperature(nx+2), b(nx+2))
allocate(residual(itermax,nx+3))


b = 0.
temperature = 50.
TW = 0
TE = 100

eq_temp(1,1) = 1
eq_temp(1,2) = TW
eq_temp(1,3) = 0
eq_temp(1,4) = 0
temperature(1) = Tw
b(1) = TW

eq_temp(nx+2,1) = 1
eq_temp(nx+2,2) = Te
eq_temp(nx+2,3) = 0
eq_temp(nx+2,4) = 0
temperature(nx+2) = TE
b(nx+2) = Te


do i = 2, nx+1
  Gamma = alpha/cp
  aw = gamma/(nodes(i,1)-nodes(i-1,1))
  ae = gamma/(nodes(i+1,1)-nodes(i,1))
  ap = aw + ae
  eq_temp(i, 1) = ap
  eq_temp(i, 2) = temperature(i)
  eq_temp(i, 3) = aw
  eq_temp(i, 4) = ae
end do

call jacobi1d(eq_temp, temperature, b, nx, residual, itermax)





end program main
