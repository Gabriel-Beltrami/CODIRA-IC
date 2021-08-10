program main

use mesh
use solver

implicit none


integer :: i, j, count

!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Ly
integer :: Nx, Ny, vector_dim
real :: stepx, stepy
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
Ly = 1.0        !Mesh's Y-size

Nx = 3       !Mesh's X-elements
Ny = 3         !Mesh's Y-elements

stepx = Lx/Nx
stepy = Ly/Ny


allocate(nodes((Nx+2) * (ny+2), 2))
vector_dim = max0((Nx+2)*(Ny+1), (Nx+1)*(Ny+2))
allocate(vectors(vector_dim, 2, 2))

call mesher(Lx, Ly, Nx, Ny, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 35
cp = 130
itermax = 500

allocate(eq_temp((nx+2)*(ny+2), 6), temperature((nx+2)*(ny+2)), b((nx+2)*(ny+2)))
allocate(residual(itermax,(nx+2)*(ny+2)+1))


b = 0.
temperature = 50.
TW = 0
TE = 0
TN = 100
TS = 0


!Southern Border
do i = 0, nx+1
  eq_temp(i*(Nx+2) + 1, 1) = 1
  eq_temp(i*(Nx+2) + 1, 2) = TS
  eq_temp(i*(Nx+2) + 1, 3) = 0
  eq_temp(i*(Nx+2) + 1, 4) = 0
  eq_temp(i*(Nx+2) + 1, 5) = 0
  eq_temp(i*(Nx+2) + 1, 6) = 0
  b(i*(Nx+2) + 1) = TS
  temperature(i*(Nx+2) + 1) = TS

  end do

!Northern Border
do i = 1, nx+2
  eq_temp(i*(Nx+2), 1) = 1
  eq_temp(i*(Nx+2), 2) = TN
  eq_temp(i*(Nx+2), 3) = 0
  eq_temp(i*(Nx+2), 4) = 0
  eq_temp(i*(Nx+2), 5) = 0
  eq_temp(i*(Nx+2), 6) = 0
  b(i*(Nx+2)) = TN
  temperature(i*(Nx+2)) = TN
  end do

!Eastern Border
  do i = 1, ny+2
    eq_temp(i, 1) = 1
    eq_temp(i, 2) = TE
    eq_temp(i, 3) = 0
    eq_temp(i, 4) = 0
    eq_temp(i, 5) = 0
    eq_temp(i, 6) = 0
    b(i) = TE
    temperature(i) = TE
    end do

!Western Border
  do i = (Nx+1)*(Ny+2) + 1, (Nx+2)*(Ny+2)
    eq_temp(i, 1) = 1
    eq_temp(i, 2) = TW
    eq_temp(i, 3) = 0
    eq_temp(i, 4) = 0
    eq_temp(i, 5) = 0
    eq_temp(i, 6) = 0
    b(i) = TW
    temperature(i) = TW
    end do


do i = 2, nx+1
  do j = 2, ny + 1
  count =  j + (Ny+2)*(i-1)

  Gamma = alpha/cp
  aw = gamma/(nodes(count,1)-nodes(count - Ny - 2,1))*stepy
  ae = gamma/(nodes(count+ ny + 2,1)-nodes(count,1))*stepy
  as = gamma/(nodes(count+1, 2) - nodes(count, 2))*stepx
  an = gamma/(nodes(count, 2) - nodes(count - 1, 2))*stepx
  ap = aw + ae + as + an

  eq_temp(count, 1) = ap
  eq_temp(count, 2) = temperature(count)
  eq_temp(count, 3) = aw
  eq_temp(count, 4) = ae
  eq_temp(count, 5) = as
  eq_temp(count, 6) = an

  end do
end do



call jacobi2d(eq_temp, temperature, b, nx, ny, residual, itermax)

open(11, file = 'results.csv', status = 'replace')
write(11,*) 'x coord, y coord, temperature'


do i = 1, nx+2
  do j = 1, ny+2
      count =  j + (Ny+2)*(i-1)
      write(11, *)  nodes(count, 1), ',', nodes(count, 2), ',', temperature(count)
end do
end do




end program main
