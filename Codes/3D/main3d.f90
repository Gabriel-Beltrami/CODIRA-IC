program main

use mesh
use solver

implicit none


integer :: i, j, k, count

!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Ly, Lz
integer :: Nx, Ny, Nz, vector_dim
real :: stepx, stepy, stepz
real, allocatable :: vectors(:,:,:), nodes(:,:)

!------------------------------------------------------------------------------!
!----------------------------Temperature Coefficient---------------------------!
!------------------------------------------------------------------------------!

real::ae, as, aw, an, ap, ab, at
real, allocatable :: Temperature(:), eq_temp(:,:), b(:), residual(:,:)
real:: alpha, cp, Gamma
real:: Tw, Te, Tn, Ts, Tb, Tt  !Initial Conditions
integer :: itermax, tracker

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!


Lx = 1.0        !Mesh's X-size
Ly = 1.0        !Mesh's Y-size
Lz = 1.0        !Mesh's Z-size

Nx = 10        !Mesh's X-elements
Ny = 10          !Mesh's Y-elements
Nz = 10          !Mesh's Z-elements

stepx = Lx/Nx
stepy = Ly/Ny
stepz = Lz/Nz

allocate(nodes((Nx+2) * (ny+2) * + (nz + 2), 3))
vector_dim=max0((Nx+2)*(Ny+1)*(Nz+2),(Nx+1)*(Ny+2)*(Nz+2),(Nz+1)*(Nx+2)*(Ny+2))
allocate(vectors(vector_dim, 3, 3))



call mesher(Lx, Ly, Lz, Nx, Ny, Nz, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 35
cp = 130
itermax = 500

allocate(eq_temp((nx+2)*(ny+2)*(nz+2), 8), temperature((nx+2)*(ny+2)*(nz+2)))
allocate(b((nx+2)*(ny+2)*(nz+2)))
allocate(residual(itermax,(nx+2)*(ny+2)*(nz+2)+1))


b = 0.
temperature = 50.
TW = 50
TE = 50
TN = 100
TS = 50
TB = 100
TT = 50



!Western Border
do i = 1, (ny+2)*(nz+2)
  eq_temp(i, 1) = 1
  eq_temp(i, 2) = TW
  eq_temp(i, 3) = 0
  eq_temp(i, 4) = 0
  eq_temp(i, 5) = 0
  eq_temp(i, 6) = 0
  eq_temp(i, 7) = 0
  eq_temp(i, 8) = 0
  b(i) = TW
  temperature(i) = TW
  end do

!Eastern Border
do i = (Nx+1)*(Ny+2)*(nz+2) + 1, (Nx+2)*(Ny+2)*(nz+2)
  eq_temp(i, 1) = 1
  eq_temp(i, 2) = TE
  eq_temp(i, 3) = 0
  eq_temp(i, 4) = 0
  eq_temp(i, 5) = 0
  eq_temp(i, 6) = 0
  eq_temp(i, 7) = 0
  eq_temp(i, 8) = 0
  b(i) = TE
  temperature(i) = TE
  end do

!Southern Border
do i = 0, (nx+1)
  do j = 1, (nz+2)
  eq_temp(i*(Ny+2)*(Nz+2) + j, 1) = 1
  eq_temp(i*(Ny+2)*(Nz+2) + j, 2) = TS
  eq_temp(i*(Ny+2)*(Nz+2) + j, 3) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j, 4) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j, 5) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j, 6) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j, 8) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j, 7) = 0
  b(i*(Ny+2)*(Nz+2) + j) = TS
  temperature(i*(Ny+2)*(Nz+2) + j) = TS
  end do
end do

!Northern Border
do i = 0, (nx+1)
  do j = 1, (nz+2)
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 1) = 1
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 2) = TN
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 3) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 4) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 5) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 6) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 8) = 0
  eq_temp(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2), 7) = 0
  b(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2)) = TN
  temperature(i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2)) = TN
  end do
end do

!Bottom Border
do i = 0, (nx+2)*(ny+2) -1
  eq_temp(i*(nz+2)+1, 1) = 1
  eq_temp(i*(nz+2)+1, 2) = TB
  eq_temp(i*(nz+2)+1, 3) = 0
  eq_temp(i*(nz+2)+1, 4) = 0
  eq_temp(i*(nz+2)+1, 5) = 0
  eq_temp(i*(nz+2)+1, 6) = 0
  eq_temp(i*(nz+2)+1, 7) = 0
  eq_temp(i*(nz+2)+1, 8) = 0
  b(i*(nz+2)+1) = TB
  temperature(i*(nz+2)+1) = TB
end do

!Top Border
do i = 1, (nx+2)*(ny+2)
  eq_temp(i*(nz+2), 1) = 1
  eq_temp(i*(nz+2), 2) = TT
  eq_temp(i*(nz+2), 3) = 0
  eq_temp(i*(nz+2), 4) = 0
  eq_temp(i*(nz+2), 5) = 0
  eq_temp(i*(nz+2), 6) = 0
  eq_temp(i*(nz+2), 7) = 0
  eq_temp(i*(nz+2), 8) = 0
  b(i*(nz+2)) = TT
  temperature(i*(nz+2)) = TT
end do




do i = 2, nx+1
  do j = 2, ny + 1
    do k = 2, nz + 1
    count = k + (j-1)*(Nz+2) + (Nz+2)*(Ny+2)*(i-1)

  Gamma = alpha/cp
  aw = gamma/(nodes(count,1)-nodes(count - (ny + 2)*(Nz+2),1))*stepy*stepz
  ae = gamma/(nodes(count+ (ny + 2)*(Nz+2),1)-nodes(count,1))*stepy*stepz
  as = gamma/(nodes(count+(Nz+2), 2) - nodes(count, 2))*stepx*stepz
  an = gamma/(nodes(count, 2) - nodes(count - (Nz+2), 2))*stepx*stepz
  ab = gamma/(nodes(count, 3) - nodes(count - 1, 3))*stepx*stepy
  at = gamma/(nodes(count+1, 3) - nodes(count, 3))*stepx*stepy

  ap = aw + ae + as + an + ab + at

  eq_temp(count, 1) = ap
  eq_temp(count, 2) = temperature(count)
  eq_temp(count, 3) = aw
  eq_temp(count, 4) = ae
  eq_temp(count, 5) = as
  eq_temp(count, 6) = an
  eq_temp(count, 7) = ab
  eq_temp(count, 8) = at


    end do
  end do
end do


call jacobi3d(eq_temp, temperature, b, nx, ny, nz, residual, itermax)

open(11, file = 'results.csv', status = 'replace')
write(11,*) 'x coord, y coord, z coord, temperature'


do i = 1, nx+2
  do j = 1, ny+2
    do k = 1, nz+2
      count = k + (j-1)*(Nz+2) + (Nz+2)*(Ny+2)*(i-1)
      write(11, *)  nodes(count, 1), ',', nodes(count, 2), ',', nodes(count, 3), ',', temperature(count)
end do
end do
end do


end program main
