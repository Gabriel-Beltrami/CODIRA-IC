program main

use mesh
use solver

implicit none


integer :: i, j, k

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

real::ae, as, aw, an, ap
real, allocatable :: Temperature(:), systemT(:,:), b(:), residual(:,:)
real:: alpha, cp, Gamma
real:: Tw, Te, Tn, Ts, Tb, Tt  !Initial Conditions
integer :: itermax, tracker

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!


Lx = 1.0        !Mesh's X-size
Ly = 1.0        !Mesh's Y-size
Lz = 1.0        !Mesh's Z-size

Nx = 200         !Mesh's X-elements
Ny = 200          !Mesh's Y-elements
Nz = 3          !Mesh's Z-elements

stepx = Lx/Nx
stepy = Ly/Ny
stepz = Lz/Nz

allocate(nodes((Nx+2) * (ny+2) * + (nz + 2), 3))
vector_dim=max0((Nx+2)*(Ny+1)*(Nz+2),(Nx+1)*(Ny+2)*(Nz+2),(Nz+1)*(Nx+2)*(Ny+2))
allocate(vectors(vector_dim, 3, 3))

end if


call mesher(Lx, Ly, Lz, Nx, Ny, Nz, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 35
cp = 130
itermax = 100







end program main
