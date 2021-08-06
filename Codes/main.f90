program main

use mesh


implicit none



!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Ly, Lz                                            !Domain dimensions
integer :: Nx, Ny, Nz, dimension, vector_dim                !Mesh's number of nodes in each direction
real, allocatable :: vectors(:,:,:), nodes(:,:)

dimension = 1

Lx = 1.0
Ly = 1.0
Lz = 1.0

Nx = 3
Ny = 3
Nz = 3

if (dimension == 1) then
  allocate (nodes(Nx+2, dimension))
  allocate (vectors(Nx+1, dimension, dimension))

else if (dimension == 2) then
  allocate(nodes((Nx+2) * (ny+2), dimension))
  vector_dim = max0((Nx+2)*(Ny+1), (Nx+1)*(Ny+2))
  allocate(vectors(vector_dim, dimension, dimension))

else if (dimension == 3) then
allocate(nodes((Nx+2) * (ny+2) * + (nz + 2), dimension))
vector_dim = max0((Nx+2)*(Ny+1)*(Nz+2), (Nx+1)*(Ny+2)*(Nz+2), (Nz+1)*(Nx+2)*(Ny+2))
allocate(vectors(vector_dim, dimension, dimension))

end if


call mesher(Lx, Ly, Lz, Nx, Ny, Nz, dimension, vectors, nodes)

write(*,*) vectors
write(*,*) Nodes

end program main
