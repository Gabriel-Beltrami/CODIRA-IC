program main

implicit none

interface
  subroutine mesher(Lx, Ly, Lz, Nx, Ny, Nz)
  real, intent(in) :: Lx, Ly, Lz
  integer, intent(in) :: Nx, Ny, Nz
end subroutine mesher
end interface

!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Ly, Lz                                            !Domain dimensions
integer :: Nx, Ny, Nz                  !Mesh's number of nodes in each direction

Lx = 1.0
Ly = 1.0
Lz = 1.0

Nx = 5
Ny = 5
Nz = 5

call mesher(Lx, Ly, Lz, Nx, Ny, Nz)


end program main
