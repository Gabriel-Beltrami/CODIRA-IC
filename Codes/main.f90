program main

use mesh
use solver

implicit none


integer :: i, j, k

!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Ly, Lz
integer :: Nx, Ny, Nz, dimension, vector_dim
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

dimension = 2   !Mesh Dimension (1, 2 or 3D)

Lx = 1.0        !Mesh's X-size
Ly = 1.0        !Mesh's Y-size
Lz = 1.0        !Mesh's Z-size

Nx = 200         !Mesh's X-elements
Ny = 200          !Mesh's Y-elements
Nz = 3          !Mesh's Z-elements

stepx = Lx/Nx
stepy = Ly/Ny
stepz = Lz/Nz

if (dimension == 1) then
  allocate (nodes(Nx+2, dimension))
  vector_dim = Nx+1
  allocate (vectors(vector_dim, dimension, dimension))

else if (dimension == 2) then
  allocate(nodes((Nx+2) * (ny+2), dimension))
  vector_dim = max0((Nx+2)*(Ny+1), (Nx+1)*(Ny+2))
  allocate(vectors(vector_dim, dimension, dimension))

else if (dimension == 3) then
allocate(nodes((Nx+2) * (ny+2) * + (nz + 2), dimension))
vector_dim=max0((Nx+2)*(Ny+1)*(Nz+2),(Nx+1)*(Ny+2)*(Nz+2),(Nz+1)*(Nx+2)*(Ny+2))
allocate(vectors(vector_dim, dimension, dimension))

end if


call mesher(Lx, Ly, Lz, Nx, Ny, Nz, dimension, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 35
cp = 130
itermax = 100

if (dimension == 1) then
  allocate(systemT(nx+2, nx+2), temperature(nx+2), b(nx+2))
  allocate(residual(itermax,nx+3))

  !Initial Conditions
  temperature = 50.
  TW = 0
  TE = 100

  !Border Conditions
  systemT(1,1) = 1
  temperature(1) = TW
  b(1) = TW

  systemT(nx+2,nx+2) = 1
  b(nx+2) = Te
  temperature(nx+2) = Te


  !System construction
  do i = 2, nx+1
    Gamma = alpha/cp
    aw = gamma/(nodes(i,1)-nodes(i-1,1))
    ae = gamma/(nodes(i+1,1)-nodes(i,1))
    ap = aw + ae
    b(i) = 0
    systemT(i, i) = ap
    systemT(i-1, i) = -aw
    systemT(i+1, i) = -ae
  end do



else if (dimension == 2) then

  allocate(systemT((Nx+2) * (ny+2), (Nx+2) * (ny+2)))
  allocate(temperature((Nx+2) * (ny+2)), b((Nx+2) * (ny+2)))
  allocate(residual(itermax,(Nx+2) * (ny+2) + 1))

  !Initial Conditions
  temperature = 50.
  TW = 0
  TE = 0
  TN = 100
  TS = 100

  !Border Conditions
  do i = 1, ny+2
    systemT(i,i) = 1
    temperature(i) = TW
    b(i) = TW
    end do

  do i = (Nx+1)*(Ny+2),(Nx+2)*(Ny+2)
    systemT(i,i) = 1
    temperature(i) = TE
    b(i) = TE
    end do

  do i = 1, ny+2
    systemT((i-1)*(Nx+2) + 1, (i-1)*(Nx+2) + 1) = 1
    temperature((i-1)*(Nx+2) + 1) = TS
    b((i-1)*(Nx+2) + 1) = TS
    end do

  do i = 1, ny+2
    systemT(i*(Nx+2), i*(Nx+2)) = 1
    temperature(i*(Nx+2)) = TN
    b(i*(Nx+2)) = TN
    end do


  !System construction
  do i = 2, nx+1
    do j = 2, ny+1
      tracker = (i-1)*(Ny+2) + j
      Gamma = alpha/cp

      aw = gamma/(nodes(tracker,1)-nodes(tracker-(Ny+2),1)) * stepy
      ae = gamma/(nodes(tracker+(Ny+2),1)-nodes(tracker,1)) * stepy
      as = gamma/(nodes(tracker,2)-nodes(tracker-1,2)) * stepx
      an = gamma/(nodes(tracker+1,2) - nodes(tracker,2)) * stepx
      ap = aw + ae + as + an

      b(i) = 0
      systemT(tracker, tracker) = ap
      systemT(tracker-(Ny+2), tracker) = -aw
      systemT(tracker+(Ny+2), tracker) = -ae
      systemT(tracker + 1, tracker) = -an
      systemT(tracker - 1, tracker) = -as
      end do
  end do
end if

call jacobi(systemT, b, temperature, residual, itermax)





end program main
