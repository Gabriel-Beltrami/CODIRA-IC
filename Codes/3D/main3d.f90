program main

use mesh
use solver

implicit none


integer :: i, j, k, count
character(len = 30) :: FMT


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
real:: alpha, cp, Gamma, q_gen, criteria, t0
real, dimension(6,2) :: T_cond
integer :: itermax

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!


Lx = 1.0        !Mesh's X-size
Ly = 1.0        !Mesh's Y-size
Lz = 1.0        !Mesh's Z-size

Nx = 3        !Mesh's X-elements
Ny = 3          !Mesh's Y-elements
Nz = 3          !Mesh's Z-elements

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
q_gen = 5
criteria = 1e-10
t0 = 5

allocate(eq_temp((nx+2)*(ny+2)*(nz+2), 8), temperature((nx+2)*(ny+2)*(nz+2)))
allocate(b((nx+2)*(ny+2)*(nz+2)))
allocate(residual(itermax,(nx+2)*(ny+2)*(nz+2)+1))


temperature(:) = t0 !Initial Temperature
T_cond(1,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Western side
T_cond(2,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Eastern side
T_cond(3,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Northern side
T_cond(4,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Southern side
T_cond(5,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Bottom side
T_cond(6,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Top side

T_cond(1,2) = 5 !Heat transfer or Fixed Temperature on Western side
T_cond(2,2) = 20 !Heat transfer or Fixed Temperature on Eastern side
T_cond(3,2) = 5 !Heat transfer or Fixed Temperature on Northern side
T_cond(4,2) = 5 !Heat transfer or Fixed Temperature on Southern side
T_cond(5,2) = 5 !Heat transfer or Fixed Temperature on Bottom side
T_cond(6,2) = 5 !Heat transfer or Fixed Temperature on Top side



!Western Border
do i = 1, (ny+2)*(nz+2)
  eq_temp(i, 1) = 1
  eq_temp(i, 3) = 0
  eq_temp(i, 5) = 0
  eq_temp(i, 6) = 0
  eq_temp(i, 7) = 0
  eq_temp(i, 8) = 0

  if (T_cond(1,1) == 1) then
    eq_temp(i, 4) = 0
    temperature(i) = T_cond(1,2)
    eq_temp(i, 2) = T_cond(1,2)
    b(i) = T_cond(1,2)

  else
    eq_temp(i, 4) = 1
    eq_temp(i, 2) = temperature(i)
    b(i) = T_cond(1,2)*(nodes(i + (ny+2)*(nz+2),1) - nodes(i,1))/alpha
    end if
end do


!Eastern Border
do i = (Nx+1)*(Ny+2)*(nz+2) + 1, (Nx+2)*(Ny+2)*(nz+2)
  eq_temp(i, 1) = 1
  eq_temp(i, 4) = 0
  eq_temp(i, 5) = 0
  eq_temp(i, 6) = 0
  eq_temp(i, 7) = 0
  eq_temp(i, 8) = 0

  if (T_cond(2,1) == 1) then
    eq_temp(i, 3) = 0
    temperature(i) = T_cond(2,2)
    eq_temp(i, 2) = T_cond(2,2)
    b(i) = T_cond(2,2)

  else
    eq_temp(i, 3) = 1
    eq_temp(i, 2) = temperature(i)
    b(i) = -T_cond(2,2)*(nodes(i,1) - nodes(i - (ny+2)*(nz+2),1))/alpha
    end if
end do

!Northern Border
do i = 0, (nx+1)
  do j = 1, (nz+2)
  count = i*(Ny+2)*(Nz+2) + j + (Ny+1)*(Nz+2)
  eq_temp(count, 1) = 1
  eq_temp(count, 3) = 0
  eq_temp(count, 4) = 0
  eq_temp(count, 5) = 0
  eq_temp(count, 8) = 0
  eq_temp(count, 7) = 0
  if (T_cond(3,1) == 1) then
    eq_temp(count,6) = 0
    b(count) = T_cond(3,2)
    temperature(count) = T_cond(3,2)
    eq_temp(count,2) = T_cond(3,2)
  else
    eq_temp(count,6) = 1
    b(count) = T_cond(3,2)*(nodes(count,2) - nodes(count-(Nz+2),2))/alpha
    eq_temp(count, 2) = temperature(count)
  end if
end do
end do


!Southern Border
do i = 0, (nx+1)
  do j = 1, (nz+2)
    count = i*(Ny+2)*(Nz+2) + j
    eq_temp(count, 1) = 1
    eq_temp(count, 3) = 0
    eq_temp(count, 4) = 0
    eq_temp(count, 5) = 0
    eq_temp(count, 8) = 0
    eq_temp(count, 7) = 0

    if (T_cond(4,1) == 1) then
      eq_temp(count,6) = 0
      b(count) = T_cond(4,2)
      temperature(count) = T_cond(4,2)
      eq_temp(count,2) = T_cond(4,2)
    else
      eq_temp(count,6) = 1
      b(count) = -T_cond(4,2)*(nodes(count + Nz+2,2) - nodes(count,2))/alpha
      eq_temp(count, 2) = temperature(count)
    end if
  end do
end do

!Bottom Border
do i = 0, (nx+2)*(ny+2) -1
  eq_temp(i*(nz+2)+1, 1) = 1
  eq_temp(i*(nz+2)+1, 3) = 0
  eq_temp(i*(nz+2)+1, 4) = 0
  eq_temp(i*(nz+2)+1, 5) = 0
  eq_temp(i*(nz+2)+1, 6) = 0
  eq_temp(i*(nz+2)+1, 7) = 0

  if (T_cond(5,1) == 1) then
    eq_temp(i*(nz+2)+1, 8) = 0
    b(i*(nz+2)+1) = T_cond(5,2)
    temperature(i*(nz+2)+1) = T_cond(5,2)
    eq_temp(i*(nz+2)+1,2) = T_cond(5,2)
  else
    eq_temp(i*(nz+2)+1, 8) = 1
    b(i*(nz+2)+1) = -T_cond(5,2)*(nodes(i*(nz+2)+2,3) - nodes(i*(nz+2)+1,3))/alpha
    eq_temp(i*(nz+2)+1,2) = temperature(i*(nz+2)+1)
  end if
end do

!Top Border
do i = 1, (nx+2)*(ny+2)
  eq_temp(i*(nz+2), 1) = 1
  eq_temp(i*(nz+2), 3) = 0
  eq_temp(i*(nz+2), 4) = 0
  eq_temp(i*(nz+2), 5) = 0
  eq_temp(i*(nz+2), 6) = 0
  eq_temp(i*(nz+2), 8) = 0

  if (T_cond(6,1) == 1) then
    eq_temp(i*(nz+2), 7) = 0
    b(i*(nz+2)) = T_cond(6,2)
    temperature(i*(nz+2)) = T_cond(6,2)
    eq_temp(i*(nz+2),2) = T_cond(6,2)

  else
    eq_temp(i*(nz+2), 7) = 1
    b(i*(nz+2)+1) = T_cond(6,2)*(nodes(i*(nz+2),3) - nodes(i*(nz+2)-1,3))/alpha
    eq_temp(i*(nz+2),2) = temperature(i*(nz+2))
  end if

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

  b(count) = q_gen*stepx*stepy*stepz
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


do i = 1, (nx+2)*(ny+2)*(nz+2)
  write(*,*) eq_temp(i,:), b(i)
  end do



call jacobi3d(eq_temp, temperature, b, nx, ny, nz, residual, itermax, criteria)

FMT = "(A, F8.6, A, F8.6, A, F8.6, A, F12.6)"

open(11, file = 'solution-steady.txt', status = 'replace')
write(11, *) '3-D Temperature Solution'
write(11, *) 'Lx = ', Lx, ' ; Ly = ', Ly, ' ;Lz = ', Lz
write(11, *) 'Nx = ', Nx, ' ; Ny = ', Ny, ' ;Nz = ', Nz
write(11, *) 'Conditions'
write(11, *) T_cond(1,:)
write(11, *) T_cond(2,:)
write(11, *) T_cond(3,:)
write(11, *) T_cond(4,:)
write(11, *) T_cond(5,:)
write(11, *) T_cond(6,:)
write(11, *) T0
write(11, *) 'Q_gen = ', Q_gen
write(11, *) '*---------------------------------*'
do k = 1, nx + 2
  do i = 1, ny+2
    do j = 1, nz+2
      count = k + (j-1)*(Nz+2) + (Nz+2)*(Ny+2)*(i-1)
      write(11, FMT) 'x = ', nodes(Count,1), ' ;y = ', nodes(count,2), ' z = ', nodes(count, 3), ' ;T = ', temperature(count)
    end do
  end do
end do
end program main
