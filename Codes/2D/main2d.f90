program main

use mesh
use solver

implicit none


integer :: i, j, count
character(len = 30) :: FMT

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
real:: alpha, cp, Gamma, criteria, q_gen
real, dimension(4,2) :: T_cond
integer :: itermax

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!


Lx = 1.0        !Mesh's X-size
Ly = 1.0        !Mesh's Y-size

Nx = 100      !Mesh's X-elements
Ny = 100       !Mesh's Y-elements

stepx = Lx/Nx
stepy = Ly/Ny


allocate(nodes((Nx+2) * (ny+2), 2))
vector_dim = max0((Nx+2)*(Ny+1), (Nx+1)*(Ny+2))
allocate(vectors(vector_dim, 2, 2))

call mesher2d(Lx, Ly, Nx, Ny, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 35
cp = 130
itermax = 200000
criteria = 1e-10
q_gen = 0     !W


allocate(eq_temp((nx+2)*(ny+2), 6), temperature((nx+2)*(ny+2)), b((nx+2)*(ny+2)))
allocate(residual(itermax,(nx+2)*(ny+2)+1))


temperature = 50. !Initial Temperature
T_cond(1,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Western side
T_cond(2,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Eastern side
T_cond(3,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Northern side
T_cond(4,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Southern side

T_cond(1,2) = 100 !Heat transfer or Fixed Temperature on Western side
T_cond(2,2) = 0 !Heat transfer or Fixed Temperature on Eastern side
T_cond(3,2) = 100 !Heat transfer or Fixed Temperature on Northern side
T_cond(4,2) = 0 !Heat transfer or Fixed Temperature on Southern side

!Southern Border
do i = 0, nx+1
  eq_temp(i*(Ny+2) + 1, 1) = 1
  eq_temp(i*(Ny+2) + 1, 3) = 0
  eq_temp(i*(Ny+2) + 1, 4) = 0
  eq_temp(i*(Ny+2) + 1, 5) = 0

  if (T_cond(4,1) == 1) then
    eq_temp(i*(Ny+2) + 1, 6) = 0
    temperature(i*(Ny+2) + 1) = T_cond(4,2)
    eq_temp(i*(Ny+2) + 1, 2) = T_cond(4,2)
    b(i*(Ny+2) + 1) = T_cond(4,2)

  else
    eq_temp(i*(Ny+2) + 1, 6) = 1
    b(i*(Ny+2) + 1) = -T_cond(4,2)*(nodes(i*(Ny+2) + 2,2) - nodes(i*(Ny+2) + 1,2))/alpha
    eq_temp(i*(Ny+2) + 1, 2) = temperature(i*(Ny+2)+1)
    end if
  end do


!Northern Border
do i = 1, nx+2
  eq_temp(i*(Ny+2), 1) = 1
  eq_temp(i*(Ny+2), 3) = 0
  eq_temp(i*(Ny+2), 4) = 0
  eq_temp(i*(Ny+2), 6) = 0

  if (T_cond(3,1) == 1) then
    eq_temp(i*(Ny+2), 5) = 0
    temperature(i*(Ny+2)) = T_cond(3,2)
    eq_temp(i*(Ny+2), 2) = T_cond(3,2)
    b(i*(Ny+2)) = T_cond(3,2)

  else
    eq_temp(i*(Ny+2), 5) = 1
    eq_temp(i*(Ny+2), 2) = temperature(i*(ny+2))
    b(i*(Ny+2)) = T_cond(3,2)*(nodes(i*(Ny+2),2) - nodes(i*(Ny+2) - 1,2))/alpha
    end if
end do


!Western Border
  do i = 1, ny+2
    eq_temp(i, 1) = 1
    eq_temp(i, 4) = 0
    eq_temp(i, 5) = 0
    eq_temp(i, 6) = 0

    if (T_cond(1,1) == 1) then
      eq_temp(i, 3) = 0
      temperature(i) = T_cond(1,2)
      eq_temp(i, 2) = T_cond(1,2)
      b(i) = T_cond(1,2)

    else
      eq_temp(i, 3) = 1
      eq_temp(i, 2) = temperature(i)
      b(i) = -T_cond(1,2)*(nodes(i,1) - nodes(i - 1,1))/alpha
      end if
  end do

!Eastern Border
  do i = (Nx+1)*(Ny+2) + 1, (Nx+2)*(Ny+2)
    eq_temp(i, 1) = 1
    eq_temp(i, 3) = 0
    eq_temp(i, 5) = 0
    eq_temp(i, 6) = 0

    if (T_cond(2,1) == 1) then
      eq_temp(i, 4) = 0
      temperature(i) = T_cond(2,2)
      eq_temp(i, 2) = T_cond(2,2)
      b(i) = T_cond(2,2)

    else
      eq_temp(i, 4) = 1
      eq_temp(i, 2) = temperature(i)
      b(i) = T_cond(2,2)*(nodes(i,1) - nodes(i - 1,1))/alpha
      end if
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
  b(count) = q_gen*stepx*stepy

  end do

end do

count = 1
do i = 1, nx+2
!write(*,*) " *****************"
  do j = 1, ny+2
    count = i + (Ny+2)*(j-1)
  !  write(*,*) eq_temp(count, :)
  !  write(*,*) b(count)
    count = count+1
    end do
    end do

call jacobi2d(eq_temp, temperature, b, nx, ny, residual, itermax, criteria)


FMT = "(A, F8.6, A, F8.6, A, F12.6)"

open(11, file = 'solution-steady.txt', status = 'replace')
write(11, *) '2-D Temperature Solution'
write(11, *) 'Lx = ', Lx, ' ; Ly = ', Ly
write(11, *) 'Nx = ', Nx, ' ; Ny = ', Ny
write(11, *) 'Conditions'
write(11, *) T_cond(1,:)
write(11, *) T_cond(2,:)
write(11, *) T_cond(3,:)
write(11, *) T_cond(4,:)
write(11, *) 'Q_gen = ', Q_gen
write(11, *) '*---------------------------------*'
do i = 1, nx+2
  do j = 1, ny+2
    count = j + (Ny+2)*(i-1)
    write(11, FMT) 'x = ', nodes(Count,1), ' ;y = ', nodes(count,2), ' ;T = ', temperature(count)
  end do
end do




end program main
