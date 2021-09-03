program main

use mesh
use solver
implicit none


integer :: i, j, count, k
character(len = 30) :: FMT

!------------------------------------------------------------------------------!
!----------------------------Domain & Mesh Settings----------------------------!
!------------------------------------------------------------------------------!
real :: Lx, Ly, Lt
integer :: Nx, Ny, vector_dim, nt, size
real :: stepx, stepy, stept
real, allocatable :: vectors(:,:,:), nodes(:,:)

!------------------------------------------------------------------------------!
!----------------------------Temperature Coefficient---------------------------!
!------------------------------------------------------------------------------!

real::ae, as, aw, an, ap, ap0
real, allocatable :: Temperature(:,:), eq_temp(:,:,:), b(:,:), residual(:,:,:)
real:: alpha, cp, Gamma, rho, q_gen, criteria, time_factor, t0
real, dimension(4,2) :: T_cond
integer :: itermax

!------------------------------------------------------------------------------!
!------------------------------Domain & Mesh Call------------------------------!
!------------------------------------------------------------------------------!


Lx = 1.0        !Mesh's X-size
Ly = 1.0        !Mesh's Y-size

Nx = 25       !Mesh's X-elements
Ny = 25         !Mesh's Y-elements

stept = 0.002

Nt = 6000

Lt = (nt)*stept

stepx = Lx/Nx
stepy = Ly/Ny

size = (nx+2)*(ny+2)

allocate(nodes(size, 2))
vector_dim = max0((Nx+2)*(Ny+1), (Nx+1)*(Ny+2))
allocate(vectors(vector_dim, 2, 2))

call mesher2d(Lx, Ly, Nx, Ny, vectors, nodes)


!------------------------------------------------------------------------------!
!-------------------------------Temperature Call-------------------------------!
!------------------------------------------------------------------------------!

alpha = 5    !W/mK
cp = 130      !J/kgK
rho = 5       !kg/m^2
q_gen = 0     !W
t0 = 0
time_factor = 1     !1 -> implicit, 0 -> explicit, 0.5 ->Crank-Nielson

itermax = 20
criteria = 1e-8

allocate(eq_temp(nt+1, size, 6), temperature(nt+1, size), b(nt+1, size))
allocate(residual(itermax, nt+1, size+1))


temperature = t0 !Initial Temperature
T_cond(1,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Western side
T_cond(2,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Eastern side
T_cond(3,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Northern side
T_cond(4,1) = 1 !1 -> Fixed Temperature, 2-> Fixed Heat Transfer, Southern side

T_cond(1,2) = 0 !Heat transfer or Fixed Temperature on Western side
T_cond(2,2) = 50 !Heat transfer or Fixed Temperature on Eastern side
T_cond(3,2) = 0 !Heat transfer or Fixed Temperature on Northern side
T_cond(4,2) = 0 !Heat transfer or Fixed Temperature on Southern side



!Southern Border
do i = 0, nx+1
  eq_temp(:,i*(Nx+2) + 1, 1) = 1
  eq_temp(:,i*(Nx+2) + 1, 3) = 0
  eq_temp(:,i*(Nx+2) + 1, 4) = 0
  eq_temp(:,i*(Nx+2) + 1, 5) = 0

  if (T_cond(4,1) == 1) then
    eq_temp(:,i*(Nx+2) + 1, 6) = 0
    temperature(:,i*(Nx+2) + 1) = T_cond(4,2)
    eq_temp(:,i*(Nx+2) + 1, 2) = T_cond(4,2)
    b(:,i*(Nx+2) + 1) = T_cond(4,2)

  else
    eq_temp(:,i*(Nx+2) + 1, 6) = 1
    b(:,i*(Nx+2) + 1) = -T_cond(4,2)*(nodes(i*(Nx+2) + 2,2) - nodes(i*(Nx+2) + 1,2))/alpha
    eq_temp(:,i*(Nx+2) + 1, 2) = temperature(:,i*(Nx+2)+1)
    end if
  end do


!Northern Border
do i = 1, nx+2
  eq_temp(:,i*(Nx+2), 1) = 1
  eq_temp(:,i*(Nx+2), 3) = 0
  eq_temp(:,i*(Nx+2), 4) = 0
  eq_temp(:,i*(Nx+2), 6) = 0

  if (T_cond(3,1) == 1) then
    eq_temp(:,i*(Nx+2), 5) = 0
    temperature(:,i*(Nx+2)) = T_cond(3,2)
    eq_temp(:,i*(Nx+2), 2) = T_cond(3,2)
    b(:,i*(Nx+2)) = T_cond(3,2)

  else
    eq_temp(:,i*(Nx+2), 5) = 1
    eq_temp(:,i*(Nx+2), 2) = temperature(:,i*(nx+2))
    b(:,i*(Nx+2)) = T_cond(3,2)*(nodes(i*(Nx+2),2) - nodes(i*(Nx+2) - 1,2))/alpha
    end if
end do


!Western Border
  do i = 1, ny+2
    eq_temp(:,i, 1) = 1
    eq_temp(:,i, 4) = 0
    eq_temp(:,i, 5) = 0
    eq_temp(:,i, 6) = 0

    if (T_cond(1,1) == 1) then
      eq_temp(:,i, 3) = 0
      temperature(:,i) = T_cond(1,2)
      eq_temp(:,i, 2) = T_cond(1,2)
      b(:,i) = T_cond(1,2)

    else
      eq_temp(:,i, 3) = 1
      eq_temp(:,i, 2) = temperature(:,i)
      b(:,i) = -T_cond(1,2)*(nodes(i,1) - nodes(i - (ny+2),1))/alpha
      end if
  end do

!Eastern Border
  do i = (Nx+1)*(Ny+2) + 1, (Nx+2)*(Ny+2)
    eq_temp(:,i, 1) = 1
    eq_temp(:,i, 3) = 0
    eq_temp(:,i, 5) = 0
    eq_temp(:,i, 6) = 0

    if (T_cond(2,1) == 1) then
      eq_temp(:,i, 4) = 0
      temperature(:,i) = T_cond(2,2)
      eq_temp(:,i, 2) = T_cond(2,2)
      b(:,i) = T_cond(2,2)

    else
      eq_temp(:,i, 4) = 1
      eq_temp(:,i, 2) = temperature(:,i)
      b(:,i) = T_cond(2,2)*(nodes(i,1) - nodes(i - (ny+2),1))/alpha
      end if
  end do


do k = 2, nt+1

do i = 2, nx+1
  do j = 2, ny + 1
  count =  j + (Ny+2)*(i-1)

  Gamma = alpha/cp
  aw = gamma/(nodes(count,1)-nodes(count - Ny - 2,1))*stepy
  ae = gamma/(nodes(count+ ny + 2,1)-nodes(count,1))*stepy
  as = gamma/(nodes(count+1, 2) - nodes(count, 2))*stepx
  an = gamma/(nodes(count, 2) - nodes(count - 1, 2))*stepx
  ap = rho*stepx*stepy/stept + time_factor*(ae + aw + an + as)
  ap0 = rho*stepx*stepy/stept

  b(k, count) = (ap0 - (1-time_factor)*(ae+aw+as+an))*temperature(k-1, count) + q_gen*stepx*stepy
  eq_temp(k, count, 1) = ap
  eq_temp(k, count, 2) = temperature(k, count)
  eq_temp(k, count, 3) = aw
  eq_temp(k, count, 4) = ae
  eq_temp(k, count, 5) = as
  eq_temp(k, count, 6) = an

  end do
end do



call solver_transiente(eq_temp, temperature, b, nx, ny, residual, k, itermax, criteria, time_factor)

end do

open(11, file = 'solution-transient.txt', status = 'replace')
write(11, *) '2-D Temperature Solution'
write(11, *) 'Lx = ', Lx, ' ; Ly = ', Ly, ' ;Lt = ', nt*stept
write(11, *) 'Nx = ', Nx, ' ; Ny = ', Ny, ' ;Nt = ', nt
write(11, *) 'Conditions'
write(11, *) T_cond(1,:)
write(11, *) T_cond(2,:)
write(11, *) T_cond(3,:)
write(11, *) T_cond(4,:)
write(11, *) T0
write(11, *) 'Q_gen = ', Q_gen
write(11, *) '*---------------------------------*'

do k = 0, nt
  FMT = "(A, F11.8)"
  write(11, FMT) 'Time = ' , k*stept
  do i = 1, nx+2
    do j = 1, ny+2
      FMT = "(A, F11.6, A, F11.6, A, F18.6)"
      count = j + (Ny+2)*(i-1)
      write(11, FMT) 'x = ', nodes(Count,1), ' ;y = ', nodes(count,2), ' ;T = ', temperature(k+1, count)
    end do
  end do
end do

end program main
