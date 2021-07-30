subroutine mesher(Lx, Ly, Lz, Nx, Ny, Nz)

implicit none
real, intent(in) :: Lx, Ly, Lz
integer, intent(in) :: Nx, Ny, Nz

real :: stepx, stepy, stepz

real, dimension(nx*ny*nz, 2, 8, 3) :: element

integer :: count, x, y, z, mesh, number_nodes, scalar, vector

!distance between nodes
stepx = Lx/Nx
stepy = Ly/Ny
stepz = Lz/Nz


do x = 1, Nx
  do y = 1, Ny
    do z = 1, Nz

      count = z + (y-1)*(Nz) + (x-1)*(Ny)*(Nz)

      !Node (Scalar) Setting

      number_nodes = 1

      element(count, 1, number_nodes, 1) = (x-0.5)*stepx
      element(count, 1, number_nodes, 2) = (y-0.5)*stepy
      element(count, 1, number_nodes, 3) = (z-0.5)*stepz

      !Triple Border

      if (x == 1 .and.  y == 1 .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = 0

      else if (x == Nx .and.  y == 1 .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = 0

      else if (x == 1 .and.  y == Ny .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = 0

      else if (x == 1 .and.  y == 1 .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = Lz

      else if (x == Nx .and.  y == 1 .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = Lz

      else if (x == 1 .and.  y == Ny .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = Lz

      else if (x == Nx .and.  y == Ny .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = 0

      else if (x == Nx .and.  y == Ny .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = Lz
      end if



      !Double Border

      if (x == 1 .and. y == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (x == Nx .and. y == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (x == 1 .and. y == Ny) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (x == Nx .and. y == Ny) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (x == 1 .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 3) = 0
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
      end if

      if (x == Nx .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 3) = 0
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
      end if

      if (x == 1 .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 3) = Lz
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
      end if

      if (x == Nx .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 3) = Lz
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
      end if

      if (y == 1 .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = 0
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
      end if

      if (y == Ny .and. z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = 0
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
      end if

      if (y == 1 .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 3) = Lz
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
      end if

      if (y == Ny .and. z == Nz) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 3) = Lz
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
      end if

      !Single Border

      if (x == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = 0
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (x == Nx) then
      number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lx
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (y == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 2) = 0
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (y == Ny) then
      number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 2) = Ly
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
        element(count, 1, number_nodes, 3) = element(count, 1, 1, 3)
      end if

      if (z == 1) then
        number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 3) = 0
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
      end if

      if (z == Nz) then
      number_nodes = number_nodes + 1
        element(count, 1, number_nodes, 1) = Lz
        element(count, 1, number_nodes, 2) = element(count, 1, 1, 2)
        element(count, 1, number_nodes, 1) = element(count, 1, 1, 1)
      end if

      !Vector Setting

      element(count, 2, 1, 1) = (x-1)*stepx
      element(count, 2, 1, 2) = (y-0.5)*stepy
      element(count, 2, 1, 3) = (z-0.5)*stepz

      element(count, 2, 2, 1) = x*stepx
      element(count, 2, 2, 2) = (y-0.5)*stepy
      element(count, 2, 2, 3) = (z-0.5)*stepz

      element(count, 2, 3, 1) = (x-0.5)*stepx
      element(count, 2, 3, 2) = (y-1)*stepy
      element(count, 2, 3, 3) = (z-0.5)*stepz

      element(count, 2, 4, 1) = (x-0.5)*stepx
      element(count, 2, 4, 2) = y*stepy
      element(count, 2, 4, 3) = (z-0.5)*stepz

      element(count, 2, 5, 1) = (x-0.5)*stepx
      element(count, 2, 5, 2) = (y-0.5)*stepy
      element(count, 2, 5, 3) = (z-1)*stepz

      element(count, 2, 6, 1) = (x-0.5)*stepx
      element(count, 2, 6, 2) = (y-0.5)*stepy
      element(count, 2, 6, 3) = z*stepz


    end do
  end do
end do



open(10, file = 'mesh.txt', status = 'replace')
write(10, *) 'Lx = ', Lx, '; Ly = ', Ly, '; Lz = ', Lz
write(10, *) 'Nx = ', Nx, '; Ny = ', Ny, '; Nz = ', Nz
do mesh = 1, count
  write(10, *) "Element", mesh, "; Center Coordinates: ", element(mesh,1,1,:)
  write(10,*) " "

  do scalar = 1,8
    if (element(mesh, 1, scalar, 1) == 0 .and. mesh /= 1) then
      if (element(mesh, 1, scalar, 2) == 0 .and. element(mesh,1,scalar,3) == 0) then
      exit
    end if
    end if

    write(10,*) "Scalar ", scalar, "; Coordinates: ", element(mesh, 1, scalar, :)

    end do

  write(10,*) " "

  do vector = 1, 6
    write(10,*) "Vector ", vector, "; Coordinates: ", element(mesh, 2, vector, :)
    end do
  write(10,*) " "
  write(10,*) "----------------------------------------------------------------"

end do
end subroutine mesher
