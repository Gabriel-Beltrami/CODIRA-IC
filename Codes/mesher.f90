module mesh

implicit none

contains

subroutine mesher(Lx, Ly, Lz, Nx, Ny, Nz, dimension, vectors, nodes)

implicit none

real :: Lx, Ly, Lz
integer :: Nx, Ny, Nz, dimension
real :: stepx, stepy, stepz
real, allocatable :: nodes (:, :), vectors(:, :, :)
integer :: count, x, y, z, mesh, number_nodes, scalar, vector_dim, i

!distance between nodes
stepx = Lx/Nx
stepy = Ly/Ny
stepz = Lz/Nz


if (dimension == 1) then


  count = 0

  do i = 1, Nx+2
    count = count + 1
    if (i == 1) then
      nodes(count, 1) = 0
      vectors(count, 1, 1) = 0
    else if (i == Nx + 2) then
      nodes(count, 1) = Lx
    else
      nodes(count, 1) = stepx*(i-1.5)
      vectors(count, 1, 1) = stepx*(i-1)
    end if
  end do

  open(10, file = 'mesh.txt', status = 'replace')

  write(10,*) "Mesh"
  write(10,*) "Dimension: ", dimension
  write(10,*) "Lx = ", Lx
  write(10,*) "Nx = ", Nx
  write(10,*) "Nodes: "
  do i = 1, size(nodes(:,1))
    write(10,*) 'Node ', i, "; X: ", nodes(i,1)
  end do
  write(10,*)"------------------------------------------------------------------"
  write(10,*) "Vectors U"
  do i = 1, size(vectors(:,1,1))
    write(10,*) 'Vector ', i, "; X: ", vectors(i,1, 1)
  end do


else if (dimension == 2) then


  count = 0
  do x = 1, Nx + 2
    do y = 1, Ny + 2
      count =  y + (Ny+2)*(x-1)

      nodes(count, 1) = (x-1.5)*stepx
      nodes(count, 2) = (y-1.5)*stepy

      if (x == 1) then
        nodes(count, 1) = 0
      else if (x == nx+2) then
        nodes(count, 1) = Lx
      end if

      if (y == 1) then
        nodes(count, 2) = 0
      else if (y == Ny+2) then
        nodes(count, 2) = Ly
      end if
      end do
    end do

    do x = 1, Nx + 1
      do y = 1, Ny + 2
        count =  y + (Ny+2)*(x-1)

        vectors(count, 1, 1) = (x-1)*stepx
        vectors(count, 1, 2) = (y-1.5)*stepy

        if (y == 1) then
          vectors(count, 1, 2) = 0
        else if (y == Ny + 2) then
          vectors(count, 1, 2) = Ly
        end if
      end do
    end do

    do x = 1, Nx + 2
      do y = 1, Ny + 1
        count = y + (Ny+1)*(x-1)

        vectors(count, 2, 1) = (x-1.5)*stepx
        vectors(count, 2, 2) = (y-1)*stepy

        if (x == 1) then
          vectors(count, 2, 1) = 0
        else if (x == Nx + 2) then
          vectors(count, 2, 1) = Lx
        end if
      end do
    end do

    open(10, file = 'mesh.txt', status = 'replace')

    write(10,*) "Mesh"
    write(10,*) "Dimension: ", dimension
    write(10,*) "Lx = ", Lx, "; Ly = ", Ly
    write(10,*) "Nx = ", Nx, "; Ny = ", Ny
    write(10,*) "Nodes: "

    do i = 1, size(nodes(:,1))
      write(10,*) 'Node ', i, "; X: ", nodes(i,1), "; Y: ", nodes(i, 2)
    end do
    write(10,*)"------------------------------------------------------------------"
    write(10,*) "Vectors U"
    do i = 1, size(vectors(:,1,1))
      write(10,*) 'Vector ', i, "; X: ", vectors(i,1, 1), "; Y: ", vectors(i, 1, 2)
    end do
    write(10,*)"------------------------------------------------------------------"
    write(10,*) "Vectors V"
    do i = 1, size(vectors(:,2,1))
      write(10,*) 'Vector ', i, "; X: ", vectors(i,2, 1), "; Y: ", vectors(i, 2, 2)
    end do

else if (dimension == 3) then



  do x = 1, Nx + 2
    do y = 1, Ny + 2
      do z = 1, Nz + 2

        count =  z +(y-1)*(Nz+2) + (Nz+2)*(Ny+2)*(x-1)

        nodes(count, 1) = (x-1.5)*stepx
        nodes(count, 2) = (y-1.5)*stepy
        nodes(count, 3) = (z-1.5)*stepz

        if (x == 1) then
          nodes(count, 1) = 0
        else if (x == nx+2) then
          nodes(count, 1) = Lx
        end if

        if (y == 1) then
          nodes(count, 2) = 0
        else if (y == Ny+2) then
          nodes(count, 2) = Ly
        end if

        if (z == 1) then
          nodes(count, 3) = 0
        else if (z == Nz+2) then
          nodes(count, 3) = Lz
        end if

        end do
      end do
    end do

  do x = 1, Nx + 1
    do y = 1, Ny + 2
      do z = 1, Nz + 2
        count =  z + (y-1)*(Nz+2) + (Ny+2)*(x-1)*(Nz+2)

        vectors(count, 1, 1) = (x-1)*stepx
        vectors(count, 1, 2) = (y-1.5)*stepy
        vectors(count, 1, 3) = (z-1.5)*stepz

        if (y == 1) then
          vectors(count, 1, 2) = 0
        else if (y == Ny + 2) then
          vectors(count, 1, 2) = Ly
        end if
        if (z == 1) then
          vectors(count, 1, 3) = 0
        else if (z == Nz + 2) then
          vectors(count, 1, 3) = Lz
        end if
      end do
    end do
  end do

  do y = 1, Ny + 1
    do x = 1, Nx + 2
      do z = 1, Nz + 2
        count =  z + (x-1)*(Nz+2) + (Nx+2)*(y-1)*(Nz+2)

        vectors(count, 2, 2) = (y-1)*stepy
        vectors(count, 2, 1) = (x-1.5)*stepx
        vectors(count, 2, 3) = (z-1.5)*stepz

        if (x == 1) then
          vectors(count, 2, 1) = 0
        else if (x == Nx + 2) then
          vectors(count, 2, 1) = Lx
        end if
        if (z == 1) then
          vectors(count, 2, 3) = 0
        else if (z == Nz + 2) then
          vectors(count, 2, 3) = Lz
        end if
      end do
    end do
  end do

  do z = 1, Nz + 1
    do y = 1, Ny + 2
      do x = 1, Nx + 2
        count =  x + (y-1)*(Nx+2) + (Ny+2)*(z-1)*(Nx+2)

        vectors(count, 3, 3) = (z-1)*stepz
        vectors(count, 3, 2) = (y-1.5)*stepy
        vectors(count, 3, 1) = (x-1.5)*stepx

        if (y == 1) then
          vectors(count, 3, 2) = 0
        else if (y == Ny + 2) then
          vectors(count, 3, 2) = Ly
        end if
        if (x == 1) then
          vectors(count, 3, 1) = 0
        else if (x == Nx + 2) then
          vectors(count, 3, 1) = Lx
        end if
      end do
    end do
  end do

  open(10, file = 'mesh.txt', status = 'replace')

  write(10,*) "Mesh"
  write(10,*) "Dimension: ", dimension
  write(10,*) "Lx = ", Lx, "; Ly = ", Ly, "; Lz = ", Lz
  write(10,*) "Nx = ", Nx, "; Ny = ", Ny, "; Nz = ", Nz
  write(10,*) "Nodes: "

  do i = 1, size(nodes(:,1))
    write(10,*) 'Node ', i, "; X: ", nodes(i,1), "; Y: ", nodes(i, 2),"; Z: ", nodes(i,3)
  end do
  write(10,*)"------------------------------------------------------------------"
  write(10,*) "Vectors U"
  do i = 1, size(vectors(:,1,1))
    write(10,*) 'Vector ', i, "; X: ", vectors(i,1, 1), "; Y: ", vectors(i, 1, 2),"; Z: ", vectors(i,1,3)
  end do
  write(10,*)"------------------------------------------------------------------"
  write(10,*) "Vectors V"
  do i = 1, size(vectors(:,2,1))
    write(10,*) 'Vector ', i, "; X: ", vectors(i,2, 1), "; Y: ", vectors(i, 2, 2),"; Z: ", vectors(i,2,3)
  end do
  write(10,*)"------------------------------------------------------------------"
  write(10,*) "Vectors W"
  do i = 1, size(vectors(:,3,1))
    write(10,*) 'Vector ', i, "; X: ", vectors(i,3,1), "; Y: ", vectors(i, 3, 2),"; Z: ", vectors(i,3,3)
  end do

end if

return
end subroutine mesher
end module mesh
