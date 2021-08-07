module mesh

implicit none

contains

subroutine mesher(Lx, Ly, Nx, Ny, vectors, nodes)

implicit none

real :: Lx, Ly
integer :: Nx, Ny
real :: stepx, stepy
real, allocatable :: nodes (:, :), vectors(:, :, :)
integer :: count, x, y, mesh, number_nodes, scalar, vector_dim, i

!distance between nodes
stepx = Lx/Nx
stepy = Ly/Ny


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
  write(10,*) "Dimension: ", 2
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


return
end subroutine mesher
end module mesh
