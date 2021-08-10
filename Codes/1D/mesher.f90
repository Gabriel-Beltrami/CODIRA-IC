module mesher

implicit none

contains

subroutine mesher1d(Lx, Nx, vectors, nodes)

implicit none

real :: Lx
integer :: Nx
real :: stepx
real, allocatable :: nodes (:, :), vectors(:, :, :)
integer :: count, x, mesh, number_nodes, scalar, vector_dim, i

!distance between nodes
stepx = Lx/Nx


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
write(10,*) "Dimension: ", 1
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


return
end subroutine mesher1d
end module mesher
