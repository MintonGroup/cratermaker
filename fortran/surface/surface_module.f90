!! Copyright 2023 - David Minton
!! This file is part of PyOOF
!! PyOOF is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! pyoof is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with pyoof. 
!! If not, see: https://www.gnu.org/licenses. 

module surface
   use globals

   type  :: surface_type 
      real(DP), dimension(:), allocatable :: elevation !! Elevation of surface mesh relative to the datum
      real(DP), dimension(:), allocatable :: node_x    !! Cartesian position of the nodes in the x-direction
      real(DP), dimension(:), allocatable :: node_y    !! Cartesian position of the nodes in the y-direction
      real(DP), dimension(:), allocatable :: node_z    !! Cartesian position of the nodes in the z-direction
      real(DP), dimension(:), allocatable :: face_x    !! Cartesian position of the faces in the x-direction
      real(DP), dimension(:), allocatable :: face_y    !! Cartesian position of the faces in the y-direction
      real(DP), dimension(:), allocatable :: face_z    !! Cartesian position of the faces in the z-direction
   contains
      procedure :: allocate   => surface_allocate   !! Allocate the allocatable components of the class
      procedure :: deallocate => surface_deallocate !! Deallocate all allocatable components of the class
      final     ::               surface_final      !! Finalizer (calls deallocate)
   end type surface_type


contains

   subroutine surface_allocate(self, n_node, n_face)
      !! author: David A. Minton
      !!
      !! Allocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_type), intent(inout), target :: self   !! Simulation object
      integer(I4B),        intent(in)    :: n_node !! Number of nodes in the surface mesh
      integer(I4B),        intent(in)    :: n_face !! Number of faces in the surface mesh
      type(c_ptr) :: cptr

      allocate(self%elevation(n_node))
      allocate(self%node_x(n_node))
      allocate(self%node_y(n_node))
      allocate(self%node_z(n_node))
      allocate(self%face_x(n_face))
      allocate(self%face_y(n_face))
      allocate(self%face_z(n_face))

      self%elevation(:) = 0.0_DP
      self%node_x(:)    = 0.0_DP
      self%node_y(:)    = 0.0_DP
      self%node_z(:)    = 0.0_DP
      self%face_x(:)    = 0.0_DP
      self%face_y(:)    = 0.0_DP
      self%face_z(:)    = 0.0_DP

      write(*,*) "Allocated elevation in Fortran with size: ", size(self%elevation)
      write(*,*) "Pointer address of elevation in Fortran: ", c_loc(self%elevation)
      write(*,*) "Allocated node_x in Fortran with size: ", size(self%node_x)
      write(*,*) "Pointer address of node_x in Fortran: ", c_loc(self%node_x)
      write(*,*) "Allocated node_y in Fortran with size: ", size(self%node_y)
      write(*,*) "Pointer address of node_y in Fortran: ", c_loc(self%node_y)
      write(*,*) "Allocated node_z in Fortran with size: ", size(self%node_z)
      write(*,*) "Pointer address of node_z in Fortran: ", c_loc(self%node_z)
      write(*,*) "Allocated face_x in Fortran with size: ", size(self%face_x)
      write(*,*) "Pointer address of face_x in Fortran: ", c_loc(self%face_x)
      write(*,*) "Allocated face_y in Fortran with size: ", size(self%face_y)
      write(*,*) "Pointer address of face_y in Fortran: ", c_loc(self%face_y)
      write(*,*) "Allocated face_z in Fortran with size: ", size(self%face_z)
      write(*,*) "Pointer address of face_z in Fortran: ", c_loc(self%face_z)


      return
   end subroutine surface_allocate


   subroutine surface_deallocate(self) 
      !! author: David A. Minton
      !!
      !! Deallocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_type), intent(inout) :: self !! Surface object

      deallocate(self%elevation)
      deallocate(self%node_x)
      deallocate(self%node_y)
      deallocate(self%node_z)
      deallocate(self%face_x)
      deallocate(self%face_y)
      deallocate(self%face_z)

      return
   end subroutine surface_deallocate


   subroutine surface_final(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the surface object
      implicit none
      ! Arguments
      type(surface_type), intent(inout) :: self

      call self%deallocate()
      return
   end subroutine surface_final


end module surface