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
      real(DP), dimension(:,:), allocatable :: elevation    !! Elevation of surface
      character(len=STRMAX)                 :: stringvar    !! A placeholder for a string component variable
   contains
      procedure :: allocate   => surface_allocate   !! Allocate the allocatable components of the class
      procedure :: deallocate => surface_deallocate !! Deallocate all allocatable components of the class
      final     ::               surface_final      !! Finalizer (calls deallocate)
   end type surface_type


contains

   subroutine surface_allocate(self, nx, ny)
      !! author: David A. Minton
      !!
      !! Allocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_type), intent(inout) :: self   !! Simulation object
      integer(I4B),        intent(in)    :: nx, ny !! Size of the grid

      allocate(self%elevation(nx,ny))

      self%elevation(:,:) = -1.0_DP
      self%stringvar = "Initialized in Fortran"

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