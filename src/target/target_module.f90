!! Copyright 2023 - David Minton
!! This file is part of PyOOF
!! PyOOF is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! pyoof is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with pyoof. 
!! If not, see: https://www.gnu.org/licenses. 

module target_body
   use globals

   type  :: target_body_type 
      real(DP), dimension(:,:), allocatable :: elevation    !! Elevation of surface mesh relative to the datum
      character(len=STRMAX)                 :: name    !! A placeholder for a string component variable
   contains
      procedure :: allocate   => target_body_allocate   !! Allocate the allocatable components of the class
      procedure :: deallocate => target_body_deallocate !! Deallocate all allocatable components of the class
      final     ::               target_body_final      !! Finalizer (calls deallocate)
   end type target_body_type


contains

   subroutine target_body_allocate(self, nx, ny)
      !! author: David A. Minton
      !!
      !! Allocate the allocatable components of the class
      implicit none
      ! Arguments
      class(target_body_type), intent(inout) :: self   !! Simulation object
      integer(I4B),        intent(in)    :: nx, ny !! Size of the grid

      allocate(self%elevation(nx,ny))

      self%elevation(:,:) = 0.0_DP
      self%name = "None"

      return
   end subroutine target_body_allocate


   subroutine target_body_deallocate(self) 
      !! author: David A. Minton
      !!
      !! Deallocate the allocatable components of the class
      implicit none
      ! Arguments
      class(target_body_type), intent(inout) :: self !! Surface object

      deallocate(self%elevation)

      return
   end subroutine target_body_deallocate


   subroutine target_body_final(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the surface object
      implicit none
      ! Arguments
      type(target_body_type), intent(inout) :: self

      call self%deallocate()
      return
   end subroutine target_body_final


end module target_body