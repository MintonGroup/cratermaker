!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module surface
   use iso_c_binding
   use globals

   type surface_tally
      real(DP), dimension(:,:), allocatable :: diam ! crater diameter
      real(SP), dimension(:,:), allocatable :: xl,yl ! Crater center 
   contains
      procedure :: allocate   => surface_tally_allocate   !! Allocate the allocatable components of the class
      procedure :: deallocate => surface_tally_deallocate !! Deallocate all allocatable components of the class
      final     ::               surface_tally_final      !! Finalizer (calls deallocate)
   end type surface_tally 

   type surface_type
      type(surface_tally), allocatable(:)   :: tally
      real(DP), dimension(:,:), allocatable :: elev  ! surface elevation
   contains
      procedure :: allocate   => surface_allocate   !! Allocate the allocatable components of the class
      procedure :: deallocate => surface_deallocate !! Deallocate all allocatable components of the class
      final     ::               surface_final      !! Finalizer (calls deallocate)
   end type surface_type


contains

   subroutine surface_allocate(self, gridsize, nlayers)
      !! author: David A. Minton
      !!
      !! Allocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_type), intent(inout) :: self     !! Surface object
      integer(I4B),        intent(in)    :: gridsize !! Size of the grid
      integer(I4B),        intent(in)    :: nlayers  !! Number of tally layers
      ! Internals
      integer(I4B) :: layer

      allocate(self%elev(gridsize,gridsize))
      allocate(self%tally(nlayers))
      do layer = 1, nlayers
         call self%tally(layer)%allocate(gridsize)
      end do

      return
   end subroutine surface_allocate


   subroutine surface_deallocate(self)
      !! author: David A. Minton
      !!
      !! Deallocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_type), intent(inout) :: self     !! Surface object

      deallocate(self%elev)
      call self%tally%deallocate()
      deallocate(self%tally)

      return
   end subroutine surface_allocate


   subroutine surface_final(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the surface object
      implicit none
      ! Arguments
      class(surface_type), intent(inout) :: self

      call self%deallocate()
      return
   end subroutine surface_final


   subroutine surface_tally_allocate(self, gridsize)
      !! author: David A. Minton
      !!
      !! Allocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_tally), intent(inout) :: self     !! Surface object
      integer(I4B),        intent(in)    :: gridsize !! Size of the grid
      ! Internals
      integer(I4B) :: layer

      allocate(self%diam(gridsize,gridsize))
      allocate(self%xl(nlayers))
      allocate(self%yl(nlayers))

      return
   end subroutine surface_tally_allocate


   subroutine surface_tally_deallocate(self)
      !! author: David A. Minton
      !!
      !! Deallocate the allocatable components of the class
      implicit none
      ! Arguments
      class(surface_type), intent(inout) :: self     !! Surface tally object

      deallocate(self%diam)
      deallocate(self%xl)
      deallocate(self%yl)

      return
   end subroutine surface_tally_allocate


   subroutine surface_tally_final(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the surface tally object
      implicit none
      ! Arguments
      class(surface_tally), intent(inout) :: self

      call self%deallocate()
      return
   end subroutine surface_tally_final


end module surface