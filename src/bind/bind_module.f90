!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module bind_module
   use iso_c_binding
   use globals
   use surface
   implicit none

contains

   type(c_ptr) function bind_surface_init(gridsize) bind(c)
      integer(I4B), value :: gridsize
      type(surface_type), pointer :: surf_ptr

      allocate(surf_ptr)
      call surf_ptr%allocate(gridsize) 
      bind_surface_init = c_loc(surf_ptr)
   end function bind_surface_init

   subroutine bind_surface_final(surf) bind(c)
      type(c_ptr), intent(in), value :: surf
      type(surface_type), pointer :: surf_ptr

      call c_f_pointer(surf, surf_ptr)
      deallocate(surf_ptr)
   end subroutine bind_surface_final


end module bind_module