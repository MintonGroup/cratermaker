! Copyright 2024 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

!! The interfaces for the Perlin noise routines.
module crater
   use globals
   use bind

   interface
      pure module subroutine crater_profile(radial_distance, reference_elevation, diameter, &
         floordepth, floor_diameter, rimheight, ejrim, elevation)
         implicit none
         real(DP),dimension(:), intent(in) :: radial_distance, reference_elevation
         real(DP), intent(in) :: diameter, floordepth, floor_diameter, rimheight, ejrim
         real(DP), dimension(:), intent(out) :: elevation
      end subroutine crater_profile
   end interface

contains
   subroutine bind_crater_profile(c_radial_distance, c_reference_elevation, num_elements, diameter, &
                                             floordepth, floor_diameter, rimheight, ejrim, c_elevation)  &
                                             bind(C)
      ! Arguments
      type(c_ptr), intent(in), value :: c_radial_distance, c_reference_elevation
      integer(I4B), intent(in), value :: num_elements
      real(DP), intent(in), value :: diameter, floordepth, floor_diameter, rimheight, ejrim
      type(c_ptr), intent(in), value :: c_elevation
      ! Internals
      real(DP), dimension(:), pointer :: radial_distance,reference_elevation,elevation
      integer(I4B) :: i
      
      if (c_associated(c_radial_distance)) then
         call c_f_pointer(c_radial_distance, radial_distance,shape=[num_elements]) 
      else
         return
      end if

      if (c_associated(c_reference_elevation)) then
         call c_f_pointer(c_reference_elevation, reference_elevation,shape=[num_elements]) 
      else
         return
      end if

      allocate(elevation(num_elements))
      if (c_associated(c_elevation)) then
         call c_f_pointer(c_elevation, elevation, shape=[num_elements]) 
      else
         return
      end if

      call crater_profile(radial_distance, reference_elevation, diameter, floordepth, &
                              floor_diameter, rimheight, ejrim, elevation)
      return
   end subroutine bind_crater_profile

end module crater