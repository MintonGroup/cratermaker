! Copyright 2024 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

module ejecta
    use globals
    use bind
 
    interface
        pure module subroutine ejecta_profile(radial_distance, diameter, ejrim, elevation)
        implicit none
        real(DP),dimension(:), intent(in) :: radial_distance
        real(DP), intent(in) :: diameter, ejrim
        real(DP), dimension(:), intent(out) :: elevation
        end subroutine ejecta_profile
    end interface

contains
    subroutine bind_ejecta_profile(c_radial_distance, num_elements, diameter, ejrim, c_elevation) bind(C)
       ! Arguments
       type(c_ptr), intent(in), value :: c_radial_distance
       integer(I4B), intent(in), value :: num_elements
       real(DP), intent(in), value :: diameter, ejrim
       type(c_ptr), intent(in), value :: c_elevation
       ! Internals
       real(DP), dimension(:), pointer :: radial_distance,elevation
       integer(I4B) :: i
       
       if (c_associated(c_radial_distance)) then
          call c_f_pointer(c_radial_distance, radial_distance,shape=[num_elements]) 
       else
          return
       end if
 
       allocate(elevation(num_elements))
       if (c_associated(c_elevation)) then
          call c_f_pointer(c_elevation, elevation, shape=[num_elements]) 
       else
          return
       end if
 
       call ejecta_profile(radial_distance, diameter, ejrim, elevation)

       return
    end subroutine bind_ejecta_profile

end module ejecta