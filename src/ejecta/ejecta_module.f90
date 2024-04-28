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
         pure module subroutine ejecta_profile(radial_distance, crater_diameter, ejrim, ejecta_thickness)
             implicit none
             real(DP),dimension(:), intent(in) :: radial_distance
             real(DP), intent(in) :: crater_diameter, ejrim
             real(DP), dimension(:), intent(out) :: ejecta_thickness
         end subroutine ejecta_profile

         module subroutine ejecta_ray_pattern(radial_distance, initial_bearing, crater_diameter, ejrim, ejecta_truncation, &
                                              ejecta_thickness)
             implicit none
             real(DP), dimension(:), intent(in) :: radial_distance, initial_bearing
             real(DP), intent(in) :: crater_diameter, ejrim, ejecta_truncation
             real(DP), dimension(:), intent(out) :: ejecta_thickness
         end subroutine ejecta_ray_pattern

    end interface

contains
    subroutine bind_ejecta_profile(c_radial_distance, num_elements, crater_diameter, ejrim, c_ejecta_thickness) bind(C)
        ! Arguments
        type(c_ptr), intent(in), value :: c_radial_distance
        integer(I4B), intent(in), value :: num_elements
        real(DP), intent(in), value :: crater_diameter, ejrim
        type(c_ptr), intent(in), value :: c_ejecta_thickness
        ! Internals
        real(DP), dimension(:), pointer :: radial_distance,ejecta_thickness
        integer(I4B) :: i
       
        if (c_associated(c_radial_distance)) then
            call c_f_pointer(c_radial_distance, radial_distance,shape=[num_elements]) 
        else
            return
        end if
 
        allocate(ejecta_thickness(num_elements))
        if (c_associated(c_ejecta_thickness)) then
            call c_f_pointer(c_ejecta_thickness, ejecta_thickness, shape=[num_elements]) 
        else
            return
        end if
 
        call ejecta_profile(radial_distance, crater_diameter, ejrim, ejecta_thickness)

        return
    end subroutine bind_ejecta_profile


    subroutine bind_ejecta_ray_pattern(c_radial_distance, c_initial_bearing, num_elements, crater_diameter, ejecta_truncation, &
                                       ejrim, c_ejecta_thickness) bind(c)
        ! Arguments
        type(c_ptr), intent(in), value :: c_radial_distance, c_initial_bearing
        integer(I4B), intent(in), value :: num_elements
        real(DP), intent(in), value :: crater_diameter, ejrim, ejecta_truncation
        type(c_ptr), intent(in), value :: c_ejecta_thickness
        ! Internals
        real(DP), dimension(:), pointer :: radial_distance, initial_bearing, ejecta_thickness

        if (c_associated(c_radial_distance)) then
            call c_f_pointer(c_radial_distance, radial_distance, shape=[num_elements]) 
        else
            return
        end if

        if (c_associated(c_initial_bearing)) then
            call c_f_pointer(c_initial_bearing, initial_bearing, shape=[num_elements]) 
        else
            return
        end if

        allocate(ejecta_thickness(num_elements))
        if (c_associated(c_ejecta_thickness)) then
            call c_f_pointer(c_ejecta_thickness, ejecta_thickness,shape=[num_elements]) 
        else
            return
        end if

        call ejecta_ray_pattern(radial_distance, initial_bearing, crater_diameter, ejrim, ejecta_truncation, ejecta_thickness)
 
        return
    end subroutine bind_ejecta_ray_pattern
end module ejecta