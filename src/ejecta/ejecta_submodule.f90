! Copyright 2023 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

!! The implementations for the Perlin noise procedures.
submodule (ejecta) s_ejecta
    use globals
contains

    pure module subroutine ejecta_profile(radial_distance, diameter, ejrim, elevation)
        !! author: David A. Minton
        !!
        !! Calculate the elevation of a crater as a function of distance from the center.
        !!
        implicit none
        ! Arguments
        real(DP),dimension(:), intent(in) :: radial_distance  
            !! Radial distance from the crater center in meters. 
        real(DP), intent(in) :: diameter 
            !! The final crater diameter in meters.
        real(DP), intent(in) :: ejrim 
            !! The final ejecta rim height in meters.
        real(DP), dimension(:), intent(out) :: elevation 
            !! The elevation of the crater relative to a reference surface.
        ! Internals
        real(DP) :: r, radius

        ! Calculate the floor radius relative to the final crater radius
        radius = diameter / 2

        where(radial_distance(:) >= radius)
            elevation(:) = elevation(:) + ejecta_profile_func(radial_distance(:)/radius)
        end where

        return         

        contains

        pure elemental function ejecta_profile_func(r) result(h)
            !! author: David A. Minton
            !!
            !! Calculate the elevation of the ejecta layer as a function of distance from the center.
            !!                                
            implicit none
            ! Arguments
            real(DP), intent(in) :: r 
                !! Radial distance from the crater center in meters. 
            real(DP) :: h
                !! Results

            ! Internals
            real(DP), parameter :: ejprofile = 3.0_DP

            h = ejrim * (r)**(-ejprofile)

        end function ejecta_profile_func

   end subroutine ejecta_profile
end submodule s_ejecta