! Copyright 2023 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

!! The implementations for the Perlin noise procedures.
submodule (morphology_module) s_morphology
    use globals_module
contains

   module pure subroutine morphology_profile(radial_distance, reference_elevation, diameter, floordepth, floordiam, &
      rimheight, ejrim, RIMDROP, elevation)
      !! author: David A. Minton
      !!
      !! Calculate the elevation of a crater as a function of distance from the center.
      !!
      implicit none
      ! Arguments
      real(DP),dimension(:), intent(in) :: radial_distance  !! Radial distance from the crater center in meters. 
      real(DP),dimension(:), intent(in) :: reference_elevation !! The elevation of the reference surface at the radial distance r.
      real(DP), intent(in) :: diameter !! The final crater diameter in meters.
      real(DP), intent(in) :: floordepth !! The final crater floor depth in meters.
      real(DP), intent(in) :: floordiam !! The final crater floor diameter in meters.
      real(DP), intent(in) :: rimheight !! The final crater rim height in meters.
      real(DP), intent(in) :: ejrim !! The final ejecta rim height in meters.
      real(DP), intent(in) :: RIMDROP !! The exponent for the ejecta rim dropoff.
      real(DP), dimension(:), intent(out) :: elevation !! The elevation of the crater relative to a reference surface.
      ! Internals
      integer(I4B) :: i, ninc
      real(DP) :: r, flrad, radius, meanref, min_elevation
      real(DP) :: c0, c1, c2, c3
      real(DP), parameter :: A = 4.0_DP / 11.0_DP
      real(DP), parameter :: B = -32.0_DP / 187.0_DP
      real(DP), parameter :: ejprofile = 3.0_DP

      ! Calculate the floor radius relative to the final crater radius
      flrad = floordiam / diameter
      radius = diameter / 2

      ! Use polynomial crater profile similar to that of Fassett et al. (2014), but the parameters are set by the crater dimensions
      c1 = (-floordepth - rimheight) / (flrad - 1.0_DP + A * (flrad**2 - 1.0_DP) + B * (flrad**3 - 1.0_DP))
      c0 = rimheight - c1 * (1.0_DP + A + B)
      c2 = A * c1
      c3 = B * c1
      ninc = count(radial_distance(:) <= radius)
      if (ninc == 0) then
         meanref = reference_elevation(minloc(radial_distance(:),dim=1))
      else
         meanref = sum(reference_elevation(:), radial_distance(:) <= radius) / ninc
      end if
      min_elevation = meanref - floordepth

      elevation(:) = reference_elevation(:) + crater_profile(radial_distance(:)/radius)
      where(radial_distance(:) >= radius)
         elevation(:) = elevation(:) + ejecta_profile(radial_distance(:)/radius)
      elsewhere(elevation(:) < min_elevation)
         elevation(:) = min_elevation
      end where

      return         

      contains

         pure elemental function crater_profile(r) result(h)
            !! author: David A. Minton
            !!
            !! Calculate the elevation of a crater as a function of distance from the center.
            !!                                
            implicit none
            ! Arguments
            real(DP), intent(in) :: r !! Radial distance from the crater center in meters. 

            ! Results
            real(DP) :: h

            if (r >= 1.0_DP) then
               h = (rimheight - ejrim) * (r**(-RIMDROP))
            else
               h = c0 + c1 * r + c2 * r**2 + c3 * r**3
            end if

         end function crater_profile

         pure elemental function ejecta_profile(r) result(h)
            !! author: David A. Minton
            !!
            !! Calculate the elevation of the ejecta layer as a function of distance from the center.
            !!                                
            implicit none
            ! Arguments
            real(DP), intent(in) :: r !! Radial distance from the crater center in meters. 

            ! Results
            real(DP) :: h

            h = ejrim * (r)**(-ejprofile)

         end function ejecta_profile

   end subroutine morphology_profile

end submodule s_morphology