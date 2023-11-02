!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module util
   use globals

   interface
      module subroutine util_perlin_noise(xx, yy, noise, dx, dy)
         implicit none
         real(DP),intent(in) :: xx,yy
         real(DP),intent(out) :: noise
         real(DP),intent(out),optional :: dx, dy
      end subroutine util_perlin_noise


      module function util_perlin_turbulence(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, noise_height, freq, pers
         integer(I4B), intent(in) :: num_octaves
         real(DP), dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function util_perlin_turbulence
   
   
      module function util_perlin_billowedNoise(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, noise_height, freq, pers
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function util_perlin_billowedNoise
         
         
      module function util_perlin_plawNoise(x, y, noise_height, freq, pers, slope, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, noise_height, freq, pers,slope
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function util_perlin_plawNoise
         
         
      module function util_perlin_ridgedNoise(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, noise_height, freq, pers
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function util_perlin_ridgedNoise
         
         
      module function util_perlin_swissTurbulence(x, y, lacunarity, gain, warp, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, lacunarity, gain, warp
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:), intent(in) :: anchor
         real(DP) :: noise
      end function util_perlin_swissTurbulence
      
         
      module function util_perlin_jordanTurbulence(x, y, lacunarity, gain0, gain, warp0, warp, damp0, damp, damp_scale,&
               num_octaves, anchor) result(noise)
         implicit none
         real(DP),intent(in) :: x, y, lacunarity, gain0, gain, warp0, warp, damp0, damp, damp_scale
         integer(I4B),intent(in) :: num_octaves
         real(DP), dimension(:,:), intent(in) :: anchor
         real(DP) :: noise
      end function util_perlin_jordanTurbulence
         
   end interface

end module util