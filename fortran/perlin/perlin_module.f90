!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module perlin
   use globals

   interface
      module pure elemental subroutine perlin_noise(xx, yy, zz, noise, dx, dy, dz)
         implicit none
         real(DP),intent(in) :: xx,yy,zz
         real(DP),intent(out) :: noise
         real(DP),intent(out),optional :: dx, dy, dz
      end subroutine perlin_noise


      module pure function perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, z, noise_height, freq, pers
         integer(I4B), intent(in) :: num_octaves
         real(DP), dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function perlin_turbulence
   
   
      module pure function perlin_billowedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, z, noise_height, freq, pers
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function perlin_billowedNoise
         
         
      module pure function perlin_plawNoise(x, y, z, noise_height, freq, pers, slope, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, z, noise_height, freq, pers,slope
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function perlin_plawNoise
         
         
      module pure function perlin_ridgedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, z, noise_height, freq, pers
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:),intent(in) :: anchor
         real(DP) :: noise
      end function perlin_ridgedNoise
         
         
      module pure function perlin_swissTurbulence(x, y, z, lacunarity, gain, warp, num_octaves, anchor) result(noise)
         implicit none
         real(DP), intent(in) ::  x, y, z, lacunarity, gain, warp
         integer(I4B), intent(in) :: num_octaves
         real(DP),dimension(:,:), intent(in) :: anchor
         real(DP) :: noise
      end function perlin_swissTurbulence
      
         
      module pure function perlin_jordanTurbulence(x, y, z, lacunarity, gain0, gain, warp0, warp, damp0, damp, damp_scale,&
               num_octaves, anchor) result(noise)
         implicit none
         real(DP),intent(in) :: x, y, z, lacunarity, gain0, gain, warp0, warp, damp0, damp, damp_scale
         integer(I4B),intent(in) :: num_octaves
         real(DP), dimension(:,:), intent(in) :: anchor
         real(DP) :: noise
      end function perlin_jordanTurbulence
         
   end interface

end module perlin