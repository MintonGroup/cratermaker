!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module perlin_module
   ! Module that contains the perlin noise functions
   use globals_module
   use bind_module

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

      module pure function perlin_noise_one(model, x, y, z, num_octaves, anchor, damp, damp0, damp_scale, freq, gain, gain0,&
                                            lacunarity, noise_height, pers, slope, warp, warp0) result(noise)
         implicit none
         character(len=*), intent(in) :: model !! The specific turbulence model to apply
         real(DP), intent(in) :: x, y, z  !! The xyz cartesian position of the noise to evaluate
         integer(I4B), intent(in) :: num_octaves
         real(DP), dimension(:,:), intent(in) :: anchor
         real(DP), intent(in) :: damp, damp0, damp_scale, freq, gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0
         real(DP) :: noise
      end function perlin_noise_one


      module subroutine perlin_noise_all(model, x, y, z, num_octaves, anchor, damp, damp0, damp_scale, freq, gain, gain0,&
                                            lacunarity, noise_height, pers, slope, warp, warp0,noise)
         implicit none
         character(len=*), intent(in) :: model !! The specific turbulence model to apply
         real(DP), dimension(:), intent(in) :: x, y, z  !! The xyz cartesian position of the noise to evaluate
         integer(I4B), intent(in) :: num_octaves
         real(DP), dimension(:,:), intent(in) :: anchor
         real(DP), intent(in) :: damp, damp0, damp_scale, freq, gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0
         real(DP), dimension(:), intent(out) :: noise
      end subroutine perlin_noise_all


         
   end interface

contains
   function bind_perlin_noise_one(c_model, x, y, z, num_octaves, c_anchor, damp, damp0, damp_scale, freq, gain, gain0, lacunarity, &
                                 noise_height, pers, slope, warp, warp0) bind(c) result(noise)
      ! Arguments
      character(kind=c_char), dimension(*), intent(in) :: c_model !! The specific turbulence model to apply
      real(DP), intent(in), value :: x, y, z  !! The xyz cartesian position of the noise to evaluate
      real(DP), intent(in), value :: damp, damp0, damp_scale, freq, gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0
      integer(I4B), intent(in), value :: num_octaves
      type(c_ptr), intent(in), value :: c_anchor
      ! Return
      real(DP) :: noise
      ! Internals
      real(DP), dimension(:,:), pointer :: anchor
      character(len=STRMAX) :: model
      

      ! Convert from C-style string to Fortran-style
      call bind_c2f_string(c_model, model)
      if (c_associated(c_anchor)) then
         call c_f_pointer(c_anchor, anchor,shape=[3,num_octaves]) 
      else
         return
      end if

      noise = perlin_noise_one(model, x, y, z, num_octaves, anchor, damp, damp0, damp_scale, freq, gain, gain0, lacunarity, &
                                 noise_height, pers, slope, warp, warp0)
   end function bind_perlin_noise_one


   subroutine bind_perlin_noise_all(c_model, c_x, c_y, c_z, num_elements, num_octaves, c_anchor, damp, damp0, damp_scale, freq, &
                                 gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0, c_noise) bind(c) 
      ! Arguments
      character(kind=c_char), dimension(*), intent(in) :: c_model !! The specific turbulence model to apply
      type(c_ptr), intent(in), value :: c_x, c_y, c_z  !! The xyz cartesian position of the noise to evaluate
      integer(I4B), intent(in), value :: num_elements, num_octaves
      real(DP), intent(in), value :: damp, damp0, damp_scale, freq, gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0
      type(c_ptr), intent(in), value :: c_anchor
      type(c_ptr), intent(in), value :: c_noise
      ! Internals
      real(DP), dimension(:,:), pointer :: anchor
      real(DP), dimension(:), pointer :: x,y,z,noise
      character(len=STRMAX) :: model
      integer(I4B) :: i
      
      ! Convert from C-style string to Fortran-style
      call bind_c2f_string(c_model, model)
      if (c_associated(c_anchor)) then
         call c_f_pointer(c_anchor, anchor,shape=[3,num_octaves]) 
      else
         return
      end if

      if (c_associated(c_x)) then
         call c_f_pointer(c_x, x,shape=[num_elements]) 
      else
         return
      end if

      if (c_associated(c_y)) then
         call c_f_pointer(c_y, y,shape=[num_elements]) 
      else
         return
      end if

      if (c_associated(c_z)) then
         call c_f_pointer(c_z, z,shape=[num_elements]) 
      else
         return
      end if

      allocate(noise(num_elements))
      if (c_associated(c_noise)) then
         call c_f_pointer(c_noise, noise,shape=[num_elements]) 
      else
         return
      end if

      call perlin_noise_all(model, x, y, z, num_octaves, anchor, damp, damp0, damp_scale, freq, gain, gain0, lacunarity, &
                                 noise_height, pers, slope, warp, warp0, noise)

      return
   end subroutine bind_perlin_noise_all

end module perlin_module