!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 


module bind_module
   !! author: David A. Minton
   !!
   !! This module defines the set of routines that connect the Cython code to the Fortran. Because Fortran derived types with
   !! allocatable arrays and type-bound procedures are not interoperable with C, you have to write a set of functions that allow
   !! you to access the Fortran structures. On the Fortran side, you can make use of the full suite of OOF features from F2003+ 
   !! e.g. classes, inheritence, polymorphism, type-bound procedures, etc., but any communication back to Python must be done 
   !! via a set of binding functions. 
   !! 
   !! The following implementation was adapted from _Modern Fortran Explained: Incorporating Fortran 2018_ by Metcalf, Reid, & 
   !! Cohen (see Fig. 19.8)
   use iso_c_binding
   use globals
   use simulation
   use target_body
   use util
   implicit none
   
contains

   type(c_ptr) function bind_body_init(ny,nx) bind(c)
      !! author: David A. Minton
      !!
      !! This function is used to initialize the target_body_type derived type object in Fortran and return a pointer to the object 
      !! that can be used as a struct in C, and ultimately to the Python class object via Cython.
      implicit none
      ! Arguments
      integer(I4B), intent(in), value   :: ny, nx !! The dimensions of the array to create. This expect C (row-major) ordering 
      ! Internals
      type(target_body_type), pointer :: f_sim  !! A pointer to the body type variable that will be passed to Cython

      nullify(f_sim)
      allocate(f_sim)
      call f_sim%allocate(nx, ny)
      bind_body_init = c_loc(f_sim)

      return
   end function bind_body_init


   subroutine bind_body_final(sim) bind(c)
      !! author: David A. Minton
      !!
      !! This subroutine is used to deallocate the pointer that links the C struct to the Fortran derived type object. 
      implicit none
      ! Arguments
      type(c_ptr), intent(in), value :: sim !! C pointer to the Fortran body object
      ! Internals
      type(target_body_type), pointer :: f_sim

      if (c_associated(sim)) then
         call c_f_pointer(sim, f_sim)
         deallocate(f_sim)
      end if

      return
   end subroutine bind_body_final


   type(c_ptr) function bind_body_get_name(c_sim) bind(c)
      !! author: David A. Minton
      !!
      !! This function is used to retrieve the string variable from the Fortran body derived-type and pass it to a to C, and
      !! ultimately to the Python string variable via Cython.
      implicit none
      ! Arguments
      type(c_ptr), value :: c_sim  !! C pointer to the Fortran body object
      ! Internals
      type(target_body_type), pointer :: f_sim  !! A pointer to the body type variable that will be passed to Cython
      character(kind=c_char), dimension(STRMAX), target :: f_str

      if (c_associated(c_sim)) then
         nullify(f_sim)
         call c_f_pointer(c_sim, f_sim)
         call bind_f2c_string(trim(f_sim%name), f_str)
         bind_body_get_name = c_loc(f_str)
      else
         write(*,*) "The c_sim pointer is NULL!"
         bind_body_get_name = c_null_ptr
      end if

      return
   end function bind_body_get_name


   subroutine bind_body_set_name(c_sim, c_string) bind(c)
      !! author: David A. Minton
      !!
      !! This subroutine is used to set value of the string variable in the Fortran body derived-type. It takes a C-style 
      !! string and converts it to a Fortran-style string 
      implicit none
      ! Arguments
      type(c_ptr), value,                   intent(in) :: c_sim  !! C pointer to the Fortran body object
      character(kind=c_char), dimension(*), intent(in) :: c_string  !! Input C-style string
      ! Internals
      character(len=STRMAX)  :: f_string
      type(target_body_type), pointer :: f_sim

      nullify(f_sim)
      if (c_associated(c_sim)) then
         call c_f_pointer(c_sim, f_sim)
         call bind_c2f_string(c_string, f_string)
         f_sim%name = f_string
         write(*,*) "The Fortran derived type has received a new string variable value: ",trim(f_sim%name)
      else
         write(*,*) "The c_sim pointer is NULL!"
      end if

      return
   end subroutine bind_body_set_name


   subroutine bind_c2f_string(c_string, f_string) 
      !! author: David A. Minton
      !!
      !! This subroutine is used to convert C style strings into Fortran. This allows one to pass Python strings as arguments to 
      !! Fortran functions.
      implicit none
      ! Arguments
      character(kind=c_char), dimension(*), intent(in)  :: c_string  !! Input C-style string
      character(len=STRMAX),                intent(out) :: f_string  !! Output Fortran-style string
      ! Internals
      integer :: i
      character(len=STRMAX,kind=c_char) :: tmp_string

      i=1
      tmp_string = ''
      do while(c_string(i) /= c_null_char .and. i <= STRMAX)
         tmp_string(i:i) = c_string(i)
         i=i+1
      end do

      if (i > 1) then
         f_string = trim(tmp_string)
      else
         f_string = ""
      end if

      return
   end subroutine bind_c2f_string


   subroutine bind_f2c_string(f_string, c_string)
      !! author: David A. Minton
      !!
      !! This subroutine is used to convert Fortran style strings to C. This allows the Python module to read strings that were 
      !! created in Fortran procedures.
      implicit none
      ! Arguments
      character(len=*), intent(in) :: f_string  !! Input Fortran-style string
      character(kind=c_char), dimension(STRMAX), target, intent(out) :: c_string  !! Output C-style string (array of C characters)
      ! Internals 
      integer :: n, i

      n = min(len(f_string), STRMAX-1)
      do i = 1, n
         c_string(i) = f_string(i:i)
      end do
      c_string(n + 1) = c_null_char
   end subroutine bind_f2c_string


   function bind_perlin_noise(c_model, x, y, z, num_octaves, c_anchor, damp, damp0, damp_scale, freq, gain, gain0, lacunarity, &
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

      select case (trim(model))
      case('turbulence')
         noise = util_perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor) 
      case('billowed')
         noise = util_perlin_billowedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor)
      case('plaw')
         noise = util_perlin_plawNoise(x, y, z, noise_height, freq, pers, slope, num_octaves, anchor)
      case('ridged')
         noise = util_perlin_ridgedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor)
      case('swiss')
         noise = util_perlin_swissTurbulence(x, y, z, lacunarity, gain, warp, num_octaves, anchor)
      case('jordan')
         noise = util_perlin_jordanTurbulence(x, y, z, lacunarity, gain0, gain, warp0, warp, damp0, damp,& 
                                                damp_scale, num_octaves, anchor) 
      case default
         noise = 0.0_DP
      end select
      return
  end function bind_perlin_noise


  

   
end module bind_module