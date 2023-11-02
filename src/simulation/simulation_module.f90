!! Copyright 2023 - David Minton
!! This file is part of PyOOF
!! PyOOF is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! pyoof is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with pyoof. 
!! If not, see: https://www.gnu.org/licenses. 

module simulation
   use globals
   use surface, only : surface_type

   type, extends(surface_type)  :: simulation_type
   contains
      final     ::               simulation_final      !! Finalizer (calls deallocate)
   end type simulation_type


contains


   subroutine simulation_final(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the simulation object
      implicit none
      ! Arguments
      type(simulation_type), intent(inout) :: self

      call self%deallocate()
      return
   end subroutine simulation_final


end module simulation