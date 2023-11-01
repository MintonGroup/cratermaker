!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module crater

   interface
      module subroutine crater_generate(user,crater,domain,prod,production_list,vdist,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in)    :: domain
      real(DP),dimension(:,:),intent(in),optional :: prod,vdist
      integer(I8B),dimension(:),intent(inout),optional  :: production_list            
      type(surftype),dimension(:,:),intent(in),optional :: surf
      end subroutine crater_generate
   end interface

end module crater