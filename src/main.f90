!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

program main
   use driver
   !! author: David A. Minton
   !!
   !! This is used as a wrapper for the cratermaker_driver subroutine so that the driver can be used as either
   !! a function call or a standalone executable program.
   implicit none
   call cratermaker_driver()
end program