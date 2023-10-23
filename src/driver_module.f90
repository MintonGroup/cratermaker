!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

module driver
   contains

   module subroutine cratermaker_driver()

      !! author: David A. Minton
      !!
      !! Driver program for the cratermaker integrators. Unlike the earlier Swift and Swifter drivers, in cratermaker all integrators 
      !!    are run from this single program. 
      !!
      !! Adapted from Swifter by David E. Kaufmann's Swifter driver programs swifter_[bs,helio,ra15,rmvs,symba,tu4,whm].f90
      !! Adapted from Hal Levison and Martin Duncan's Swift driver programs
      implicit none

      write(*,*) "Makin' craters!"
      return
   end subroutine cratermaker_driver

end module driver
