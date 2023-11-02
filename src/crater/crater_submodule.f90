!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 


submodule (crater) s_crater
contains

   module subroutine crater_generate()
   implicit none
   end subroutine crater_generate

! subroutine crater_generate(user,crater,domain,prod,production_list,vdist,surf)
!    use globals
!    implicit none

!    ! Arguments
!    type(usertype),intent(in) :: user
!    type(cratertype),intent(inout) :: crater
!    type(domaintype),intent(in) :: domain
!    real(DP),dimension(:,:),intent(in),optional :: prod,vdist
!    integer(I8B),dimension(:),intent(inout),optional :: production_list
!    type(surftype),dimension(:,:),intent(in),optional :: surf

!    ! Internal variables
!    real(DP),dimension(6)    :: rn        
!    real(DP)                 :: nmark,frac,limp
!    real(DP)                 :: lnmark,lprod1,lprod1p,lprod2,lprod2p
!    real(DP)                 :: dburial,trfin
!    integer(I4B)             :: k,khi,klo,Nk
!    integer(I8B)             :: random_index,numremaining,nabove,nabovep1

!    ! Get all six random numbers we need in one call
!    if (.not.domain%initialize) then
!       ! Initialize the random number generator with the current value of the seeds
!       call random_seed(put=crater%seedarr)
!       call random_number(rn)
!       ! Save the current value of the seeds for the next time we need a new crater. This ensures we can get a repeatable population population of craters, 
!       ! regardless of whether other procedures that use the random number generator are used or not
!       call random_seed(get=crater%seedarr)
!    end if

!    ! Find crater center position
!    if (domain%initialize) then
!       crater%xl = real(domain%side * 0.5_DP, kind=SP)
!       crater%yl = real(domain%side * 0.5_DP, kind=SP)
!    else if (user%testflag) then
!       crater%xl = real(domain%side * 0.5_DP + user%testxoffset, kind=SP)
!       crater%yl = real(domain%side * 0.5_DP + user%testyoffset, kind=SP)
!    else
!       crater%xl = real(domain%side * rn(1), kind=SP)
!       crater%yl = real(domain%side * rn(2), kind=SP)
!    end if

!    crater%xlpx = min(max(nint(crater%xl / user%pix),1),user%gridsize)
!    crater%ylpx = min(max(nint(crater%yl / user%pix),1),user%gridsize)

!    ! Get impactor size from the distribution
!    ! Reverted back to old code because the new one was way too slow
!    if (.not.domain%initialize) then
!       if (user%testflag) then
!          crater%imp = user%testimp
!       else 
!          if (domain%pnum == 1) then
!             crater%imp = prod(1,1)
!          else ! Draw a random impactor from the production SFD
!             numremaining = sum(production_list)
!             random_index = max(int(numremaining * rn(3)),1)
!             khi = domain%pnum
!             klo = domain%smallest_impactor_index
!             k = klo + (khi - klo) / 2
!             do
!                nabove = sum(production_list(k:domain%pnum))
!                nabovep1 = sum(production_list(k+1:domain%pnum))
!                if ((random_index <= nabove).and.(random_index > nabovep1)) exit
!                if (random_index > nabove) then
!                   khi = k
!                   k = klo + (khi - klo) / 2
!                else  
!                   klo = k
!                   k = klo + (khi - klo) / 2
!                end if
!                if (klo == khi) then
!                   klo = klo - 1
!                   khi = khi + 1
!                end if 
!                if (k > khi) then
!                   write(*,*)
!                   write(*,*) 'Error in crater_generate: Went past the end of the SFD'
!                   write(*,*) random_index
!                   stop
!                   exit
!                end if
!             end do
             
!             production_list(k) = production_list(k) - 1
!             frac = rn(4)
!             lprod1 = prod(3,k) 
!             lprod1p = prod(3,k + 1) 
!             lprod2 = prod(4,k) 
!             lprod2p = prod(4,k + 1) 
!             limp = lprod1 + (frac * (lprod1p - lprod1))
!             crater%imp = exp(limp)
!          end if
!       end if
!       crater%imp = crater%imp * (1._DP + 1.0e-3_DP*rn(3)) ! Some user-input SFDs can result in many craters having identical 
!                                                           ! diameters. This random number prevents more than one crater from having 
!                                                           ! exactly the same diameter, as diameter is used as identification.
!    end if
                                                       

!    crater%impmass = 4*THIRD*PI*user%prho*(crater%imp*0.5_DP)**3

!    ! Determine which strength model to use


!    !  find impact angle
!    if (.not.domain%initialize) then
!       if (user%testflag) then
!          crater%sinimpang = sin(user%testang * DEG2RAD)
!       else
!          if (user%doangle) then
!             crater%sinimpang = sqrt(rn(5))
!          else
!             crater%sinimpang = 1._DP ! Vertical impact only
!          end if
!       end if
!    end if

!    if (.not.domain%initialize) then
!       if (user%testflag) then
!          crater%impvel = user%testvel
!       else
!          if (domain%vnum == 1) then
!             crater%impvel = vdist(1,1)
!          else 
!             !  Draw impact velocity from the velocity distribution
!             nmark = vdist(3,domain%vlo) + (rn(6) * (vdist(3,domain%vhi) - vdist(3,domain%vlo)))
!             ! Make a guess as to where in the velocity distribution the impactor might be. 
!             ! This could speed up the searching if the velocity distribution has a lot of elements in it.
!             klo = domain%vlo
!             khi = domain%vhi
!             Nk = 1 + khi - klo
!             k = domain%vlo + int(Nk * rn(6) * (vdist(3,khi) - vdist(3,klo)))
!             call util_search(vdist,3,domain%vnum-1,nmark,k)
!             if (k == 0) then
!                crater%impvel = vdist(1,1)
!             else if (k == domain%vnum) then
!                crater%impvel = vdist(1,domain%vnum)
!             else
!                frac = (nmark-vdist(3,k))/(vdist(3,k+1)-vdist(3,k))
!                crater%impvel = vdist(1,k) + (frac*(vdist(1,k+1)-vdist(1,k)))
!             end if
!          end if
!       end if
!    end if

!    !  scale to crater size
!    if (.not.domain%initialize) crater%strflag = 0 ! Begin with regolith strength
!    call crater_scale(user,crater%imp,crater%rad,crater%grad,crater%strflag,crater%sinimpang,crater%impvel)
!    if (.not.domain%initialize) then 
!       dburial = EXFAC * crater%rad  
!       if (dburial > surf(crater%xlpx,crater%ylpx)%ejcov) then 
!          crater%strflag = 1 ! Use bedrock strength 
!          call crater_scale(user,crater%imp,crater%rad,crater%grad,crater%strflag,crater%sinimpang,crater%impvel) 
!       end if 
!    end if 


!    trfin = 2 * TRSIM * crater%rad
!    if (trfin <= crater%cxtran) then
!       crater%morphtype = "SIMPLE"
!       crater%fcrat = trfin   ! Simple Crater
!    else  ! Complex crater
!       crater%morphtype = "COMPLEX"
!       crater%fcrat = trfin * ((trfin/crater%cxtran)**crater%cxexp) ! Complex Crater
!    end if
!    if (crater%imp >= user%basinimp) then
!       crater%morphtype = "MULTIRING" 
!       ! This section is temporary until a better basin scaling law can be implemented. Potter et al. (2012) GRL v39. p 18203
!       if ((crater%xl < (0.5_DP * domain%side) ).or.(user%testflag).or.(domain%initialize)) then ! pick TP1 for left-hand hemisphere
!          crater%fcrat = 0.1354e3_DP * (0.5_DP * trfin * 1e-3_DP)**(1.389_DP)  ! TP2
!          crater%fcrat = 2 * 1.56_DP * crater%fcrat
!       else
!          crater%fcrat = 0.0718e3_DP * (0.5_DP * trfin * 1e-3_DP)**(1.613_DP)  ! TP1
!          crater%fcrat = 2 * 1.2_DP * crater%fcrat
!       end if
!    end if

!    crater%frad = 0.5_DP  * crater%fcrat
   
!    ! Get pixel space values
!    crater%fcratpx = nint(crater%fcrat / user%pix)
!    crater%fe = user%fe
!    return
! end subroutine crater_generate

end submodule s_crater
