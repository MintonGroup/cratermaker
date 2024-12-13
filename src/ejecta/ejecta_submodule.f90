! Copyright 2023 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

submodule (ejecta) s_ejecta
    use globals
    integer(I4B), parameter :: Nraymax = 5 
        !! Maximum number of rays in a given pattern
    integer(I4B), parameter :: Npatt = 8 
        !! Number of different ray patterns
    real(DP), parameter :: frayreduction = 0.5_DP
        !! Factor by which to reduce the relative thickness of the ray for each subsequent pattern
    real(DP), parameter :: ejprofile = 3.0_DP 
        !! Ejecta profile power law index
contains

    module subroutine ejecta_distribution(radial_distance, initial_bearing, crater_diameter, ejrim, ejecta_truncation, dorays, &
                                            ejecta_thickness)
        !! author: David A. Minton
        !!
        !! Calculate the spatial distribution of ejecta in distal rays given a radial distance and initial bearing array.
        !! If dorays is true, calculate the ejecta distribution using the ray pattern. If false, use a homogenous power law profile.
        implicit none
        ! Arguments
        real(DP), dimension(:), intent(in) :: radial_distance
            !! Radial distance from the crater center in meters. 
        real(DP), dimension(:), intent(in) :: initial_bearing
            !! The initial bearing of the ray in radians.
        real(DP), intent(in) :: crater_diameter
            !! The final crater diameter in meters.
        real(DP), intent(in) :: ejrim 
            !! The final ejecta rim height in meters.
        real(DP), intent(in) :: ejecta_truncation
            !! The distance relative to the crater radius at which to truncate the ejecta pattern 
        logical(LGT), intent(in) :: dorays
            !! If true, calculate the ejecta distribution using the ray pattern. If false, use a homogenous power law profile.
        real(DP), dimension(:), intent(out) :: ejecta_thickness
            !! The thickness of the ejecta layer at each radial_distance, initial_bearing pair
        ! Internals
        real(DP) :: crater_radius
        integer(I4B) :: i, N

        N = size(radial_distance)

        if (dorays) then
            call ejecta_ray_intensity(radial_distance, initial_bearing, crater_diameter, ejecta_truncation, ejecta_thickness)
            crater_radius = crater_diameter / 2
            do concurrent (i = 1:N, radial_distance(i) >= crater_radius)
                ejecta_thickness(i) = ejecta_thickness(i) * ejecta_profile_func(radial_distance(i) / crater_radius, ejrim)
            end do
        else
            call ejecta_profile(radial_distance, crater_diameter, ejrim, ejecta_thickness)
        end if

        return
    end subroutine ejecta_distribution


    pure module subroutine ejecta_profile(radial_distance, crater_diameter, ejrim, ejecta_thickness)
        !! author: David A. Minton
        !!
        !! Calculate a homongenous power law profile of an ejecta blanket as a function of distance from the crater center.
        !! Follows the power law formula given by  McGetchin et al. (1973). 
        implicit none
        ! Arguments
        real(DP),dimension(:), intent(in) :: radial_distance  
            !! Radial distance from the crater center in meters. 
        real(DP), intent(in) :: crater_diameter 
            !! The final crater diameter in meters.
        real(DP), intent(in) :: ejrim 
            !! The final ejecta rim height in meters.
        real(DP), dimension(:), intent(out) :: ejecta_thickness 
            !! The thicknes of the ejecta relativ
        ! Internals
        real(DP) :: r, crater_radius

        ! Calculate the floor radius relative to the final crater radius
        crater_radius = crater_diameter / 2

        where(radial_distance(:) >= crater_radius)
            ejecta_thickness(:) = ejecta_profile_func(radial_distance(:) / crater_radius, ejrim)
        elsewhere
            ejecta_thickness(:) = 0.0_DP
        end where

        return         
    end subroutine ejecta_profile


    pure elemental function ejecta_profile_func(r, ejrim) result(h)
        !! author: David A. Minton
        !!
        !! Calculate the thickness of the ejecta layer as a function of distance from the center.
        !!                                
        implicit none
        ! Arguments
        real(DP), intent(in) :: r 
            !! Radial distance from the crater centerin meters. 
        real(DP), intent(in) :: ejrim
            !! The final ejecta rim height 
        real(DP) :: h
            !! The ejecta thickness at the given radial distance 

        h = ejrim * (r)**(-ejprofile)

    end function ejecta_profile_func


    module subroutine ejecta_ray_intensity(radial_distance, initial_bearing, crater_diameter, ejecta_truncation, intensity)
        !! author: David A. Minton
        !!
        !! Calculate the spatial distribution of ejecta in distal rays given a radial distance and initial bearing array.
        implicit none
        ! Arguments
        real(DP), dimension(:), intent(in) :: radial_distance
            !! Radial distance from the crater center in meters. 
        real(DP), dimension(:), intent(in) :: initial_bearing
            !! The initial bearing of the ray in radians.
        real(DP), intent(in) :: crater_diameter
            !! The final crater diameter in meters.
        real(DP), intent(in) :: ejecta_truncation
            !! The distance relative to the crater radius at which to truncate the ejecta pattern 
        real(DP), dimension(:), intent(out) :: intensity
            !! The intensity value of the ray pattern at each radial_distance, initial_bearing pair

        ! Internals
        real(DP) :: minray, rmax, rmin, theta, crater_radius
        real(DP), dimension(Nraymax) :: thetari
        real(DP), dimension(Npatt)   :: rn, ray_pattern_intensity
        real(DP) :: r_pattern
        integer(I4B) :: i, j, N

        N = size(radial_distance)
        if (size(initial_bearing) /= N) then
            print *, "Error: initial_bearing and radial_distance arrays must be the same size."
            stop
        end if
        crater_radius = crater_diameter / 2
        rmax = ejecta_truncation
        rmin = 1.0_DP 
        minray = 2.348_DP * crater_radius**(0.006_DP)  ! The continuous ejecta blanket radius relative to the crater radius

        ! Distribute ray patterns evenly around the crater
        do concurrent (i = 1:Nraymax)
            thetari(i) = 2 * pi * i / Nraymax
        end do
        call shuffle(thetari) ! randomize the ray pattern so that the ray lengths are randomly distributed
        call random_number(rn) ! randomize the ray orientation

        where(radial_distance(:) < crater_radius)
            intensity(:) = 0.0_DP
        endwhere

        do concurrent (i = 1:N, radial_distance(i) >= crater_radius)
            ! Layer the ray patterns on top of each other
            do concurrent (j = 1:Npatt)
                theta = mod(initial_bearing(i) + rn(j) * 2 * PI, 2 * PI)

                ! Make a slight random adjustment to the r values so that the indivdual rays on different patterns 
                r_pattern = radial_distance(i) / crater_radius - rn(j)
                ray_pattern_intensity(j) = frayreduction**(j-1) * ejecta_ray_intensity_func(r_pattern,theta,rmin,rmax,thetari, &
                                                                                                minray)
            end do
            intensity(i) = sum(ray_pattern_intensity(:)) 
        end do

        ! Normalize the intensity values so they span from 0 to 1
        where(radial_distance(:) >= crater_radius)
            intensity(:) = intensity(:) / maxval(intensity(:))
        end where

        contains
            subroutine shuffle(a)
                real(DP), intent(inout) :: a(:)
                integer :: i, randpos
                real(DP) :: r,temp
             
                do i = size(a), 2, -1
                    call random_number(r)
                    randpos = int(r * i) + 1
                    temp = a(randpos)
                    a(randpos) = a(i)
                    a(i) = temp
                end do
                return
           end subroutine shuffle

    end subroutine ejecta_ray_intensity


    pure function ejecta_ray_intensity_func(r,theta,rmin,rmax,thetari,minray) result(ans)
        !! author: David A. Minton
        !!
        !! Calculate the spatial distribution of ejecta in distal rays.
        implicit none
        ! Arguments
        real(DP),intent(in) :: r,theta,rmin,rmax
        real(DP),dimension(:),intent(in) :: thetari
        real(DP),intent(in) :: minray
        ! Result
        real(DP) :: ans
        ! Internals
        real(DP) :: thetar,rw, w
        real(DP) :: length,minray,FF
        integer(I4B) :: n,i
        real(DP) :: tmp
        real(DP), parameter :: rayp = 4.0_DP ! Factor that controls how the lengths of rays are distributed
     
        if (r > rmax) then
            ans = 0._DP
        else if (r < 1.0_DP) then
            ans = 1.0_DP
        else
            tmp = (Nraymax**rayp - (Nraymax**rayp - 1) * log(r/minray) / log(rmax/minray))
            if (tmp < 0.0_DP) then
               n = Nraymax ! "Nrays" in Minton et al. (2019)
            else ! Exponential decay of ray number with distance
               n = max(min(floor((Nraymax**rayp - (Nraymax**rayp - 1) * log(r/minray) / log(rmax/minray))**(1._DP/rayp)),Nraymax),1)
            end if
            ans = 0._DP
            do i = 1,Nraymax
                length = minray * exp(log(rmax/minray) * ((Nraymax - i + 1)**rayp - 1.0_DP) / ((Nraymax**rayp - 1)))
                if (r > length) cycle ! Don't add any material beyond the length of the ray
                w = (rmax / length) ** (1.0_DP)
                ! equation 41 Minton et al. 2019
                rw = PI / (w * Nraymax) * (rmin / r) * (1._DP - (1.0_DP - w / rmin) * exp(1._DP - (r / rmin)**2)) 
                tmp = ejecta_ray_func(theta,thetari(i),r,n,rw) 
                if (tmp > epsilon(ans)) ans = ans + tmp  ! Ensure that we don't get an underflow
           end do
        end if

        return
    end function ejecta_ray_intensity_func 

    
    pure function ejecta_ray_func(theta,thetar,r,n,w) result(ans)
        !! author: David A. Minton
        !!
        !! Calculate the spatial distribution of ejecta in a single ray. The intensity function is a Gaussian centered on the 
        !! midline of the ray.
        implicit none
        real(DP) :: ans
        real(DP),intent(in) :: theta,thetar,r,w
        integer(I4B),intent(in) :: n
        real(DP) :: thetap,thetapp,a,b,c,dtheta,logval
    
        c = w / r
        b = thetar 
        dtheta = min(2*pi - abs(theta - b),abs(theta - b))
        logval = -dtheta**2 / (2 * c**2)
        if (logval < LOGVSMALL) then
            ans = 0.0_DP
        else
            a = sqrt(2 * pi) / (n * c * erf(pi / (2 *sqrt(2._DP) * c))) 
            ans = a * exp(logval)
        end if
    
        return
    end function ejecta_ray_func
end submodule s_ejecta