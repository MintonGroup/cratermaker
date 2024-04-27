! Copyright 2023 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

!! The implementations for the Perlin noise procedures.
submodule (ejecta) s_ejecta
    use globals
    integer(I4B), parameter :: Nraymax = 11
        !! Maximum number of rays in a given pattern
    real(DP), parameter :: rayfmult = (Nraymax)**(-4.0_DP / (1.2_DP)) 
        !! The factor by which to multiply the ray intensity
    real(DP), parameter :: fpeak = 8000.0_DP 
        !! Factor that sets the peak intensity of a ray 
    real(DP), parameter :: rayp = 2.0_DP 
        !! Another ray intensity factor
    integer(I4B), parameter :: Npatt = 8 
        !! Number of different ray patterns
    real(DP), parameter :: frayreduction = 0.5_DP
        !! Factor by which to reduce the relative thickness of the ray for each subsequent pattern
    real(DP), parameter :: ejprofile = 3.0_DP 
        !! Ejecta profile power law index
contains

    pure module subroutine ejecta_profile(radial_distance, crater_diameter, ejrim, ejecta_thickness)
        !! author: David A. Minton
        !!
        !! Calculate the thickness of an ejecta blanket as a function of distance from the crater center.
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


    module subroutine ejecta_ray_pattern(radial_distance, initial_bearing, crater_diameter, ejrim, ejecta_truncation, &
                                        ejecta_thickness)
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
        real(DP), intent(in) :: ejrim 
            !! The final ejecta rim height in meters.
        real(DP), intent(in) :: ejecta_truncation
            !! The distance relative to the crater radius at which to truncate the ejecta pattern 
        real(DP), dimension(:), intent(out) :: ejecta_thickness
            !! The thickness of the ejecta layer at each radial_distance, initial_bearing pair

        ! Internals
        real(DP) :: l1, rmax, rmin, theta, crater_radius
        real(DP), dimension(Nraymax) :: thetari
        real(DP), dimension(Npatt)   :: rn, ray_pattern_thickness
        integer(I4B) :: i, j, N

        N = size(radial_distance)
        if (size(initial_bearing) /= N) then
            print *, "Error: initial_bearing and radial_distance arrays must be the same size."
            stop
        end if
        crater_radius = crater_diameter / 2
        rmax = ejecta_truncation
        rmin = 2.348_DP * crater_radius**(0.006_DP)  ! The continuous ejecta blanket thickness relative to the crater radius
        l1 = (5.32_DP*(crater_radius/1000)**1.27)/(crater_radius/1000)

        do concurrent (i = 1:Nraymax)
            thetari(i) = 2 * pi * i / Nraymax
        end do
        call shuffle(thetari) ! randomize the ray pattern
        call random_number(rn) ! randomize the ray orientation

        where(radial_distance(:) < crater_radius)
            ejecta_thickness(:) = 0.0_DP
        endwhere

        do concurrent (i = 1:N, radial_distance(i) >= crater_radius)
            do concurrent (j = 1:Npatt)
                theta = mod(initial_bearing(i) + rn(j) * 2 * PI, 2 * PI)
                ray_pattern_thickness(j) = frayreduction**(j-1) * ejecta_ray_pattern_func(radial_distance(i) / crater_radius,&
                                                                                            theta,rmin,rmax,thetari,l1)
            end do
            ejecta_thickness(i) = sum(ray_pattern_thickness(:))
            ejecta_thickness(i) = ejecta_thickness(i) * ejecta_profile_func(radial_distance(i) / crater_radius, ejrim)
        end do

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

    end subroutine ejecta_ray_pattern


    pure function ejecta_ray_pattern_func(r,theta,rmin,rmax,thetari,l1) result(ans)
        !! author: David A. Minton
        !!
        !! Calculate the spatial distribution of ejecta in distal rays.
        implicit none
        ! Arguments
        real(DP),intent(in) :: r,theta,rmin,rmax
        real(DP),dimension(:),intent(in) :: thetari
        real(DP),intent(in) :: l1
        ! Result
        real(DP) :: ans
        ! Internals
        real(DP) :: a,c
        real(DP) :: thetar,rw,rw0,rw1
        real(DP) :: rtrans,length,minray,FF
        integer(I4B) :: n,i
        real(DP) :: tmp
        real(DP) :: width_factor
     
        minray = l1
     
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
            rtrans = r - 1.0_DP
            rw1 = 2 * pi / Nraymax
            do i = 1,Nraymax
                length = minray * exp(log(rmax/minray) * ((Nraymax - i + 1)**rayp - 1.0_DP) / ((Nraymax**rayp - 1)))
                width_factor = log(rmax / (0.99_DP * length)) / log(2 * rmax/rmin)
                rw0 = (rmin * pi / Nraymax) * width_factor
                rw = rw0 * (1._DP - (1.0_DP - rw1 / rw0) * exp(1._DP - (r / rmin)**2)) ! equation 40 Minton et al. 2019
                c = rw / r
                a = sqrt(2 * pi) / (n * c * erf(pi / (2 *sqrt(2._DP) * c))) !equation 39 Minton et al., 2019
                if (r > length) cycle ! Don't add any material beyond the length of the ray
                tmp = ejecta_ray_func(theta,thetari(i),a,n,rw)
                if (tmp > epsilon(ans)) ans = ans + tmp  ! Ensure that we don't get an underflow
           end do
        end if

        return
    end function ejecta_ray_pattern_func    


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