! Copyright 2023 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

!! The implementations for the Perlin noise procedures.
submodule (perlin_module) s_perlin
    use globals_module
contains
   pure module elemental subroutine perlin_noise(xx,yy,zz,noise,dx,dy,dz)
      ! Perlin noise with derivatives. Adapted from Ken Perlin's original code, with derivatives
      ! that are used in the noise functions by Giliam de Carpentier
      implicit none
      real(DP),intent(in) :: xx,yy,zz
      real(DP),intent(out) :: noise
      real(DP),intent(out),optional :: dx, dy, dz
   
      ! Internal variables
      integer(I4B) :: xi,yi,zi,A,B,AA,BB,AB,BA
      real(DP) :: x,y,z
   
      ! Intermediate variables used to calculate noise
      real(DP) :: u,v,w
      real(DP) :: GA1,GB1,GA2,GB2,GA3,GB3,GA4,GB4
      real(DP) :: L1,L2,L3,L4,LL1,LL2
   
      ! Intermediate variables used to calculate the derivative wrt x
      real(DP) :: ux,vx,wx
      real(DP) :: GA1x,GB1x,GA2x,GB2x,GA3x,GB3x,GA4x,GB4x
      real(DP) :: L1x,L2x,L3x,L4x,LL1x,LL2x
   
      ! Intermediate variables used to calculate the derivative wrt y
      real(DP) :: uy,vy,wy
      real(DP) :: GA1y,GB1y,GA2y,GB2y,GA3y,GB3y,GA4y,GB4y
      real(DP) :: L1y,L2y,L3y,L4y,LL1y,LL2y
   
      ! Intermediate variables used to calculate the derivative wrt z
      real(DP) :: uz,vz,wz
      real(DP) :: GA1z,GB1z,GA2z,GB2z,GA3z,GB3z,GA4z,GB4z
      real(DP) :: L1z,L2z,L3z,L4z,LL1z,LL2z
   
      real(DP),parameter :: hcorrection = 1.25_DP 
   
      integer(I4B),dimension(0:511),parameter :: p = (/ 151,160,137,91,90,15, &
      131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,&
      190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,&
      88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,&
      77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,&
      102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,&
      135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,&
      5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,&
      223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,&
      129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,&
      251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,&
      49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,&
      138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,&
      151,160,137,91,90,15,&
      131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,&
      190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,&
      88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,&
      77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,&
      102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,&
      135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,&
      5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,&
      223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,&
      129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,&
      251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,&
      49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,&
      138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180 /) 
   
      x = xx
      y = yy
      z = zz
      xi = iand(floor(x),255)
      yi = iand(floor(y),255)
      zi = iand(floor(z),255)
      x = x - real(floor(x), kind=DP)
      y = y - real(floor(y), kind=DP)
      z = z - real(floor(z), kind=DP)
      u = fade(x)
      v = fade(y)
      w = fade(z)
   
      A = p(xi) + yi
      B = p(xi + 1) + yi
   
      AA = p(A) + zi
      BA = p(B) + zi
   
      AB = p(A + 1) + zi
      BB = p(B + 1) + zi
   
      GA1 = grad(p(AA    ), x        , y        , z)
      GB1 = grad(p(BA    ), x - 1._DP, y        , z)
      GA2 = grad(p(AB    ), x        , y - 1._DP, z)
      GB2 = grad(p(BB    ), x - 1._DP, y - 1._DP, z)
      GA3 = grad(p(AA + 1), x        , y        , z - 1._DP)
      GB3 = grad(p(BA + 1), x - 1._DP, y        , z - 1._DP)
      GA4 = grad(p(AB + 1), x        , y - 1._DP, z - 1._DP)
      GB4 = grad(p(BB + 1), x - 1._DP, y - 1._DP, z - 1._DP)

   
      L1 = lerp(u, GA1, GB1)
      L2 = lerp(u, GA2, GB2)
      L3 = lerp(u, GA3, GB3)
      L4 = lerp(u, GA4, GB4)

      LL1 = lerp(v, L1, L2)
      LL2 = lerp(v, L3, L4)
     
      noise = hcorrection * lerp(w, LL1, LL2)

      if (present(dx)) then
         ! Calculate derivatives wrt x
         ux = dfade(x)
         vx = 0.0_DP
         wx = 0.0_DP
         GA1x = grad(p(AA    ), 1._DP, 0._DP, 0._DP)
         GB1x = grad(p(BA    ), 1._DP, 0._DP, 0._DP)
         GA2x = grad(p(AB    ), 1._DP, 0._DP, 0._DP)
         GB2x = grad(p(BB    ), 1._DP, 0._DP, 0._DP)
         GA3x = grad(p(AA + 1), 1._DP, 0._DP, 0._DP)
         GB3x = grad(p(BA + 1), 1._DP, 0._DP, 0._DP)
         GA4x = grad(p(AB + 1), 1._DP, 0._DP, 0._DP)
         GB4x = grad(p(BB + 1), 1._DP, 0._DP, 0._DP)
   
         L1x = dlerp(u, GA1, GB1, ux, GA1x, GB1x)
         L2x = dlerp(u, GA2, GB2, ux, GA2x, GB2x)
         L3x = dlerp(u, GA3, GB3, ux, GA3x, GB3x)
         L4x = dlerp(u, GA4, GB4, ux, GA4x, GB4x)
    
         LL1x = dlerp(v, L1, L2, vx, L1x, L2x)
         LL2x = dlerp(v, L3, L4, vx, L3x, L4x)
    
         dx = hcorrection * dlerp(w, LL1, LL2, wx, LL1x, LL2x)
      end if

      if (present(dy)) then
         ! Calculate derivatives wrt y
         uy = 0.0_DP
         vy = dfade(y)
         wy = 0.0_DP
         GA1y = grad(p(AA    ), 0._DP, 1._DP, 0._DP)
         GB1y = grad(p(BA    ), 0._DP, 1._DP, 0._DP)
         GA2y = grad(p(AB    ), 0._DP, 1._DP, 0._DP)
         GB2y = grad(p(BB    ), 0._DP, 1._DP, 0._DP)
         GA3y = grad(p(AA + 1), 0._DP, 1._DP, 0._DP)
         GB3y = grad(p(BA + 1), 0._DP, 1._DP, 0._DP)
         GA4y = grad(p(AB + 1), 0._DP, 1._DP, 0._DP)
         GB4y = grad(p(BB + 1), 0._DP, 1._DP, 0._DP)
   
         L1y = dlerp(u, GA1, GB1, uy, GA1y, GB1y)
         L2y = dlerp(u, GA2, GB2, uy, GA2y, GB2y)
         L3y = dlerp(u, GA3, GB3, uy, GA3y, GB3y)
         L4y = dlerp(u, GA4, GB4, uy, GA4y, GB4y)
    
         LL1y = dlerp(v, L1, L2, vy, L1y, L2y)
         LL2y = dlerp(v, L3, L4, vy, L3y, L4y)
    
         dy = hcorrection * dlerp(w, LL1, LL2, wy, LL1y, LL2y)
      end if

      if (present(dz)) then
         ! Calculate derivatives wrt z
         uz = 0.0_DP
         vz = 0.0_DP
         wz = dfade(z)
         GA1z = grad(p(AA    ), 0._DP, 0._DP, 1._DP)
         GB1z = grad(p(BA    ), 0._DP, 0._DP, 1._DP)
         GA2z = grad(p(AB    ), 0._DP, 0._DP, 1._DP)
         GB2z = grad(p(BB    ), 0._DP, 0._DP, 1._DP)
         GA3z = grad(p(AA + 1), 0._DP, 0._DP, 1._DP)
         GB3z = grad(p(BA + 1), 0._DP, 0._DP, 1._DP)
         GA4z = grad(p(AB + 1), 0._DP, 0._DP, 1._DP)
         GB4z = grad(p(BB + 1), 0._DP, 0._DP, 1._DP)
   
         L1z = dlerp(u, GA1, GB1, uz, GA1z, GB1z)
         L2z = dlerp(u, GA2, GB2, uz, GA2z, GB2z)
         L3z = dlerp(u, GA3, GB3, uz, GA3z, GB3z)
         L4z = dlerp(u, GA4, GB4, uz, GA4z, GB4z)
    
         LL1z = dlerp(v, L1, L2, vz, L1z, L2z)
         LL2z = dlerp(v, L3, L4, vz, L3z, L4z)
    
         dz = hcorrection * dlerp(w, LL1, LL2, wz, LL1z, LL2z)
      end if
   
      return

      contains
      
         pure function fade(t) result(f)
            implicit none
            real(DP),intent(in) :: t
            real(DP) :: f
      
            f = t * t * t * (t * (t * 6._DP - 15._DP) + 10._DP)
            return
         end function fade
      
         pure function dfade(t) result(df)
            implicit none
            real(DP),intent(in) :: t
            real(DP) :: df
      
            df = t * t * (t * (t * 30._DP - 60._DP) + 30._DP)
            return
         end function dfade
      
         pure function lerp(t,a,b) result(l)
            implicit none
            real(DP),intent(in) :: t,a,b
            real(DP) :: l
      
            l = a + t * (b - a)
            return
         end function lerp
   
         pure function dlerp(t,a,b,dt,da,db) result(dl)
            implicit none
            real(DP),intent(in) :: t,a,b
            real(DP),intent(in) :: dt,da,db
            real(DP) :: dl
      
            dl = da + t * (db - da) + dt * (b - a)

            return
         end function dlerp
   
         pure function grad(hash,x,y,z) result(g)
            implicit none
            integer(I4B),intent(in) :: hash
            real(DP),intent(in) :: x,y,z
            real(DP) :: g,u,v
            integer(I4B) :: h
      
            h = iand(hash,15) 
            if (h < 8) then 
               u = x
            else
               u = y
            end if
            if (h < 4) then
               v = y
            else if ((h == 12) .or. (h == 14)) then
               v = x
            else
               v = z
            end if
            if (iand(h,1 ) == 0) then
               g = u
            else 
               g = -u
            end if
            if (iand(h, 2) == 0) then
               g = g + v
            else
               g = g - v
            end if
            return
         end function grad
   
   end subroutine perlin_noise


   pure module function perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, z, noise_height, freq, pers
      integer(I4B), intent(in) :: num_octaves
      real(DP), dimension(:,:),intent(in) :: anchor
      real(DP) :: noise,xnew,ynew,znew,norm
      
      integer(I4B) :: i
      real(DP) :: spatial_fac, noise_mag,dn
      
      noise = 0.0_DP
      norm = 0.5_DP
      do i = 1, num_octaves
         spatial_fac = freq ** (i - 1)
         noise_mag = pers ** (i - 1)
         norm = norm + 0.5_DP * noise_mag
         xnew = (x + anchor(1,i)) * spatial_fac
         ynew = (y + anchor(2,i)) * spatial_fac
         znew = (z + anchor(3,i)) * spatial_fac
         call perlin_noise(xnew,ynew,znew,dn)
         noise = noise + dn * noise_mag
      end do
      
      noise = noise * noise_height / norm
      
      return
      
   end function perlin_turbulence


   pure module function perlin_billowedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, z, noise_height, freq, pers
      integer(I4B), intent(in) :: num_octaves
      real(DP),dimension(:,:),intent(in) :: anchor
      real(DP) :: noise
      
      noise = abs(perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor))
      
      return
   end function perlin_billowedNoise
      
      
   pure module function perlin_plawNoise(x, y, z, noise_height, freq, pers, slope, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, z, noise_height, freq, pers,slope
      integer(I4B), intent(in) :: num_octaves
      real(DP),dimension(:,:),intent(in) :: anchor
      real(DP) :: noise
      
      noise = perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor)
      noise = sign(abs(noise)**(slope),noise)
      
      return
   end function perlin_plawNoise
      
      
   pure module function perlin_ridgedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, z, noise_height, freq, pers
      integer(I4B), intent(in) :: num_octaves
      real(DP),dimension(:,:),intent(in) :: anchor
      real(DP) :: noise
      
      noise = noise_height - abs(perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor))
      
      return
   end function perlin_ridgedNoise
      
      
   pure module function perlin_swissTurbulence(x, y, z, lacunarity, gain, warp, num_octaves, anchor) result(noise)
      implicit none
      ! Arguments
      real(DP), intent(in) ::  x, y, z, lacunarity, gain, warp
      integer(I4B), intent(in) :: num_octaves
      real(DP),dimension(:,:), intent(in) :: anchor
      ! Result
      real(DP) :: noise
      ! Internals
      real(DP) :: freq, amp,xnew,ynew,znew, norm
      real(DP),dimension(3) :: dsum
      real(DP),dimension(4) :: n
      integer(I4B) :: i
   
      noise =  0.0_DP
      freq = 1.0_DP
      amp = 1.0_DP
      dsum(:) = 0._DP
   
      norm = 0.5_DP
         
      do i = 1, num_octaves 
         xnew = (x + anchor(1,i) + warp * dsum(1)) * freq
         ynew = (y + anchor(2,i) + warp * dsum(2)) * freq
         znew = (z + anchor(3,i) + warp * dsum(3)) * freq
         call perlin_noise(xnew,ynew,znew,n(1),dx = n(2),dy = n(3), dz = n(4))
         noise = noise + amp * (1._DP - abs(n(1)))
         norm = norm + amp
         dsum(:) = dsum(:) + amp * n(2:4) * (-n(1))
         freq = freq * lacunarity;
         amp = amp * gain * max(min(noise,1.0_DP),0.0_DP)
      end do
      noise = noise  
      return
   
   end function perlin_swissTurbulence
   
      
   pure module function perlin_jordanTurbulence(x, y, z, lacunarity, gain0, gain, warp0, warp, damp0, damp, damp_scale,&
            num_octaves, anchor) result(noise)
      ! Fortran implementation of noise function by Giliam de Carpentier
      implicit none
      real(DP),intent(in) :: x, y, z, lacunarity, gain0, gain, warp0, warp, damp0, damp, damp_scale
      integer(I4B),intent(in) :: num_octaves
      real(DP), dimension(:,:), intent(in) :: anchor
      ! Result variable
      real(DP) :: noise
   
      ! Internal variables
      real(DP),dimension(4) :: n,n2
      real(DP),dimension(3) :: dsum_warp,dsum_damp
      real(DP) :: freq,damped_amp,xnew,ynew,znew, amp
      integer(I4B) :: i
   
      xnew = x + anchor(1,1)
      ynew = y + anchor(2,1)
      znew = z + anchor(3,1)
      call perlin_noise(xnew,ynew,znew, n(1),dx = n(2),dy = n(3), dz = n(4))
      n2(:) = n(:) * n(1)
      noise = n2(1)
      dsum_warp(:) = warp0 * n2(2:4)
      dsum_damp(:) = damp0 * n2(2:4)
   
      freq = lacunarity
      amp = gain0
      damped_amp = amp * gain
   
      do i = 2, num_octaves
         xnew = x * freq + dsum_warp(1) + anchor(1,i)
         ynew = y * freq + dsum_warp(2) + anchor(2,i)
         znew = z * freq + dsum_warp(3) + anchor(3,i)
         call perlin_noise(xnew,ynew,znew,n(1),dx = n(2),dy = n(3),dz = n(4))
         n2(:) = n(:) * n(1)
         noise = noise + damped_amp * n2(1)
         dsum_warp(:) = dsum_warp(:) + warp * n2(2:4)
         dsum_damp(:) = dsum_damp(:) + damp * n2(2:4)
         freq = freq * lacunarity
         amp = amp * gain
         damped_amp = amp * (1._DP - damp_scale / (1._DP + dot_product(dsum_damp(:),dsum_damp(:))))
      end do
      return 
   end function perlin_jordanTurbulence


   pure module function perlin_noise_one(model, x, y, z, num_octaves, anchor, damp, damp0, damp_scale, freq, gain, gain0,&
                                          lacunarity, noise_height, pers, slope, warp, warp0) result(noise)
      implicit none
      character(len=*), intent(in) :: model !! The specific turbulence model to apply
      real(DP), intent(in) :: x, y, z  !! The xyz cartesian position of the noise to evaluate
      integer(I4B), intent(in) :: num_octaves
      real(DP), dimension(:,:), intent(in) :: anchor
      real(DP), intent(in) :: damp, damp0, damp_scale, freq, gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0
      real(DP) :: noise

      select case (trim(model))
      case('turbulence')
         noise = perlin_turbulence(x, y, z, noise_height, freq, pers, num_octaves, anchor) 
      case('billowed')
         noise = perlin_billowedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor)
      case('plaw')
         noise = perlin_plawNoise(x, y, z, noise_height, freq, pers, slope, num_octaves, anchor)
      case('ridged')
         noise = perlin_ridgedNoise(x, y, z, noise_height, freq, pers, num_octaves, anchor)
      case('swiss')
         noise = perlin_swissTurbulence(x, y, z, lacunarity, gain, warp, num_octaves, anchor)
      case('jordan')
         noise = perlin_jordanTurbulence(x, y, z, lacunarity, gain0, gain, warp0, warp, damp0, damp,& 
                                                damp_scale, num_octaves, anchor) 
      case default
         noise = 0.0_DP
      end select
      return
   end function perlin_noise_one


   module subroutine perlin_noise_all(model, x, y, z, num_octaves, anchor, damp, damp0, damp_scale, freq, gain, gain0,&
                                          lacunarity, noise_height, pers, slope, warp, warp0, noise)
      implicit none
      character(len=*), intent(in) :: model !! The specific turbulence model to apply
      real(DP), dimension(:), intent(in) :: x, y, z  !! The xyz cartesian position of the noise to evaluate
      integer(I4B), intent(in) :: num_octaves
      real(DP), dimension(:,:), intent(in) :: anchor
      real(DP), intent(in) :: damp, damp0, damp_scale, freq, gain, gain0, lacunarity, noise_height, pers, slope, warp, warp0
      real(DP), dimension(:), intent(out) :: noise
      integer :: i, N

      N = size(x)

      select case (trim(model))
      case('turbulence')
         do concurrent(i = 1:N)
            noise(i) = perlin_turbulence(x(i), y(i), z(i), noise_height, freq, pers, num_octaves, anchor) 
         end do
      case('billowed')
         do concurrent(i = 1:N)
            noise(i) = perlin_billowedNoise(x(i), y(i), z(i), noise_height, freq, pers, num_octaves, anchor)
         end do
      case('plaw')
         do concurrent(i = 1:N)
            noise(i) = perlin_plawNoise(x(i), y(i), z(i), noise_height, freq, pers, slope, num_octaves, anchor)
         end do
      case('ridged')
         do concurrent(i = 1:N)
            noise(i) = perlin_ridgedNoise(x(i), y(i), z(i), noise_height, freq, pers, num_octaves, anchor)
         end do
      case('swiss')
         do concurrent(i = 1:N)
            noise(i) = perlin_swissTurbulence(x(i), y(i), z(i), lacunarity, gain, warp, num_octaves, anchor)
         end do
      case('jordan')
         do concurrent(i = 1:N)
            noise(i) = perlin_jordanTurbulence(x(i), y(i), z(i), lacunarity, gain0, gain, warp0, warp, damp0, damp,& 
                                                damp_scale, num_octaves, anchor) 
         end do
      case default
         noise(:) = 0.0_DP
      end select
      return
   end subroutine perlin_noise_all

end submodule s_perlin