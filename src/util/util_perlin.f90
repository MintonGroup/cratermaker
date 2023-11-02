!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (util) s_util_perlin
   use globals
contains

   module subroutine util_perlin_noise(xx,yy,noise, dx,dy)
      ! Perlin noise with derivatives. Adapted from Ken Perlin's original code, with derivatives
      ! that are used in the noise functions by Giliam de Carpentier
      implicit none
      real(DP),intent(in) :: xx,yy
      real(DP),intent(out) :: noise
      real(DP),intent(out),optional :: dx, dy
   
      ! Internal variables
      integer(I4B) :: xi,yi,A,B,AA,BB,AB,BA
      real(DP) :: x,y
   
      ! Intermediate variables used to calculate noise
      real(DP) :: u,v
      real(DP) :: GA1,GB1,GA2,GB2
      real(DP) :: L1,L2
   
      ! Intermediate variables used to calculate the derivative wrt x
      real(DP) :: ux,vx
      real(DP) :: GA1x,GB1x,GA2x,GB2x
      real(DP) :: L1x,L2x
   
      ! Intermediate variables used to calculate the derivative wrt y
      real(DP) :: uy,vy
      real(DP) :: GA1y,GB1y,GA2y,GB2y
      real(DP) :: L1y,L2y
   
      real(DP),parameter :: hcorrection = 1.25 !2.0_DP / SQRT2
   
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
      xi = iand(floor(x),255)
      yi = iand(floor(y),255)
      x = x - real(floor(x), kind=DP)
      y = y - real(floor(y), kind=DP)
      u = fade(x)
      v = fade(y)
   
      A = p(xi) + yi
      B = p(xi + 1) + yi
   
      AA = p(A) 
      BA = p(B) 
   
      AB = p(A + 1) 
      BB = p(B + 1) 
   
      GA1 = grad(p(AA), x        , y        )
      GB1 = grad(p(BA), x - 1._DP, y        )
      GA2 = grad(p(AB), x        , y - 1._DP)
      GB2 = grad(p(BB), x - 1._DP, y - 1._DP)
   
      L1 = lerp(u, GA1, GB1)
      L2 = lerp(u, GA2, GB2)
     
      noise = hcorrection * lerp(v, L1, L2)
   
      if (present(dx).and.present(dy)) then
   
         ! Calculate derivatives wrt x
         ux = dfade(x)
         vx = 0.0_DP
         GA1x = grad(p(AA),1._DP,0._DP)
         GB1x = grad(p(BA),1._DP,0._DP)
         GA2x = grad(p(AB),1._DP,0._DP)
         GB2x = grad(p(BB),1._DP,0._DP)
   
         L1x = dlerp(u, GA1, GB1, ux, GA1x, GB1x)
         L2x = dlerp(u, GA2, GB2, ux, GA2x, GB2x)
   
         dx = hcorrection * dlerp(v, L1, L2, vx, L1x, L2x)
   
         ! Calculate derivatives wrt y
         uy = 0.0_DP
         vy = dfade(y)
         GA1y = grad(p(AA),0._DP,1._DP)
         GB1y = grad(p(BA),0._DP,1._DP)
         GA2y = grad(p(AB),0._DP,1._DP)
         GB2y = grad(p(BB),0._DP,1._DP)
   
         L1y = dlerp(u, GA1, GB1, uy, GA1y, GB1y)
         L2y = dlerp(u, GA2, GB2, uy, GA2y, GB2y)
   
         dy = hcorrection * dlerp(v, L1, L2, vy, L1y, L2y)
   
      end if
   
      return

      contains
      
         function fade(t) result(f)
            implicit none
            real(DP),intent(in) :: t
            real(DP) :: f
      
            f = t * t * t * (t * (t * 6._DP - 15._DP) + 10._DP)
            return
            end function fade
      
            function dfade(t) result(df)
            implicit none
            real(DP),intent(in) :: t
            real(DP) :: df
      
            df = t * t * (t * (t * 30._DP - 60._DP) + 30._DP)
            return
            end function dfade
      
            function lerp(t,a,b) result(l)
            implicit none
            real(DP),intent(in) :: t,a,b
            real(DP) :: l
      
            l = a + t * (b - a)
            return
         end function lerp
   
         function dlerp(t,a,b,dt,da,db) result(dl)
            implicit none
            real(DP),intent(in) :: t,a,b
            real(DP),intent(in) :: dt,da,db
            real(DP) :: dl
      
            dl = da + t * (db - da) + dt * (b - a)

            return
         end function dlerp
   
   
         function grad(hash,x,y) result(g)
            implicit none
            integer(I4B),intent(in) :: hash
            real(DP),intent(in) :: x,y
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
               v = 0._DP
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
   
   
   end subroutine util_perlin_noise


   module function util_perlin_turbulence(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, noise_height, freq, pers
      integer(I4B) :: num_octaves
      real(DP), dimension(3,num_octaves),intent(in) :: anchor
      real(DP) :: noise,xnew,ynew,norm,rx,ry,thetai
      
      integer(I4B) :: i
      real(DP) :: spatial_fac, noise_mag,dn,dx,dy
      
      noise = 0.0_DP
      norm = 0.5_DP
      do i = 1, num_octaves
         spatial_fac = freq ** (i - 1)
         noise_mag = pers ** (i - 1)
         norm = norm + 0.5_DP * noise_mag
         thetai = anchor(3,i)
         rx = x * cos(thetai) - y * sin(thetai)
         ry = x * sin(thetai) + y * cos(thetai)
         xnew = (rx + anchor(1,i)) * spatial_fac
         ynew = (ry + anchor(2,i)) * spatial_fac
         call util_perlin_noise(xnew,ynew,dn)
         noise = noise + dn * noise_mag
      end do
      
      noise = noise * noise_height / norm
      
      return
      
   end function util_perlin_turbulence


   module function util_perlin_arrayinput(x, y, num_octaves, Sarr, Aarr, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y
      integer(I4B) :: num_octaves
      real(DP), dimension(num_octaves),intent(in) :: Sarr, Aarr
      real(DP), dimension(3,num_octaves),intent(in) :: anchor
      real(DP) :: noise
      
      !Internal variables
      real(DP) :: xnew,ynew,norm,rx,ry,thetai
      integer(I4B) :: i
      real(DP) :: spatial_fac, noise_mag,dn,dx,dy
      
      noise = 0.0_DP
      norm = 0.5_DP
      do i = 1, num_octaves
         spatial_fac = Sarr(i)
         noise_mag = Aarr(i)
         norm = norm + 0.5_DP * noise_mag
         thetai = anchor(3,i)
         rx = x * cos(thetai) - y * sin(thetai)
         ry = x * sin(thetai) + y * cos(thetai)
         xnew = (rx + anchor(1,i)) * spatial_fac
         ynew = (ry + anchor(2,i)) * spatial_fac
         call util_perlin_noise(xnew,ynew,dn)
         noise = noise + dn * noise_mag / norm
      end do
      
      noise = noise !/ norm
      
      return
   end function util_perlin_arrayinput   


   module function util_perlin_billowedNoise(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, noise_height, freq, pers
      integer(I4B) :: num_octaves
      real(DP),dimension(3,num_octaves),intent(in) :: anchor
      real(DP) :: noise
      
      noise = abs(util_perlin_turbulence(x, y, noise_height, freq, pers, num_octaves, anchor))
      
      
      return
      
   end function util_perlin_billowedNoise
      
      
   module function util_perlin_plawNoise(x, y, noise_height, freq, pers, slope, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, noise_height, freq, pers,slope
      integer(I4B) :: num_octaves
      real(DP),dimension(3,num_octaves),intent(in) :: anchor
      real(DP) :: noise
      
      noise = util_perlin_turbulence(x, y, noise_height, freq, pers, num_octaves, anchor)
      noise = sign(abs(noise)**(slope),noise)
      
      return
   end function util_perlin_plawNoise
      
      
   module function util_perlin_ridgedNoise(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, noise_height, freq, pers
      integer(I4B) :: num_octaves
      real(DP),dimension(3,num_octaves),intent(in) :: anchor
      real(DP) :: noise
      
      noise = noise_height - abs(util_perlin_turbulence(x, y, noise_height, freq, pers, num_octaves, anchor))
      
      return
   end function util_perlin_ridgedNoise
      
      
   module function util_perlin_swissTurbulence(x, y, lacunarity, gain, warp, num_octaves, anchor) result(noise)
      implicit none
      real(DP), intent(in) ::  x, y, lacunarity, gain, warp
      integer(I4B) :: num_octaves
      real(DP),dimension(3,num_octaves) :: anchor
      real(DP) :: noise
      real(DP) :: freq, amp,newx,newy,norm
      real(DP),dimension(2) :: dsum
      real(DP),dimension(3) :: n
      integer(I4B) :: i
   
      noise =  0.0_DP
      freq = 1.0_DP
      amp = 1.0_DP
      dsum(:) = 0._DP
   
      norm = 0.5_DP
         
      do i = 1, num_octaves 
         newx = (x + anchor(1,i) + warp * dsum(1)) * freq
         newy = (y + anchor(2,i) + warp * dsum(2)) * freq
         call util_perlin_noise(newx,newy,n(1),dx = n(2),dy = n(3))
         noise = noise + amp * (1._DP - abs(n(1)))
         norm = norm + amp
         dsum(:) = dsum(:) + amp * n(2:3) * (-n(1))
         freq = freq * lacunarity;
         amp = amp * gain * max(min(noise,1.0_DP),0.0_DP)
      end do
      noise = noise  !/ norm
      return
   
   end function util_perlin_swissTurbulence
   
      
   module function util_perlin_jordanTurbulence(x, y, lacunarity, gain1, gain, warp0, warp, damp0, damp, damp_scale,&
            num_octaves, anchor) result(noise)
      ! Fortran implementation of noise function by Giliam de Carpentier
      implicit none
      real(DP),intent(in) :: x, y, lacunarity, gain1, gain, warp0, warp, damp0, damp, damp_scale
      integer(I4B),intent(in) :: num_octaves
      real(DP), dimension(3,num_octaves), intent(in) :: anchor
      ! Result variable
      real(DP) :: noise
   
      ! Internal variables
      real(DP),dimension(3) :: n,n2
      real(DP),dimension(2) :: dsum_warp,dsum_damp
      real(DP) :: amp,freq,damped_amp,xnew,ynew
      integer(I4B) :: i
   
      xnew = x + anchor(1,1)
      ynew = y + anchor(2,1)
      call util_perlin_noise(xnew,ynew,n(1),dx = n(2),dy = n(3))
      n2(:) = n(:) * n(1)
      noise = n2(1)
      dsum_warp(:) = warp0 * n2(2:3)
      dsum_damp(:) = damp0 * n2(2:3)
   
      amp = gain1
      freq = lacunarity
      damped_amp = amp * gain
   
      do i = 2, num_octaves
         xnew = x * freq + dsum_warp(1) + anchor(1,i)
         ynew = y * freq + dsum_warp(2) + anchor(2,i)
         call util_perlin_noise(xnew,ynew,n(1),dx = n(2),dy = n(3))
         n2(:) = n(:) * n(1)
         noise = noise + damped_amp * n2(1)
         dsum_warp(:) = dsum_warp(:) + warp * n2(2:3)
         dsum_damp(:) = dsum_damp(:) + damp * n2(2:3)
         freq = freq * lacunarity
         amp = amp * gain
         damped_amp = amp * (1._DP - damp_scale / (1._DP + dot_product(dsum_damp(:),dsum_damp(:))))
      end do
      return 
   end function util_perlin_jordanTurbulence
      
      

end submodule s_util_perlin
