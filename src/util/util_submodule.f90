!! Copyright 2023 - David Minton
!! This file is part of Cratermaker
!! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Cratermaker. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (util) s_util
   use globals
contains

   module function util_perlin_noise(xx,yy,zz) result(noise)
      !! author: David A. Minton
      !!
      !! Modern Fortran implementation of Ken Perlin's 2002 noise generator from https://mrl.cs.nyu.edu/~perlin/noise/
      implicit none

      ! Arguments
      real(DP),intent(in) :: xx,yy
      real(DP),intent(in),optional :: zz
      ! Result
      real(DP) :: noise
      ! Internal variables
      integer(I4B) :: xi,yi,zi,A,B,AA,BB,AB,BA
      real(DP) :: x,y,z,u,v,w
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
      if (present(zz)) then
         z = zz
      else
         z = 0.0_DP 
      end if
      xi = iand(floor(x),255)
      yi = iand(floor(y),255)
      zi = iand(floor(z),255)
      x = x - floor(x)
      y = y - floor(y)
      z = z - floor(z)
      u = fade(x)
      v = fade(y)
      w = fade(z)

      A = p(xi) + yi
      B = p(xi + 1) + yi

      AA = p(A) + zi
      BA = p(B) + zi

      AB = p(A + 1) + zi
      BB = p(B + 1) + zi
      
      noise = lerp(w,   lerp(v,  lerp(u,  grad(p(AA    ), x    , y,     z    ),&
                                          grad(p(BA    ), x - 1, y,     z    )),& 
                                 lerp(u,  grad(p(AB    ), x    , y - 1, z    ),&  
                                          grad(p(BB    ), x - 1, y - 1, z    ))),&
                        lerp(v,  lerp(u,  grad(p(AA + 1), x    , y    , z - 1),&
                                          grad(p(BA + 1), x - 1, y    , z - 1)),& 
                                 lerp(u,  grad(p(AB + 1), x    , y - 1, z - 1),&  
                                          grad(p(BB + 1), x - 1, y - 1, z - 1))))
      return
      
      contains
         function fade(t) result(f)
            implicit none
            real(DP),intent(in) :: t
            real(DP) :: f

            f = t * t * t * (t * (t * 6._DP - 15._DP) + 10_DP)
            return
            end function fade

            function lerp(t,a,b) result(l)
            implicit none
            real(DP),intent(in) :: t,a,b
            real(DP) :: l

            l = a + t * (b - a)
            return
         end function lerp

         function grad(hash,x,y,z) result(g)
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

   end function util_perlin_noise

end submodule s_util
