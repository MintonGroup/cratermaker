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




type configtype
   ! Required input variables
   integer(I4B)      :: gridsize  ! Resolution
   integer(I4B)      :: numlayers ! Number of perched layers
   real(DP)          :: pix       ! Pixel size (m)
   real(DP)          :: mu        ! Crater scaling exponential constant (ignored for basins)
   real(DP)          :: kv        ! Crater scaling linear constant
   integer(I4B)      :: seed      ! Random number generator seed (only used in non-IDL driven mode)
   real(DP)          :: targetDensity      ! Target surface density
   real(DP)          :: Ybar      ! Targetstrength (Pa)
   real(DP)          :: bodyGravityAccel  ! Gravitational accel at target
   real(DP)          :: targetBodyRadius 
   real(DP)          :: projectileDensity   
   character(STRMAX) :: sfdfile   ! Name of size distribution file
   character(STRMAX) :: velfile   ! Name of velocity distribution file
   real(DP) :: seisk,cohaccel     ! seismic keff, cohesion breaking acceleration
   
   ! Optional input variables
   logical           :: dorealistic ! Set to T to enable realistic crater morphology. Default is F.
   logical           :: docollapse ! Set T to use the slope collapse model (turning off speeds up the code for testing)
   logical           :: doangle    ! Set to F to only do vertical impacts, otherwise do range of angles (default is T)
   logical           :: doporosity ! Porosity on/off flg. Set to F to turn the model off. Default F. 
   logical           :: domixing   ! Set to F to turn off regolith mixing (useful for test craters when you don't want to simulate gardening). Default is T.
   logical           :: doquasimc  ! set to T for quasi-MC run. Default F.
   real(DP)          :: basinimp  ! Impactor size to switch to multiring basin
   real(DP)          :: maxcrat   ! fraction that maximum crater can be relative to grid
   real(DP)          :: deplimit  ! complex crater depth limit

   ! Seismic input variables 
   logical ::  doseismic   ! Set to T if you want to do the seismic shaking model
   real(DP) :: seisq       ! Seismic energy attenuation quality factor (Q)
   real(DP) :: neff        ! impact seismic energy efficiency factor
   real(DP) :: tvel        ! target P-wave (body wave) speed (m/s)
   real(DP) :: tfrac       ! mean free path for seismic wave scattering in medium
   real(DP) :: regcoh      ! target surface regolith layer cohesion

   ! Crater diffusion input parameters
   real(DP) :: Kd1 ! Degradation function coefficient (from Minton et al. (2019))
   real(DP) :: psi ! Degradation function exponent (from Minton et al. (2029))
   real(DP) :: psi2 ! Degradation function large size exponent (from Minton et al. (2020))
   real(DP) :: rbreak ! Degradation function break in exponent (from Minton et al. (2020))
   real(DP) :: fe  ! Scale factor for size of degradation region (from Minton et al. (2019))

   ! Ejecta softening variables
   logical           :: dosoftening  ! Set T to use the extra crater softening model
   real(DP)          :: ejecta_truncation ! Set the number of crater radii to truncate the ejecta
   logical           :: dorays       ! Set T to use ray model
   logical           :: superdomain  ! Set T to include the superdomain

   ! Regolith tracking variables
   logical           :: doregotrack ! Set T to use the regolith tracking model (EXPERIMENTAL)
   ! Scouring variables
   logical           :: doscour  ! Set T to use the ejecta scouring model (EXPERIMENTAL)
   ! Crustal thinning variables
   logical           :: docrustal_thinning ! Set T to use the crustal thinning model (EXPERIMENTAL)
   ! Diffusion variables
   logical           :: dotopodiffusion ! set T to do subpixel diffusion; default T; should be F for melt distribution runs

   logical           :: killatmaxcrater ! Set T to end the run when a crater exceeds the maximum allowable size
   
   ! Test input variables
   logical           :: testflag     ! Set to T if you want a single crater centered on the grid
   real(DP)          :: testimp      ! Test impactor size
   real(DP)          :: testvel      ! Test impactor velocity
   real(DP)          :: testang      ! Test impactor angle
   real(DP)          :: testxoffset  ! Offset of test crater from center in x direction (m)  
   real(DP)          :: testyoffset  ! Offset of test crater from center in y direction (m)   
   logical           :: tallyonly    ! Only run the tally routine (don't generate any new craters)
   logical           :: testtally    ! Set to T to count all non-cookie cut craters, regardless of score
   real(DP)          :: rctime       ! time (in interval units) for emplacement of quasi-MC crater

   ! IDL driver variables
   character(STRMAX) :: impfile      ! Name of impactor size distribution file (impacts per m^2 per y)
   real(DP)          :: interval     ! Length of interval between outputs (y)
   integer(I4B)      :: numintervals ! Total number of intervals
   character(STRMAX) :: runtype      ! Type of run: single or statistical
   logical           :: restart      ! Set to T restart an old run
   logical           :: popupconsole ! Pop up console window every output interval 
   logical           :: saveshaded   ! Output shaded relief images  
   logical           :: saverego     ! Output regolith map images 
   logical           :: savecomp     ! Output composition map images
   logical           :: savepres     ! Output simplified console display images (presentation-compatible images) 
   logical           :: savetruelist ! Save the true cumulative crater distribution for each interval (large file size)
   real(DP)          :: shadedminh   ! Minimum height for shaded relief map (m)
   real(DP)          :: shadedmaxh   ! Maximum height for shaded relief map (m)
   character(STRMAX) :: sfdcompare   ! Type of run: 0 for normal, 1 for statistical (domain is reset between intervals)
   character(STRMAX) :: realcraterlist ! This is only included here so the terminal doesn't return "Unknown parameter"
end type configtype

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