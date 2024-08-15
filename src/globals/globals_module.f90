! Copyright 2023 - David Minton
! This file is part of Cratermaker
! Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with cratermaker. 
! If not, see: https://www.gnu.org/licenses. 

!! Type definitions and global parameters for Cratermaker project. All types are defined using ISO_C_BINDING types to ensure that
!! the types are compatible with the Cython interface.
module globals
    !! author: David A. Minton
    !!
    !! Basic parameters, definitions, and global type definitions used throughout the Cratermaker project
    use, intrinsic :: iso_c_binding  ! Use the intrinsic kind definitions
    implicit none
    public

    integer, parameter :: I8B = c_int_least64_t !! Symbolic name for kind types of 8-byte integers
    integer, parameter :: I4B = c_int_least32_t !! Symbolic name for kind types of 4-byte integers
    integer, parameter :: I2B = c_int_least16_t !! Symbolic name for kind types of 2-byte integers
    integer, parameter :: I1B = c_int_least8_t  !! Symbolic name for kind types of 1-byte integers

    integer, parameter :: SP = c_float  !! Symbolic name for kind types of single-precision reals
    integer, parameter :: DP = c_double  !! Symbolic name for kind types of double-precision reals
    integer, parameter :: QP = c_long_double !! Symbolic name for kind types of quad-precision reals
    integer, parameter :: LGT = c_bool !! Symbolic name for kind types of logicals

    integer(I4B), parameter :: STRMAX = 512 !! Maximum size of character strings 

    ! Frequently used mathematical constants 
    real(DP), parameter :: PI       = ACOS(-1.0_DP) 
    real(DP), parameter :: SQRT2    = SQRT(2.0_DP) 
    real(DP), parameter :: LOGSQRT2 = LOG(SQRT2) 
    real(DP), parameter :: SQRT3    = SQRT(3.0_DP) 
    real(DP), parameter :: THIRD    = 1.0_DP / 3.0_DP 
    real(DP), parameter :: SIXTH    = 1.0_DP / 6.0_DP 
    real(DP), parameter :: DEG2RAD  = PI / 180.0_DP 

    real(DP),parameter :: VSMALL  = 10*tiny(1._DP)    ! Very small number
    real(DP),parameter :: LOGVSMALL = log(VSMALL)  ! log of a very small number
    real(DP),parameter :: VBIG    = huge(1._DP)    ! Very big number
    real(DP),parameter :: SMALLFAC = 1e-5_DP       ! Smallest unit of measurement proportional to pixel size   
end module globals
