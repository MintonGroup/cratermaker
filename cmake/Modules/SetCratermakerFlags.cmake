# Copyright 2025 - David Minton
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 


######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(SetCompileFlag)

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

SET(BUILD_TYPE_MSG "Choose the type of build, options are DEBUG, RELEASE, PROFILE, or TESTING.")

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(BT STREQUAL "PROFILE")
    SET (CMAKE_BUILD_TYPE PROFILE CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      ${BUILD_TYPE_MSG}
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid! ${BUILD_TYPE_MSG}")
ENDIF(BT STREQUAL "RELEASE")

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG AND CMAKE_Fortran_FLAGS_PROFILE )
    RETURN ()
ENDIF()

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    SET(WINOPT True)
ELSE ()
    SET(WINOPT False)
ENDIF ()
#####################
### GENERAL FLAGS ###
#####################

# Free form
IF (COMPILER_OPTIONS STREQUAL "GNU")
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
        Fortran "-ffree-form" # GNU
        ) 
    # Don't add underscores in symbols for C-compatability
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
        Fortran "-fno-underscoring" # GNU
        ) 
    # Compile code assuming that IEEE signaling NaNs may generate user-visible traps during floating-point operations. 
    # Setting this option disables optimizations that may change the number of exceptions visible with signaling NaNs. 
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
        Fortran "-fsignaling-nans " # GNU
        ) 
    # Allows for lines longer than 80 characters without truncation
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
        Fortran "-ffree-line-length-512" # GNU (gfortran)
        )
    # Sets the dialect standard
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
        Fortran "-std=f2018" 
        )
ELSEIF (COMPILER_OPTIONS STREQUAL "Intel")
    # Disables right margin wrapping in list-directed output
    IF (WINOPT)
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran  "/wrap-margin-"   # Intel Windows    
        )
        # Aligns a variable to a specified boundary and offset
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "/align:all /align:array64byte" # Intel
        )
        # Enables changing the variable and array memory layout
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "/Qpad" # Intel Windows
        )
        #Enables byte record length
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "/assume:byterecl" # Intel Windows
        )
        #Enables heap arrays
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "/heap-arrays" 
        )
    ELSE ()
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran  "-no-wrap-margin" # Intel
        )
        # Aligns a variable to a specified boundary and offset
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "-align all -align array64byte" # Intel
        )
        # Enables changing the variable and array memory layout
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "-pad" # Intel Windows
        )
        #Enables byte record length
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "-assume byterecl" 
        )
        #Enables heap arrays
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
            Fortran "-heap-arrays" 
        )
    ENDIF ()
ENDIF ()

IF (NOT BUILD_SHARED_LIBS AND NOT WINOPT)
    SET_COMPILE_FLAG(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS}"
        Fortran "-fPIC"
        )
    IF (COMPILER_OPTIONS STREQUAL "Intel")
        # Use static Intel libraries
        SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
            Fortran "-static-intel"  # Intel
        )
        # Use static Intel MPI libraries
        SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
            Fortran "-static_mpi"  # Intel
        )
        IF (USE_OPENMP)
            SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran "-qopenmp-link=static"  # Intel
            )
        ENDIF (USE_OPENMP)
    ELSEIF (COMPILER_OPTIONS STREQUAL "GNU") 
        IF (NOT BUILD_SHARED_LIBS) 
            # Set GNU static libraries
            SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran  "-static-libgfortran" 
            )
            SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran  "-static-libgcc" 
            )
            SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran  "-static-libstdc++" 
            )
            SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran  "-static-libquadmath" 
            )
        ENDIF ()
        IF (USE_OPENMP)
            SET_COMPILE_FLAG(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS}"
                Fortran "-lomp"  
                        "-lgomp"  
            )
        ENDIF (USE_OPENMP)
    ENDIF ()
ENDIF ()
IF (USE_SIMD)
    IF (COMPILER_OPTIONS STREQUAL "Intel")
        SET(MACHINE_CODE_VALUE "Host" CACHE STRING "Tells the compiler which processor features it may target, including which instruction sets and optimizations it may generate.")

        IF (MACHINE_CODE_VALUE STREQUAL "generic")
            SET(MACHINE_CODE_VALUE "SSE2" CACHE STRING "SSE2 is the safest option when compiling for non-host compatibility" FORCE)
        ENDIF()

        # Enables OpenMP SIMD compilation when OpenMP parallelization is disabled. 
        IF (NOT USE_OPENMP)
            IF (WINOPT) 
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "/Qopenmp- /Qopenmp-simd" # Intel
            )
            ELSE ()
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-qno-openmp -qopenmp-simd" # Intel
            )
            ENDIF ()     
        ENDIF (NOT USE_OPENMP)

        IF (WINOPT)
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "/Qx${MACHINE_CODE_VALUE}" # Intel
            )
            # Generate an extended set of vector functions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "/Qvecabi:cmdtarget" # Intel Windows
            )
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "/Qsimd-honor-fp-model" # Intel
            )
        ELSE ()
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-x${MACHINE_CODE_VALUE}" # Intel
            )
            # Generate an extended set of vector functions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-vecabi=cmdtarget" # Intel
            )
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-qsimd-honor-fp-model" # Intel
            )
        ENDIF ()

    ELSEIF (COMPILER_OPTIONS STREQUAL "GNU")
        # Enables OpenMP SIMD compilation when OpenMP parallelization is disabled. 
        IF (NOT USE_OPENMP)
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-fno-openmp -fopenmp-simd" # GNU
            )     
        ENDIF (NOT USE_OPENMP)

        SET(MACHINE_CODE_VALUE "native" CACHE STRING "Tells the compiler which processor features it may target, including which instruction sets and optimizations it may generate.")
        IF (MACHINE_CODE_VALUE STREQUAL "Host" OR MACHINE_CODE_VALUE STREQUAL "native")
            SET(MARCH_VALUE "native" CACHE STRING "value passed to -march=")
            SET(MTUNE_VALUE "native" CACHE STRING "value passed to -mtune=")
        ELSEIF (MACHINE_CODE_VALUE STREQUAL "generic")
            IF (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
                SET(MARCH_VALUE "x86-64" CACHE STRING "value passed to -march=")
                SET(MTUNE_VALUE "generic" CACHE STRING "value passed to -mtune=")
            ELSEIF (CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "aarch64")
                SET(MARCH_VALUE "armv8-a" CACHE STRING "value passed to -march=")
                SET(MTUNE_VALUE "generic" CACHE STRING "value passed to -mtune=")
            ENDIF ()
        ELSE () 
            SET(MARCH_VALUE "${MACHINE_CODE_VALUE}" CACHE STRING "value passed to -march=" )
            SET(MTUNE_VALUE "${MACHINE_CODE_VALUE}" CACHE STRING "value passed to -mtune=")
        ENDIF ()
        
        IF (APPLE)
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-mtune=${MTUNE_VALUE}" 
            )
        ELSE ()
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-march=${MARCH_VALUE}" 
            )
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                Fortran "-mtune=${MTUNE_VALUE}" 
            )
        ENDIF ()
    ENDIF ()
ENDIF (USE_SIMD)    

###################
### DEBUG FLAGS ###
###################
IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "TESTING" )
    # Disable optimizations
    IF (COMPILER_OPTIONS STREQUAL "Intel")
        IF (WINOPT)
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran REQUIRED "/Od" # Intel Windows
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C REQUIRED "/0d"
            )
            # Turn on all warnings 
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/warn:all" # Intel Windows
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "/warn:all" # Intel Windows
            )
            # Tells the compiler to issue compile-time messages for nonstandard language elements (Fortran 2018).    
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/stand:f18"  # Intel Windows
            )  
            # Traceback
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/traceback"   # Intel Windows
            )
            # Check everything
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/check:all"      # Intel Windows
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "/check:all"      # Intel Windows
            )
            # Initializes matrices/arrays with NaN values
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/Qinit:snan,arrays" # Intel Windows
            )
            # Does not generate an interface block for each routine in a source file
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/nogen-interfaces" # Intel Windows
            )

            # Does not set denormal results from floating-point calculations to zero
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/Qftz-"  # Intel Windows
            )
            # Enables floating-point invalid, divide-by-zero, and overflow exceptions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/fpe-all:0"         # Intel Windows
            )
            # Enables floating-point invalid, divide-by-zero, and overflow exceptions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran  "/fpe:0" # Intel Windows
            )
            # Enables debug info
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/debug:all" # Intel Windows
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "/debug:all" # Intel Windows
            )
            # Disables additional interprocedural optimizations for a single file compilation
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/Qip-"  # Intel Windows
            )
            # Disables prefetch insertion optimization
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "/Qopt-prefetch-"   # Intel Windows
            )
        ELSE ()
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran REQUIRED "-O0" # All compilers not on Windows
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C REQUIRED "-O0" # All compilers not on Windows
            )
            # Turn on all warnings 
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-warn all" # Intel
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "-Wall" # Intel
            )
            # Tells the compiler to issue compile-time messages for nonstandard language elements (Fortran 2018).    
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-stand f18"  # Intel
            )  
            # Traceback
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-traceback"   # Intel Group
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "-traceback"   # Intel Group
            )
            # Check everything
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-check all"       # Intel
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "-check=conversions,stack,uninit"       # Intel
            )
            # Initializes matrices/arrays with NaN values
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-init=snan,arrays"  # Intel
            )
            # Does not generate an interface block for each routine in a source file
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-nogen-interfaces" # Intel
            )
            # Does not generate aposition independent executable
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-no-pie" # Intel
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "-no-pie" # Intel
            )
            # Does not set denormal results from floating-point calculations to zero
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-no-ftz" # Intel
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "-no-ftz" # Intel
            )
            # Enables floating-point invalid, divide-by-zero, and overflow exceptions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-fpe-all=0"         # Intel
            )
            # Enables floating-point invalid, divide-by-zero, and overflow exceptions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-fpe0"  # Intel 
            )
            # Enables debug info
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-debug all" # Intel
            )
            SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
                C "-debug all" # Intel
            )
            # Disables additional interprocedural optimizations for a single file compilation
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-no-ip" # Intel
            )
            # Disables prefetch insertion optimization
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                Fortran "-qno-opt-prefetch" # Intel
            )
        ENDIF ()
    ELSEIF (COMPILER_OPTIONS STREQUAL "GNU")
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran REQUIRED "-Og" # GNU (gfortran)
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C REQUIRED "-Og" # GNU (gfortran)
        )
        # Turn on all warnings 
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-Wall"     # GNU
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-Wall"     # GNU
        )
        # This enables some extra warning flags that are not enabled by -Wall
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-Wextra" # GNU
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-Wextra" # GNU
        )
        # Disable the warning that arrays may be uninitialized, which comes up due to a known bug in gfortran
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-Wno-maybe-uninitialized" # GNU
        )
        # Disable the warning about unused dummy arguments. These primarily occur due to interface rules for type-bound procedures used in extendable types.
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-Wno-unused-dummy-argument" # GNU
        )
        # Traceback
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-fbacktrace"  # GNU (gfortran)
        )
        # Check everything
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-fcheck=all" # GNU 
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-fcheck=all" # GNU 
        )
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-fstack-check" # GNU 
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-fstack-check" # GNU 
        )
        # Initializes matrices/arrays with NaN values
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-finit-real=snan"   # GNU
        )
        # Generates non position-independent code
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-fno-PIE" # GNU
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-fno-PIE" # GNU
        )
        # Enables floating-point invalid, divide-by-zero, and overflow exceptions
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-ffpe-trap=zero,overflow,underflow" # GNU
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-ffpe-trap=zero,overflow,underflow" # GNU
        )
        # List of floating-point exceptions, whose flag status is printed to ERROR_UNIT when invoking STOP and ERROR STOP
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-ffpe-summary=all" # GNU
        )
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran "-fstack-check" # GNU 
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
            C "-fstack-check" # GNU 
        )
    ENDIF ()
ENDIF ()

#####################
### TESTING FLAGS ###
#####################

IF (CMAKE_BUILD_TYPE STREQUAL "TESTING" )
    # Optimizations
    IF (WINOPT)
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran REQUIRED "/O3" # Intel Windows
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_TESTING "${CMAKE_C_FLAGS_DEBUG}"
            C REQUIRED "/O3" # Intel Windows
        )
    ELSE ()
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_DEBUG}"
            Fortran REQUIRED "-O3" # All compilers not on Windows
        )
        SET_COMPILE_FLAG(CMAKE_C_FLAGS_TESTING "${CMAKE_C_FLAGS_DEBUG}"
            C REQUIRED "-O3" # All compilers not on Windows
        )
    ENDIF ()
ENDIF ()

#####################
### RELEASE FLAGS ###
#####################
# NOTE: agressive optimizations (-O3) are already turned on by default

IF (CMAKE_BUILD_TYPE STREQUAL "RELEASE" OR CMAKE_BUILD_TYPE STREQUAL "PROFILE")
    IF (COMPILER_OPTIONS STREQUAL "Intel")
        IF (WINOPT)
            # Unroll loops
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/Qunroll"    # Intel Windows
            )
            # Inline functions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/inline"       # Intel Windows
            )
            # Calls the Matrix Multiply library
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/Qopt-matmul" # Intel Windows
            )
            # Aligns a variable to a specified boundary and offset
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/align:all" # Intel Windows
            )
            # No floating-point exceptions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/fp:except-"       # Intel Windows
            )
            # Generate fused multiply-add instructions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/Qfma" # Intel Windows
            )
            # Tells the compiler to link to certain libraries in the Intel oneAPI Math Kernel Library (oneMKL). 
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/Qmkl:cluster" # Intel Windows
                        "/Qmkl"     # Intel Windows
            ) 
            # Enables additional interprocedural optimizations for a single file compilation
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/Qip" # Intel Windows
            )
        ELSE ()
            # Unroll loops
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-unroll"    # Intel
            )
            # Inline functions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-inline"    # Intel
            )

            # Calls the Matrix Multiply library
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-qopt-matmul" # Intel
            )
            # Aligns a variable to a specified boundary and offset
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-align all" # Intel
            )
            # No floating-point exceptions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-fp-model no-except" # Intel
            )

            # Generate fused multiply-add instructions
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-fma"  # Intel
            )
    ENDIF ()

    ELSEIF(COMPILER_OPTIONS STREQUAL "GNU")
        # Unroll loops
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
            Fortran "-funroll-loops" # GNU
        )
        # Inline functions
        SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
            Fortran "-finline-functions" # GNU
        )
    ENDIF ()
ENDIF ()
 
#####################
### MATH FLAGS ###
#####################
IF (COMPILER_OPTIONS STREQUAL "Intel")
    IF (WINOPT)
        # Some subroutines require more strict floating point operation optimizations for repeatability
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "/fp:precise" # Intel Windows 
        )
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "/Qprec-div" # Intel Windows 
        ) 
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "/Qprec-sqrt" # Intel Windows 
        )
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "/assume:protect-parens" # Intel Windows 
        ) 
        # Improves floating-point precision and consistency
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "/Qprec" # Intel Windows
        ) 
        # Most subroutines can use aggressive optimization of floating point operations without problems.       
        SET_COMPILE_FLAG(FASTMATH_FLAGS "${FASTMATH_FLAGS}"
            Fortran "/fp:fast"       # Intel Windows
        )
    ELSE ()
        # Some subroutines require more strict floating point operation optimizations for repeatability
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "-fp-model=precise" # Intel 
        )
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "-prec-div" # Intel 
        ) 
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "-prec-sqrt" # Intel 
        )
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "-assume protect-parens" # Intel
        ) 
        # Improves floating-point precision and consistency
        SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
            Fortran "-mp1" # Intel Windows
        ) 
        # Most subroutines can use aggressive optimization of floating point operations without problems.       
        SET_COMPILE_FLAG(FASTMATH_FLAGS "${FASTMATH_FLAGS}"
            Fortran "-fp-model=fast"       # Intel Windows
        )
    ENDIF ()    
ELSEIF (COMPILER_OPTIONS STREQUAL "GNU")
    # Some subroutines require more strict floating point operation optimizations for repeatability
    SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
        Fortran "-fno-unsafe-math-optimizations" # GNU
    )
    # Disable transformations and optimizations that assume default floating-point rounding behavior. 
    SET_COMPILE_FLAG(STRICTMATH_FLAGS "${STRICTMATH_FLAGS}"
        Fortran "-frounding-math"
    )
    # Most subroutines can use aggressive optimization of floating point operations without problems.       
    SET_COMPILE_FLAG(FASTMATH_FLAGS "${FASTMATH_FLAGS}"
        Fortran "-ffast-math"    # GNU
    )
ENDIF ()

# Debug mode always uses strict math
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" 
    Fortran ${STRICTMATH_FLAGS}
)

#####################
### PROFILE FLAGS ###
#####################
IF (CMAKE_BUILD_TYPE STREQUAL "PROFILE")
    IF (COMPILER_OPTIONS STREQUAL "Intel")
    # Enables the optimization reports to be generated
        IF (WINOPT)
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "/O2 /Qopt-report:5 /traceback /Z7"    # Intel Windows
            )
        ELSE ()
            SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE}"
                Fortran "-O2 -pg -qopt-report=5 -traceback -p -g3" # Intel
            )
        ENDIF ()
    ELSEIF (COMPILER_OPTIONS STREQUAL "GNU")
    # Enables the optimization reports to be generated
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE}"
        Fortran "-O2 -pg -fbacktrace"          # GNU
        )
    ENDIF ()
ENDIF ()
