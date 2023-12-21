# Copyright 2023 - David Minton
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

# CMake project file for Cratermaker
##################################################
# Define the project and the depencies that it has
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.6.0...3.27.1)
SET(SKBUILD_PROJECT_NAME "cratermaker" CACHE STRING "Name of project set by scikit-build")
SET(SKBUILD_SCRIPTS_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE STRING "Install location of binary executable")

# Get version stored in text file
FILE(READ "version.txt" VERSION_STRING)

# Regular expression to extract major, minor, and patch numbers
STRING(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ ${VERSION_STRING})

# Set the project version
SET(VERSION "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")


# Check if DOC_BUILD environment variable is set and skip build if it is
IF (DEFINED ENV{DOC_BUILD})
    MESSAGE(STATUS "DOC_BUILD is set. Skipping build process.")
    PROJECT(${SKBUILD_PROJECT_NAME} VERSION ${VERSION})
ELSE()
    PROJECT(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran VERSION ${VERSION})

    # Use the old method to get Python packages, as that's what scikit-build uses
    IF (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.27")
        CMAKE_POLICY(SET CMP0148 OLD)
    ENDIF ()

    # Set some options the user may choose
    OPTION(USE_COARRAY "Use Coarray Fortran for parallelization of test particles" OFF)
    OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
    OPTION(USE_SIMD "Use SIMD vectorization" ON)
    OPTION(BUILD_SHARED_LIBS "Build using shared libraries" ON)

    INCLUDE(GNUInstallDirs)
    IF (SKBUILD)
        SET(INSTALL_PREFIX ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME})
        SET(INSTALL_BINDIR ${INSTALL_PREFIX}/bin) 
        SET(INSTALL_LIBDIR ${INSTALL_PREFIX}/lib) 
        SET(INSTALL_INCLUDEDIR ${INSTALL_LIBDIR}) 
        SET(INSTALL_CYTHON ${INSTALL_PREFIX}/cython)
        IF (APPLE)
            SET(CMAKE_INSTALL_RPATH "@loader_path/../lib")
        ELSEIF (LINUX)
            SET(CMAKE_INSTALL_RPATH $ORIGIN/../lib)
        ENDIF ()
    ELSE ()
        SET(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
        SET(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
        SET(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
        SET(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})
        SET(INSTALL_CYTHON ${INSTALL_LIBDIR})
    ENDIF ()

    # Have the .mod files placed in the include folder
    SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

    # Add our local modules to the module ldpath
    FILE(TO_CMAKE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules" LOCAL_MODULE_PATH)
    LIST(APPEND CMAKE_MODULE_PATH ${LOCAL_MODULE_PATH})

    IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
        SET(COMPILER_OPTIONS "Intel" CACHE STRING "Compiler identified as Intel")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        SET(COMPILER_OPTIONS "GNU" CACHE STRING "Compiler identified as gfortran")
    ELSE ()
        MESSAGE(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER_ID} not recognized!") 
    ENDIF ()

    # The following section is modified from Numpy f2py documentation
    IF(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    MESSAGE(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
    ENDIF()

    # Ensure scikit-build modules
    FIND_PACKAGE(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)


    SET(CYTHON_SRC "${CMAKE_SOURCE_DIR}/cratermaker/cython")
    SET(SRC "${CMAKE_SOURCE_DIR}/src")

    # Make sure paths are correct for Unix or Windows style
    FILE(TO_CMAKE_PATH ${SRC} SRC)
    FILE(TO_CMAKE_PATH ${CYTHON_SRC} CYTHON_SRC)

    # Set the name of the cratermaker library
    SET(CRATERMAKER_LIBRARY "${SKBUILD_PROJECT_NAME}")

    # Determine compiler options
    IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
        MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
    ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    INCLUDE(SetParallelizationLibrary)

    IF (COMPILER_OPTIONS STREQUAL "Intel" AND NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
        FIND_PACKAGE(MKL REQUIRED)
    ENDIF ()
    INCLUDE(SetCratermakerFlags) 

    # First build the external dependencies
    INCLUDE(externalproject)

    ExternalProject_Add(jigsaw
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/jigsaw
        BINARY_DIR ${CMAKE_BINARY_DIR}/jigsaw-build
        CMAKE_ARGS 
            -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}
            -DCMAKE_INSTALL_BINDIR=${INSTALL_BINDIR}
            -DCMAKE_INSTALL_LIBDIR=${INSTALL_LIBDIR}
            -DCMAKE_BUILD_TYPE=Release
        BUILD_ALWAYS 1
    )

    # Compile the project source code
    ADD_SUBDIRECTORY(${SRC}) 
    ADD_SUBDIRECTORY(${CYTHON_SRC})

ENDIF(DEFINED ENV{DOC_BUILD})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_SOURCE_DIR}/distclean.cmake"
)
