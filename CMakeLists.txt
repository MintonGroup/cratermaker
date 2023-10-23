# Copyright 2023 - David Minton
# This file is part of Cratermaker
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

# CMake project file for Cratermaker
##################################################
# Define the project and the depencies that it has
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.6.0...3.27.1)

# Get version stored in text file
FILE(READ "version.txt" VERSION)
PROJECT(CTEM LANGUAGES C Fortran VERSION ${VERSION})

IF (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    SET(COMPILER_OPTIONS "Intel" CACHE STRING "Compiler identified as Intel")
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(COMPILER_OPTIONS "GNU" CACHE STRING "Compiler identified as gfortran")
ELSE ()
    MESSAGE(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER_ID} not recognized!") 
ENDIF ()

# Set some options the user may choose
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
OPTION(USE_SIMD "Use SIMD vectorization" ON)
OPTION(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# The following section is modified from Numpy f2py documentation
IF(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  MESSAGE(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
ENDIF()

# Ensure scikit-build modules
FIND_PACKAGE(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Add our local modules to the module path
FILE(TO_CMAKE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules" LOCAL_MODULE_PATH)
LIST(APPEND CMAKE_MODULE_PATH ${LOCAL_MODULE_PATH})

# Define some directories that are important to the build
SET(SRC "${CMAKE_SOURCE_DIR}/src")
SET(LIB "${CMAKE_SOURCE_DIR}/lib")
SET(BIN "${CMAKE_SOURCE_DIR}/bin")
SET(MOD "${CMAKE_SOURCE_DIR}/include")
SET(PY "${CMAKE_SOURCE_DIR}/cratermaker")

# Make sure paths are correct for Unix or Windows style
FILE(TO_CMAKE_PATH ${SRC} SRC)
FILE(TO_CMAKE_PATH ${LIB} LIB)
FILE(TO_CMAKE_PATH ${BIN} BIN)
FILE(TO_CMAKE_PATH ${MOD} MOD)
FILE(TO_CMAKE_PATH ${PY} PY)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

# Set the name of the Cratermaker library
SET(CRATERMAKER_LIBRARY cratermaker)

# The source for the CTEM binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRC} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_SOURCE_DIR}/distclean.cmake"
)
