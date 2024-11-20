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
CMAKE_MINIMUM_REQUIRED(VERSION 3.23.1...3.30.2)
IF (SKBUILD)
    SET(VERSION ${SKBUILD_PROJECT_VERSION})
    SET(PROJECT_NAME ${SKBUILD_PROJECT_NAME})
ELSE ()
    SET(VERSION "0.0.0")
    SET(PROJECT_NAME "cratermaker" CACHE STRING "Name of project")
ENDIF()

SET(SKBUILD_SCRIPTS_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE STRING "Install location of binary executable")
IF (SKBUILD)
    set(ENV_PREFIX "$ENV{PREFIX}")
    string(REPLACE "<ENV_PREFIX>" "${ENV_PREFIX}" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
ENDIF ()
IF(DEFINED ENV{CIBUILDWHEEL})
    SET(CIBUILDWHEEL ON)
ENDIF()

# Check if DOC_BUILD environment variable is set and skip build if it is
IF (DEFINED ENV{DOC_BUILD})
    MESSAGE(STATUS "DOC_BUILD is set. Skipping build process.")
    PROJECT(${PROJECT_NAME} VERSION ${VERSION})
ELSE()
    PROJECT(${PROJECT_NAME} LANGUAGES C CXX Fortran VERSION ${VERSION})

    # Use the old method to get Python packages, as that's what scikit-build uses
    IF (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.27")
        CMAKE_POLICY(SET CMP0148 OLD)
    ENDIF ()

    # Set some options the user may choose
    OPTION(USE_COARRAY "Use Coarray Fortran for parallelization of test particles" OFF)
    OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
    OPTION(USE_SIMD "Use SIMD vectorization" ON)
    OPTION(BUILD_SHARED_LIBS "Build using shared libraries" ON)

    SET(SRC "${CMAKE_SOURCE_DIR}/src")
    SET(PY "${CMAKE_SOURCE_DIR}/cratermaker")
    SET(FORTRAN_BINDINGS_SRC "${CMAKE_SOURCE_DIR}/cratermaker/fortran_bindings")
    FILE(TO_CMAKE_PATH ${SRC} SRC)
    FILE(TO_CMAKE_PATH ${PY} PY)
    FILE(TO_CMAKE_PATH ${FORTRAN_BINDINGS_SRC} FORTRAN_BINDINGS_SRC)

    INCLUDE(GNUInstallDirs)
    IF (SKBUILD)
        SET(CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/deps/usr/lib")
        SET(INSTALL_INCLUDEDIR ${SKBUILD_HEADERS_DIR}) 
        SET(INSTALL_PYPROJ ${SKBUILD_PLATLIB_DIR}/${PROJECT_NAME}) 
        SET(INSTALL_BINDIR ${INSTALL_PYPROJ})
        # The following is needed so that the shared libraries get installed into the correct place. cibuildwheel behaves differently 
        # than just running pip install directly. When running pip install directly, the libraries will be installed to the environment's
        # lib directory by specifing ${SKBUILD_DATA_DIR}/lib as the library directory. 
        # However, when running cibuildwheel, the libraries are repaired and placed in a special directory of all libraries
        # that are packed up into the wheel. Currently, cibuildwheel cannot find the compiled libraries if they are installed to there,
        # so we need to make sure the rpath is set to the original build directory so that the library repair step can reset the rpaths
        # when running cibuildwheel.
        IF (CIBUILDWHEEL)
            SET(INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)
            SET(CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/bin")
        ELSE()
            SET(INSTALL_LIBDIR ${INSTALL_BINDIR})
            IF (APPLE)
                SET(CMAKE_INSTALL_RPATH "@loader_path;${INSTALL_LIBDIR}")
            ELSEIF (LINUX)
                SET(CMAKE_INSTALL_RPATH "$ORIGIN;${INSTALL_LIBDIR}")
            ENDIF ()
        ENDIF()
    ELSE ()
        SET(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})
        SET(INSTALL_PYPROJ ${PY})
        SET(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
        SET(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
    ENDIF ()
    SET(INSTALL_FORTRAN_BINDINGS ${INSTALL_PYPROJ})

    MESSAGE(STATUS "INSTALL_BINDIR: ${INSTALL_BINDIR}")
    MESSAGE(STATUS "INSTALL_LIBDIR: ${INSTALL_LIBDIR}")
    MESSAGE(STATUS "INSTALL_INCLUDEDIR: ${INSTALL_INCLUDEDIR}")
    MESSAGE(STATUS "INSTALL_PYPROJ: ${INSTALL_PYPROJ}")
    MESSAGE(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # Have the .mod files placed in the include folder
    SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

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

    # Make sure paths are correct for Unix or Windows style
    FILE(TO_CMAKE_PATH ${SRC} SRC)
    FILE(TO_CMAKE_PATH ${FORTRAN_BINDINGS_SRC} FORTRAN_BINDINGS_SRC)

    # Set the name of the cratermaker library
    SET(CRATERMAKER_LIBRARY "${PROJECT_NAME}")

    # Determine compiler options
    IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
        MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
    ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    INCLUDE(SetParallelizationLibrary)

    INCLUDE(SetCratermakerFlags) 

    # Compile the project source code
    ADD_SUBDIRECTORY(${SRC}) 
    ADD_SUBDIRECTORY(${FORTRAN_BINDINGS_SRC})

ENDIF(DEFINED ENV{DOC_BUILD})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_SOURCE_DIR}/distclean.cmake"
)
