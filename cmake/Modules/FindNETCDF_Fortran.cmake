# Copyright 2023 - David Minton
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 


# - Finds the NetCDF libraries 
# Tries to find the cmake config files first. Otherwise, try to find the libraries and headers by hand

IF (NOT netCDF-Fortran_DIR)
    IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
       FILE(GLOB LIBDIRS "C:/Program Files*/NC4F")
       LIST(SORT LIBDIRS)
       LIST(GET LIBDIRS -1 LIBPREFIX)
       SET(netCDF-Fortran_DIR "${LIBPREFIX}/lib/cmake/netCDF" CACHE PATH "Location of provided netCDF-FortranConfig.cmake file")
    ELSE()
       SET(netCDF-Fortran_DIR "/usr/local/lib/cmake/netCDF" CACHE PATH "Location of provided netCDF-FortranConfig.cmake file")
    ENDIF ()
ENDIF()

MESSAGE(STATUS "Looking for netCDF-FortranConfig.cmake in ${netCDF-Fortran_DIR}")
FIND_PACKAGE(netCDF-Fortran QUIET)
IF (netCDF-Fortran_FOUND) 
   MESSAGE(STATUS "Found the NetCDF library via cmake configuration files")
ELSE ()
   MESSAGE(STATUS "No cmake configuration files detected. Searching for library files manually.")
   IF (CMAKE_SYSTEM_NAME STREQUAL "Linux")
      SET(NFPREFIX_DIR "/usr/local" CACHE PATH "Location of provided NetCDF-Fortran dependencies")
      SET(NFINCLUDE_DIR "/usr/local/include" CACHE PATH "Location of provided netcdf.mod")
      IF (NOT BUILD_SHARED_LIBS)
         SET(NCPREFIX_DIR "/usr/local" CACHE PATH "Location of provided NetCDF-C dependencies")
         SET(H5PREFIX_DIR "/usr/local" CACHE PATH "Location of provided HDF5 dependencies")
         SET(ZPREFIX_DIR  "/usr/local" CACHE PATH "Location of provided zlib dependencies")
      ENDIF ()
   ELSEIF (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      IF (DEFINED ENV{HOMEBREW_PREFIX})
         SET(LIBPREFIX "$ENV{HOMEBREW_PREFIX}")
      ELSE ()
         SET(LIBPREFIX "/usr/local")
      ENDIF()

      SET(NFPREFIX_DIR "${LIBPREFIX}" CACHE PATH "Location of provided NetCDF-Fortran dependencies")
      SET(NFINCLUDE_DIR "${LIBPREFIX}/include" CACHE PATH "Location of provided netcdf.mod")
      IF (NOT BUILD_SHARED_LIBS)
         SET(NCPREFIX_DIR "${LIBPREFIX}" CACHE PATH "Location of provided NetCDF-C dependencies")
         SET(H5PREFIX_DIR "${LIBPREFIX}" CACHE PATH "Location of provided HDF5 dependencies")
         SET(ZPREFIX_DIR  "${LIBPREFIX}" CACHE PATH "Location of provided zlib dependencies")
      ENDIF ()
   ELSEIF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
      FILE(GLOB LIBDIRS "C:/Program Files*/NC4F")
      LIST(SORT LIBDIRS)
      LIST(GET LIBDIRS -1 LIBPREFIX)
      SET(NFPREFIX_DIR "${LIBPREFIX}" CACHE PATH "Location of provided NetCDF-Fortran dependencies")
      SET(NFINCLUDE_DIR "${LIBPREFIX}/include" CACHE PATH "Location of provided netcdf.mod")
      IF (NOT BUILD_SHARED_LIBS)
         # Assumes that the dependency libraries are packaged with NetCDF-C.
         FILE(GLOB LIBDIRS "C:/Program Files*/netCDF*")
         LIST(SORT LIBDIRS)
         LIST(GET LIBDIRS -1 LIBPREFIX)
         SET(NCPREFIX_DIR "${LIBPREFIX}" CACHE PATH "Location of provided NetCDF-C dependencies")
         SET(H5PREFIX_DIR "${LIBPREFIX}" CACHE PATH "Location of provided HDF5 dependencies")
         SET(ZPREFIX_DIR  "${LIBPREFIX}" CACHE PATH "Location of provided zlib dependencies")
      ENDIF ()
   ENDIF ()

   IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
      FIND_FILE(NFBIN
      NAMES nf-config
      HINTS 
         NFPREFIX_DIR
         ENV NETCDF_FORTRAN_HOME
         ENV PATH
      PATH_SUFFIXES
         bin
      )

      IF (NFBIN)
         SET(CMD "${NFBIN}")
         LIST(APPEND CMD "--includedir")
         MESSAGE(STATUS "Searching for NetCDF-Fortran include directory using ${CMD}")
         EXECUTE_PROCESS(COMMAND ${CMD} OUTPUT_VARIABLE NFINCLUDE_DIR ERROR_VARIABLE ERR RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
         IF (NFINCLUDE_DIR)
            MESSAGE(STATUS "Found in ${NFINCLUDE_DIR}")
         ELSE ()
            MESSAGE(STATUS "Cannot execute ${CMD}")
            MESSAGE(STATUS "OUTPUT: ${NFINCLUDE_DIR}")
            MESSAGE(STATUS "RESULT: ${RES}")
            MESSAGE(STATUS "ERROR : ${ERR}")
         ENDIF ()

         SET(CMD "${NFBIN}")
         LIST(APPEND CMD "--prefix")
         MESSAGE(STATUS "Searching for NetCDF-Fortran library directory using ${CMD}")
         EXECUTE_PROCESS(COMMAND ${CMD} OUTPUT_VARIABLE NFPREFIX_DIR ERROR_VARIABLE ERR RESULT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
         IF (NFPREFIX_DIR)
            MESSAGE(STATUS "Found in ${NFPREFIX_DIR}")
         ELSE ()
            MESSAGE(STATUS "Cannot execute ${CMD}")
            MESSAGE(STATUS "OUTPUT: ${NFPREFIX_DIR}")
            MESSAGE(STATUS "RESULT: ${RES}")
            MESSAGE(STATUS "ERROR : ${ERR}")
         ENDIF ()
      ENDIF()
   ENDIF()

   FIND_PATH(NETCDF_FORTRAN_INCLUDE_DIR 
      NAMES netcdf.mod 
      HINTS 
         ${NFINCLUDE_DIR}
         ENV NETCDF_INCLUDE
         ENV NETCDF_FORTRAN_HOME
         ENV CPATH
      PATH_SUFFIXES
         include
         modules
         mod
      REQUIRED
   )

   IF (BUILD_SHARED_LIBS) 
      SET(NETCDFF "netcdff")
   ELSE ()
      IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
         SET(NETCDFF "netcdff.lib")
         SET(NETCDF  "netcdf.lib")
         SET(HDF5    "libhdf5.lib")
         SET(HDF5_HL "libhdf5_hl.lib")
         SET(ZLIB    "zlib.lib")
      ELSE ()
         SET(NETCDFF "libnetcdff.a")
         SET(NETCDF  "libnetcdf.a")
         SET(HDF5    "libhdf5.a")
         SET(HDF5_HL "libhdf5_hl.a")
         SET(ZLIB    "libz.a")
      ENDIF()
   ENDIF()

   FIND_LIBRARY(NFLIB
      NAMES ${NETCDFF} 
      PATHS
         ${NFPREFIX_DIR}
         ENV NETCDF_FORTRAN_HOME
         ENV NETCDF_HOME
         ENV LD_LIBRARY_PATH
      PATH_SUFFIXES
         lib
         ${CMAKE_LIBRARY_ARCHITECTURE} 
      REQUIRED
   )
   ADD_LIBRARY(netCDF::netcdff UNKNOWN IMPORTED PUBLIC)
   IF (CMAKE_SYSTEM_NAME STREQUAL "Windows" AND BUILD_SHARED_LIBS)
      # Get the DLL added in
      FIND_FILE(NFDLL
         NAMES "netcdff.dll"
         HINTS 
            NFPREFIX_DIR
            ENV NETCDF_FORTRAN_HOME
            ENV PATH
         PATH_SUFFIXES
            bin
      )
      SET_TARGET_PROPERTIES(netCDF::netcdff PROPERTIES 
                           IMPORTED_IMPLIB "${NFLIB}"
                           IMPORTED_LOCATION "${NFDLL}"
                           INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_FORTRAN_INCLUDE_DIR}"
                           )
      MESSAGE(STATUS "NetCDF-Fortran dll: ${NFDLL}")
   ELSE ()
      SET_TARGET_PROPERTIES(netCDF::netcdff PROPERTIES 
                           IMPORTED_LOCATION "${NFLIB}"
                           INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_FORTRAN_INCLUDE_DIR}"
                           )
   ENDIF()
   
   MESSAGE(STATUS "NetCDF-Fortran library: ${NFLIB}")
   MESSAGE(STATUS "NetCDF-Fortran include directory: ${NETCDF_FORTRAN_INCLUDE_DIR}")
ENDIF ()


SET(NETCDF_FORTRAN_FOUND TRUE)
MARK_AS_ADVANCED(NETCDF_FORTRAN_LIBRARY NETCDF_FORTRAN_INCLUDE_DIR)