# Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

SET(CRATERMAKER_BINDINGS _bindings)

# Find the Cython executable, but don't look in the cmake root directory (due how cmake and cython are installed on the RCAC system)
FIND_PROGRAM(CYTHON 
            NAMES "cython" 
            NO_CMAKE_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            NO_CMAKE_FIND_ROOT_PATH
            )
MESSAGE(STATUS "Cython executable path: ${CYTHON}")
SET(CYTHON_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${CRATERMAKER_BINDINGS}.pyx" "--output-file" "${CMAKE_CURRENT_BINARY_DIR}/${CRATERMAKER_BINDINGS}.c")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF (BT STREQUAL "DEBUG")
    LIST(APPEND CYTHON_ARGS "--gdb")
endif ()
ADD_CUSTOM_COMMAND(
  OUTPUT "${CRATERMAKER_BINDINGS}.c"
  DEPENDS "${CRATERMAKER_BINDINGS}.pyx"
  VERBATIM
  COMMAND "${CYTHON}" ${CYTHON_ARGS} )

PYTHON_ADD_LIBRARY(${CRATERMAKER_BINDINGS} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${CRATERMAKER_BINDINGS}.c" WITH_SOABI)

IF (NOT BUILD_SHARED_LIBS) 
    SET_PROPERTY(TARGET ${CRATERMAKER_BINDINGS} PROPERTY POSITION_INDEPENDENT_CODE)
ENDIF ()
TARGET_LINK_LIBRARIES(${CRATERMAKER_BINDINGS} PUBLIC ${CRATERMAKER_LIBRARY} netCDF::netcdff HDF5::HDF5)
TARGET_INCLUDE_DIRECTORIES(${CRATERMAKER_BINDINGS} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Define the install locations
INSTALL(TARGETS ${CRATERMAKER_BINDINGS} LIBRARY DESTINATION ${INSTALL_LIBDIR})