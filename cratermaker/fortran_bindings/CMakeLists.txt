# Copyright 2023 - David Minton
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

SET(REALISTIC_BIND realistic)
SET(CRATER_BIND crater)
SET(EJECTA_BIND ejecta)

# Find the Cython executable, but don't look in the cmake root directory (due how cmake and cython are installed on the RCAC system)
FIND_PROGRAM(CYTHON 
            NAMES "cython" 
            NO_CMAKE_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            NO_CMAKE_FIND_ROOT_PATH
            )
MESSAGE(STATUS "Cython executable path: ${CYTHON}")
SET(REALISTIC_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${REALISTIC_BIND}.pyx" "--output-file" "${CMAKE_CURRENT_BINARY_DIR}/${REALISTIC_BIND}.c")
SET(CRATER_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${CRATER_BIND}.pyx" "--output-file" "${CMAKE_CURRENT_BINARY_DIR}/${CRATER_BIND}.c")
SET(EJECTA_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${EJECTA_BIND}.pyx" "--output-file" "${CMAKE_CURRENT_BINARY_DIR}/${EJECTA_BIND}.c")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF (BT STREQUAL "DEBUG")
    LIST(APPEND REALISTIC_ARGS "--gdb")
    LIST(APPEND CRATER_ARGS "--gdb")
    LIST(APPEND EJECTA_ARGS "--gdb")
endif ()
ADD_CUSTOM_COMMAND(
  OUTPUT "${REALISTIC_BIND}.c"
  DEPENDS "${REALISTIC_BIND}.pyx"
  VERBATIM
  COMMAND "${CYTHON}" ${REALISTIC_ARGS}
)

ADD_CUSTOM_COMMAND(
  OUTPUT "${CRATER_BIND}.c"
  DEPENDS "${CRATER_BIND}.pyx"
  VERBATIM
  COMMAND "${CYTHON}" ${CRATER_ARGS}
)  

ADD_CUSTOM_COMMAND(
    OUTPUT "${EJECTA_BIND}.c"
    DEPENDS "${EJECTA_BIND}.pyx"
    VERBATIM
    COMMAND "${CYTHON}" ${EJECTA_ARGS}
)
    

PYTHON_ADD_LIBRARY(${REALISTIC_BIND} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${REALISTIC_BIND}.c" WITH_SOABI)
PYTHON_ADD_LIBRARY(${CRATER_BIND} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${CRATER_BIND}.c" WITH_SOABI)
PYTHON_ADD_LIBRARY(${EJECTA_BIND} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${EJECTA_BIND}.c" WITH_SOABI)

EXECUTE_PROCESS(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

TARGET_LINK_LIBRARIES(${REALISTIC_BIND} PUBLIC ${CRATERMAKER_LIBRARY})
TARGET_LINK_LIBRARIES(${CRATER_BIND} PUBLIC ${CRATERMAKER_LIBRARY})
TARGET_LINK_LIBRARIES(${EJECTA_BIND} PUBLIC ${CRATERMAKER_LIBRARY})
TARGET_INCLUDE_DIRECTORIES(${REALISTIC_BIND} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${NUMPY_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${CRATER_BIND} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${NUMPY_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${EJECTA_BIND} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${NUMPY_INCLUDE_DIR})

# Define the install locations
INSTALL(TARGETS ${REALISTIC_BIND} ${CRATER_BIND} ${EJECTA_BIND} LIBRARY DESTINATION ${INSTALL_PYPROJ})