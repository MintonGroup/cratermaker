# Copyright 2023 - David Minton
# This file is part of Cratermaker
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

# Communicate version number and other CMake build variables to the source code
set(GLOBAL_MODULE_IN ${SRC}/globals/globals_module.f90.in)
set(GLOBAL_MODULE_OUT ${SRC}/globals/globals_module.f90)
CONFIGURE_FILE(${GLOBAL_MODULE_IN} ${GLOBAL_MODULE_OUT})

# Add the source files
SET(SRC_FILES
    ${SRC}/globals/globals_module.f90
    ${SRC}/bind/bind_module.f90
    ${SRC}/simulation/simulation_module.f90
    ${SRC}/target/target_module.f90
    ${SRC}/util/util_module.f90
    ${SRC}/util/util.f90
)

SET(DRIVER_src ${SRC}/main.f90)

# Combine all source files 
set(CRATERMAKER_src ${SRC_FILES})

# Turn preprocessor on for all files
SET_SOURCE_FILES_PROPERTIES(${CRATERMAKER_src} ${DRIVER_src} PROPERTIES Fortran_PREPROCESS ON)

# Set Build Type variable
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "DEBUG")
   ADD_DEFINITIONS(-DDEBUG)
ELSEIF(BT STREQUAL "PROFILE")
   ADD_DEFINITIONS(-DPROFILE)
ENDIF()

#####################################################
# Add the needed libraries 
#####################################################
# Create a library from the source files, except the driver
ADD_LIBRARY(${CRATERMAKER_LIBRARY} ${CRATERMAKER_src})
IF (NOT BUILD_SHARED_LIBS) 
    SET_PROPERTY(TARGET ${CRATERMAKER_LIBRARY} PROPERTY POSITION_INDEPENDENT_CODE)
ENDIF ()


IF(USE_OPENMP OR USE_SIMD)
    SET_PROPERTY(TARGET ${CRATERMAKER_LIBRARY} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET ${CRATERMAKER_LIBRARY} APPEND_STRING PROPERTY LINK_FLAGS  "${OpenMP_Fortran_FLAGS} ")
ENDIF()

IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
   SET_PROPERTY(TARGET ${CRATERMAKER_LIBRARY} APPEND_STRING PROPERTY LINK_FLAGS  "/NODEFAULTLIB")
ENDIF()

# Check to see if the compiler allows for local-spec in do concurrent statements. Set a preprocessor variable if it does
SET(TESTFILE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}")
SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranDoConcurrentLoc.f90")
FILE(WRITE "${TESTFILE}"
"
program TestDoConcurrentLoc
integer :: i
real,dimension(10) :: a
do concurrent(i = 1:10) shared(a)
   a(i) = i 
end do
end program TestDoConcurrentLoc
")
TRY_COMPILE(DOCONLOC_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE} COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
IF (DOCONLOC_WORKS)
    MESSAGE(STATUS "DO CONCURRENT supports locality-spec")
    TARGET_COMPILE_DEFINITIONS(${CRATERMAKER_LIBRARY} PRIVATE -DDOCONLOC)
ELSE ()
    MESSAGE(STATUS "DO CONCURRENT does not support locality-spec")
ENDIF (DOCONLOC_WORKS)

# Check to see if quad precision is supported
SET(TESTFILE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}")
SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranQuadPrecisionReal.f90")
FILE(WRITE "${TESTFILE}"
"
program TestQuadPrecisionReal
integer, parameter :: QP = selected_Real_kind(30) 
real(QP) :: x
end program TestQuadPrecisionReal
")
TRY_COMPILE(QUADPREC ${CMAKE_BINARY_DIR} ${TESTFILE} COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
IF (QUADPREC)
    MESSAGE(STATUS "Quad precision real is supported")
    TARGET_COMPILE_DEFINITIONS(${CRATERMAKER_LIBRARY} PRIVATE -DQUADPREC)
ELSE ()
    MESSAGE(STATUS "Quad precision real is not supported")
ENDIF ()


#####################################
# Tell how to install this executable
#####################################
IF(MSVC)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
    FILE(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX} CMAKE_INSTALL_PREFIX)
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(MSVC)
INSTALL(TARGETS ${CRATERMAKER_LIBRARY}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

