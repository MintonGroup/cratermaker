# Copyright 2023 - David Minton
# This file is part of Cratermaker
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

# Add the source files for the Fortran library
SET(SRC_FILES
    ${SRC}/globals/globals_module.f90
    ${SRC}/bind/bind_module.f90
    ${SRC}/surface/surface_module.f90
    ${SRC}/perlin/perlin_module.f90
    ${SRC}/perlin/perlin_submodule.f90
)

SET(DRIVER_src ${SRC}/main.f90)

# Combine all source files 
SET(CRATERMAKER_src ${SRC_FILES})

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

# Now create the mesh conversion tools borrowed from MPAS
FIND_PACKAGE(HDF5 COMPONENTS C HL REQUIRED)
FIND_PACKAGE(NetCDF REQUIRED)
SET(MESH_CONVERTER mesh_converter_exe)
SET(MESH_CONVERTER_SRC ${SRC}/mesh_conversion_tools_netcdf_c)
SET(MESH_CONVERTER_FILES ${MESH_CONVERTER_SRC}/mpas_mesh_converter.cpp)

ADD_EXECUTABLE(${MESH_CONVERTER} ${MESH_CONVERTER_FILES})
SET_PROPERTY(TARGET ${MESH_CONVERTER} PROPERTY OUTPUT_NAME MpasMeshConverter.x)

TARGET_COMPILE_FEATURES(${MESH_CONVERTER} PRIVATE cxx_std_11)
TARGET_INCLUDE_DIRECTORIES(${MESH_CONVERTER} PRIVATE ${netCDF_INCLUDE_DIR} ${MESH_CONVERTER_SRC})
TARGET_LINK_LIBRARIES(${MESH_CONVERTER} netCDF::netcdf)

SET(CELL_CULLER cell_culler_exe)
SET(CELL_CULLER_SRC ${SRC}/mesh_conversion_tools_netcdf_c)
SET(CELL_CULLER_FILES ${CELL_CULLER_SRC}/mpas_cell_culler.cpp)

ADD_EXECUTABLE(${CELL_CULLER} ${CELL_CULLER_FILES})
SET_PROPERTY(TARGET ${CELL_CULLER} PROPERTY OUTPUT_NAME MpasCellCuller.x)

TARGET_COMPILE_FEATURES(${CELL_CULLER} PRIVATE cxx_std_11)
TARGET_INCLUDE_DIRECTORIES(${CELL_CULLER} PRIVATE ${netCDF_INCLUDE_DIR} ${CELL_CULLER_SRC})
TARGET_LINK_LIBRARIES(${CELL_CULLER} netCDF::netcdf)

# Install all the targets
INSTALL(TARGETS ${CRATERMAKER_LIBRARY} ${MESH_CONVERTER} ${CELL_CULLER} 
    LIBRARY DESTINATION ${INSTALL_LIBDIR}
    RUNTIME DESTINATION ${INSTALL_BINDIR}
    ARCHIVE DESTINATION ${INSTALL_LIBDIR}
    INCLUDES DESTINATION ${INSTALL_INCLUDEDIR}
)
