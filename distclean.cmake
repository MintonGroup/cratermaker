# Copyright 2025 - David Minton
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

# This CMake script will delete build directories and files to bring the
# package back to it's distribution state

# We want to start from the top of the source dir, so if we are in build
# we want to start one directory up
GET_FILENAME_COMPONENT(BASEDIR ${CMAKE_SOURCE_DIR} NAME)
IF(${BASEDIR} STREQUAL "build")
    SET(TOPDIR "${CMAKE_SOURCE_DIR}/..")
ELSE()
    SET(TOPDIR "${CMAKE_SOURCE_DIR}")
ENDIF()

CMAKE_POLICY(SET CMP0009 NEW) # Don't follow symlinks
MACRO(GET_PARENT_DIRECTORIES search_string return_list grandparents)
    FILE(GLOB_RECURSE new_list ${search_string})
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        # Remove an extra directory component to return grandparent
        IF(${grandparents})
            # Tack on a fake extension to trick CMake into removing a second
            # path component
            SET(dir_path "${dir_path}.tmp")
            GET_FILENAME_COMPONENT(dir_path ${dir_path} PATH)
        ENDIF(${grandparents})
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

SET(DOC_DIR "docs/_build" "docs/_static/fortran_docs")

# Find directories and files that we will want to remove
FILE(GLOB_RECURSE CMAKECACHE "${TOPDIR}/*CMakeCache.txt")
FILE(GLOB_RECURSE CMAKEINSTALL "${TOPDIR}/*cmake_install.cmake"
                               "${TOPDIR}/*install_manifest.txt")
FILE(GLOB_RECURSE CMAKETESTFILES "${TOPDIR}/*CTestTestfile.cmake")
SET(TOPDIRECTORIES "${TOPDIR}/lib" 
                   "${TOPDIR}/libexec"
                   "${TOPDIR}/bin"
                   "${TOPDIR}/include"
                   "${TOPDIR}/share"
)

# CMake has trouble finding directories recursively, so locate these
# files and then save the parent directory of the files
GET_PARENT_DIRECTORIES(Makefile.cmake CMAKEFILES 0)
GET_PARENT_DIRECTORIES(LastTest.log CMAKETESTING 1)

# Place these files and directories into a list
SET(DEL ${TOPDIRECTORIES}
        ${CMAKECACHE}
        ${CMAKEINSTALL}
        ${MAKEFILE}
        ${CMAKEFILES}
        ${CMAKETESTING}
        ${CMAKETESTFILES}
        ${DOC_DIR}
)

# If we are not in the build dir, delete that as well
IF(NOT (${BASEDIR} STREQUAL "build"))
    FILE(GLOB BUILD "${TOPDIR}/build")
    SET(DEL ${DEL} ${BUILD})
ENDIF()

# Loop over the directories and delete each one
FOREACH(D ${DEL})
    IF(EXISTS ${D})
        FILE(REMOVE_RECURSE ${D})
    ENDIF()
ENDFOREACH()
