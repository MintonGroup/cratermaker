#!/bin/bash
# Sets environment variables for Linux or Mac. This script is designed so that the build system can be customized by the user simply
# by defining any environment variables prior to running this script.
# 
#
# The following variables are set by this script:
# 
# PREFIX: The installation directory for the dependencies
# NCDIR: The installation directory for the netcdf-c library
# NETCDF_FORTRAN_HOME: The installation directory for the netcdf-fortran library
# NFDIR: Alternate name for the installation directory for the netcdf-fortran library
# NETCDF_FORTRAN_INCLUDE: The include directory for the netcdf-fortran library
# ZLIB_ROOT: The installation directory for the zlib library
# SZIP_ROOT: The installation directory for the szip library
# BZ2_ROOT: The installation directory for the bzip2 library
# ZSTD_ROOT: The installation directory for the zstd library
# HDF5_ROOT: The installation directory for the hdf5 library
# HDF5_LIBDIR: The library directory for the hdf5 library
# HDF5_INCLUDE_DIR: The include directory for the hdf5 library
# HDF5_PLUGIN_PATH: The component directory for the hdf5 library
# HDF5_DIR: The cmake directory for the hdf5 library
# SHTOOLS_HOME: The installation directory for the SHTOOLS library
# LD_LIBRARY_PATH: The library path for the dependencies
# CPATH: The include path for the dependencies
# PATH: The path for the dependency binaries
# CMAKE_INSTALL_LIBDIR: The cmake install directory name (default is "lib"). This is used to override the behavior in Linux systems that use lib64 instead of lib in order that dependency libraries are installed in a consistent location.
# MACOSX_DEPLOYMENT_TARGET: The macOS deployment target (for macOS only)
# HOMEBREW_PREFIX: The homebrew prefix (for macOS only)
# DYLD_LIBRARY_PATH: The dynamic library path (for macOS only)
# LDFLAGS: The linker flags (defined for macOS only)
# CFLAGS: The C compiler flags (defined for macOS only)
# FCFLAGS: The Fortran compiler flags (defined for macOS only)
# FFLAGS: The Fortran 77 compiler flags (defined for macOS only)
# CXXFLAGS: The C++ compiler flags (defined for macOS only)
# CPPFLAGS: The C preprocessor flags (defined for macOS only)
# FC: The Fortran compiler
# F77: The Fortran 77 compiler
# F95: The Fortran 95 compiler
# CC: The C compiler (defined for macOS only)
# CXX: The C++ compiler (defined for macOS only)
# 
# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)
OS=$(uname -s)

set -a
PREFIX="${PREFIX:-/usr/local}"
DEPENDENCY_DIR="${PREFIX}"
LD_LIBRARY_PATH="${PREFIX}/lib"
CPATH="${PREFIX}/include"
PATH="${PREFIX}/bin:${PATH}"
CMAKE_INSTALL_LIBDIR="lib"
NPROC=$(nproc)

FC="$(${SCRIPT_DIR}/get_gfortran_path.sh)"
F77="${FC}"
F95="${FC}"
GFORTRAN_VERSION="$(${SCRIPT_DIR}/get_gfortran_version.sh)"

if [ $OS = "Darwin" ]; then
    CC="/usr/bin/clang"
    CXX="/usr/bin/clang++"
    MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-"$(${SCRIPT_DIR}/get_macosx_deployment_target.sh)"}
    ARCH="$(uname -m)"
    HOMEBREW_PREFIX=${HOMEBREW_PREFIX:-"$(brew --prefix)"}
    SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
    CMAKE_OSX_SYSROOT="${SDKROOT}"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${HOMEBREW_PREFIX}/lib:${HOMEBREW_PREFIX}/lib/gcc/${GFORTRAN_VERSION}"
    DYLD_LIBRARY_PATH="${LD_LIBRARY_PATH}"
    LDFLAGS="-Wl,-no_compact_unwind -L${PREFIX}/lib -L${HOMEBREW_PREFIX}/lib"
    CPATH="${CPATH}:${HOMEBREW_PREFIX}/include"
    CPPFLAGS="-isystem ${PREFIX}/include -Xclang -fopenmp"
    LIBS="-lomp -lquadmath"
    CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} -Wno-deprecated-non-prototype -arch ${ARCH}"
    FCFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
    FFLAGS="${FCFLAGS}"
    CXXFLAGS="${CFLAGS}"
    PATH="${HOMEBREW_PREFIX}/opt/coreutils/libexec/gnubin:${HOMEBREW_PREFIX}/bin:${PATH}"
else
    CC="$(command -v cc)"
    CXX="$(command -v c++)"
    LIBS="-lgomp"
    CFLAGS="-Wa,--noexecstack -fPIC"
    FCFLAGS="${CFLAGS}"
    FFLAGS="${FCFLAGS}"
    CXXFLAGS="${CFLAGS}"
    PATH="${MPI_HOME}/bin:${PATH}"
fi
set +a
