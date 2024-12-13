#!/bin/bash
# This script will download the correct OpenMP library for a given MacOS deployment target
# 
# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

# Determine the platform and architecture
OS=$(uname -s)
# If it is not Darwin then exit
if [ $OS != "Darwin" ]; then
   echo "This script is only for MacOS"
   exit 1
fi

SCRIPT_DIR=$(realpath $(dirname $0))
ROOT_DIR=$(realpath ${SCRIPT_DIR}/..)

set -e
cd $ROOT_DIR
. ${SCRIPT_DIR}/set_environment.sh

TARGET_MAJOR=`echo $MACOSX_DEPLOYMENT_TARGET | cut -d. -f1`
TARGET_MINOR=`echo $MACOSX_DEPLOYMENT_TARGET | cut -d. -f2`
TARGET_REV=`echo $MACOSX_DEPLOYMENT_TARGET | cut -d. -f3`

#Figure out which version to get
case $TARGET_MAJOR in
   14|15)
      OMPVER="16.0.4"
      DVER="20"
      ;;
   13)
      OMPVER="14.0.6"
      DVER="20"
      ;;
   12)
      if ((TARGET_MINOR>=5)); then
         OMPVER="14.0.6"
         DVER="20"
      else
         OMPVER="13.0.0"
         DVER="21"
      fi
      ;;
   11)
      if ((TARGET_MINOR>=3)); then
         OMPVER="12.0.1"
         DVER="20"
      else
         OMPVER="11.0.1"
         DVER="20"
      fi
      ;;
   10)
      DVER="17"
      case $TARGET_MINOR in
         15)
            case $TARGET_REV in
               4)
                  OMPVER="10.0.0"
                  ;;
               2)
                  OMPVER="9.0.1"
                  ;;
            esac
            ;;
         14)
            case $TARGET_REV in
               4)
                  OMPVER="8.0.1"
                  ;;
               3)
                  OMPVER="7.1.0"
                  ;;
            esac
            ;;
         *)
            OMPVER="7.1.0"
            ;;
      esac
      ;;
esac

printf "*********************************************************\n"
printf "*             FETCHING OPENMP LIBRARY                   *\n"
printf "*********************************************************\n"
LOMP_DIR="${PREFIX}/../.."
printf "Copying files to ${LOMP_DIR}\n"
mkdir -p ${DEPENDENCY_DIR}

filename="openmp-${OMPVER}-darwin${DVER}-Release.tar.gz"
#Download and install the libraries
printf "Downloading ${filename}\n"
if [ -w "${LOMP_DIR}" ]; then
   curl -L https://mac.r-project.org/openmp/${filename} | tar xvz -C ${LOMP_DIR}
else
   sudo curl -L https://mac.r-project.org/openmp/${filename} | tar xvz -C ${LOMP_DIR}
fi
