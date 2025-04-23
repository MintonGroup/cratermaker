#!/bin/bash
#
# This gets the current gfortran version number
# 
# Copyright 2025 - The Minton Group at Purdue University
# This file is part of Cratermaker.
# Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Cratermaker. 
# If not, see: https://www.gnu.org/licenses. 

FC="$(command -v gfortran-14 || command -v gfortran-13 || command -v gfortran-12 || command -v gfortran)"
echo "$(${FC} --version | head -n1 | awk '{print $NF}' | cut -d. -f1)"
set +a
