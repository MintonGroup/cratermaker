#!/bin/bash
#
# This gets the current system's MacOS version number and converts it from its full version to just major version. For instance, if
# the MacOS version is 14.4.1 then this script will return 14.0. 
# 
# Copyright 2024 - The Minton Group at Purdue University
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

OS=$(uname -s)
if [ $OS = "Darwin" ]; then
    echo "$(sw_vers -productVersion | cut -d. -f1).0" 
else
    echo ""
fi
set +a
