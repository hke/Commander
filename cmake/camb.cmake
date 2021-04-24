#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
# 
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================
# Author: Maksym Brilenkov
#================================================================================
# Description: This script installs CAMB as Commander3 subproject. 
#================================================================================


#------------------------------------------------------------------------------
# Getting CAMB from source
#------------------------------------------------------------------------------
# We need to build CAMB as only either static or shared library for this to work
ExternalProject_Add(camb
	DEPENDS required_libraries 
					curl
					cfitsio
					healpix	
	GIT_REPOSITORY "${camb_git_url}"
	GIT_TAG "${camb_git_tag}"
	# PREFIX should be present, otherwise it will pull it into "build" dir
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/camb"
	DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
	SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/camb/src/camb"
	INSTALL_DIR "${CMAKE_INSTALL_PREFIX}" 
	LOG_DIR "${CMAKE_LOG_DIR}"
	LOG_DOWNLOAD ON
	LOG_CONFIGURE OFF
	LOG_BUILD OFF
	LOG_INSTALL OFF
	# commadns to build the project
	CMAKE_ARGS
		-DCMAKE_BUILD_TYPE=Release
		# Specifying installations paths for binaries and libraries
		-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		-DBUILD_SHARED_LIBS:BOOL=OFF
		# Specifying compilers
		-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DBUILD_SHARED_LIBS:BOOL=OFF
		# Check submodules during build
		-DGIT_SUBMODULE:BOOL=ON
		# CFitsIO paths
		#-DCFITSIO_FOUND:BOOL=TRUE
		-DCFITSIO_LIBRARIES:FILEPATH=${CFITSIO_LIBRARIES}
		-DCFITSIO_INCLUDE_DIRS:PATH=${CFITSIO_INCLUDE_DIRS}
		# HEALPix Paths
		-DHEALPIX_LIBRARIES:FILEPATH=${HEALPIX_LIBRARIES}
		-DHEALPIX_INCLUDE_DIRS:PATH=${HEALPIX_INCLUDE_DIRS}
	)
