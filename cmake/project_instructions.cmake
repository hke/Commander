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
# Description: This script contains general instructions on how to fetch and build 
# Commander3 and all its dependencies. It is split into three parts, each containing 
# its set of instaructions/variables. It is done for easier maintenance. 
#================================================================================

# TODO:
# [x] Split this file into several different files, containing corresponding instructions (versions, variables, toolchains etc.);
# [ ] Change URL_MD5 to URL_HASH of every project;
# [x] Change compiler variables from list to string (but leave APPEND); <= doesn't work this way
# [ ] Remove include_directory() and use target_include_directory() instead (for commander3 target);
# [x] Add one variable which will force all libraries to be recompiled;

#------------------------------------------------------------------------------
# including compiler definitions
include(project_toolchains)
# including project defined variables
include(project_variables)
#------------------------------------------------------------------------------

#==============================================================================
# MAIN PROJECT DEPENDENCIES
#==============================================================================

unset(projects)
# project names <= order matters
list(APPEND projects 
	required_libraries
	zlib
	libaec
	mbedtls
	libssh2
	curl
	cfitsio
	blas # blas-lapack module 
	##sharp2
	fftw
	hdf5
	doxygen
	healpix
	camb
	commander3
	)
#==============================================================================
# PROJECTS' URL SOURCES, MD5 HASHES AND CONFIGURE COMMANDS
#==============================================================================
# ZLib -- required by HDF5, cURL and others. HDF Group provides this one, but
# we are going to use the one from official website. They should be identical.
# Official website:
# https://zlib.net/
# Mirrors on Sourforge:
# https://sourceforge.net/projects/libpng/files/zlib/1.2.11/
set(zlib_url "https://sourceforge.net/projects/libpng/files/zlib/1.2.11/zlib-1.2.11.tar.gz")
set(zlib_md5 "1c9f62f0778697a09d36121ead88e08e")
#------------------------------------------------------------------------------
# SZip -- required by HDF5. HDF Group provides this one as well and we are going
#------------------------------------------------------------------------------
# to use their implementation as it contains cmake scripts and easier to compile.
# TODO: Add checking for SZip on the system and if not it will install AEC from source. 
# Probably it is a good idea to write your own FindSZIP.cmake?
# To get the code: https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
# Note: from 1.10.7 HDF5 is shipped with Free and Open Source alternative to SZip -- 
# Adaptive Entropy Coding library (libaec). So we are going to use this from now on.
# Official Git Repo is: https://gitlab.dkrz.de/k202009/libaec
# Version to use now: https://gitlab.dkrz.de/k202009/libaec/-/tree/v1.0.4
set(libaec_git_url "https://gitlab.dkrz.de/k202009/libaec.git")
set(libaec_git_tag "0c0453a0e463da9c2183f46d0255f05645e0e5ef")
#------------------------------------------------------------------------------
# MbedTLS -- needed by LibSSH2 and cURL.
#------------------------------------------------------------------------------
# Using git hash instead of tag as described in CMake website.
# mbedtls used:
# https://github.com/ARMmbed/mbedtls/tree/v2.16.9
set(mbedtls_git_url "https://github.com/ARMmbed/mbedtls.git")
set(mbedtls_git_tag "3fac0bae4a50113989b3d015cd2d948f51a6d9ac")
#------------------------------------------------------------------------------
# LibSSH2 -- needed by cURL
#------------------------------------------------------------------------------
# Using git hash instead of tag as described in CMake website.
# libssh2 used:
# https://github.com/libssh2/libssh2/releases/tag/libssh2-1.9.0
set(libssh2_git_url "https://github.com/libssh2/libssh2.git")
set(libssh2_git_tag "42d37aa63129a1b2644bf6495198923534322d64")
#------------------------------------------------------------------------------
# cURL - needed by CFitsio and HEALPix (if cfitsio is built with cURL support)
#------------------------------------------------------------------------------
# need to specify command separately, othewise it won't work
# TODO: switch to github version 7.74 <= probably better for security?
# It seems that 7.74 is much better in terms of CMake support than version 7.69.
#set(curl_url "https://github.com/curl/curl/releases/download/curl-7_69_0/curl-7.69.0.zip")
set(curl_git_url "https://github.com/curl/curl.git")
set(curl_git_tag "e052859759b34d0e05ce0f17244873e5cd7b457b")
#------------------------------------------------------------------------------
# OpenBLAS -  Open Source Implementation of BLAS and LAPACK
#------------------------------------------------------------------------------
set(blas_url "https://github.com/xianyi/OpenBLAS/releases/download/v0.3.12/OpenBLAS-0.3.12.tar.gz")
set(blas_md5 "baf8c58c0ef6ebe0f9eb74a5c4acd662")
#------------------------------------------------------------------------------
# FFTW
#------------------------------------------------------------------------------
#set(fftw_url "http://fftw.org/fftw-3.3.8.tar.gz")
set(fftw_url "http://fftw.org/fftw-3.3.9.tar.gz")
#set(fftw_md5 "8aac833c943d8e90d51b697b27d4384d")
set(fftw_md5 "50145bb68a8510b5d77605f11cadf8dc")
#------------------------------------------------------------------------------
# HDF5
#------------------------------------------------------------------------------
# TODO: Think about inclusion of SZip and whether you need to download CMake version?
set(hdf5_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz")
set(hdf5_md5 "9e22217d22eb568e09f0cc15fb641d7c")
# This version is CMake prepared by HDf Group
#set(hdf5_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/CMake-hdf5-1.12.0.tar.gz")
#set(hdf5_md5 "33ab3d5b9019ca468364d226e0ccdea6")
#------------------------------------------------------------------------------
# LibSharp2
#------------------------------------------------------------------------------
# Note: libsharp2 comes as a native part of Healpix 3.70 and thus, we do not
# need to compile it independently. But, I will leave this just in case we need
# to switch to older versions (for whatever reason).
#set(sharp2_url "https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz")#"https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz") #"https://github.com/Libsharp/libsharp/archive/v1.0.0.tar.gz")
#------------------------------------------------------------------------------
# CFitsio
#------------------------------------------------------------------------------
#set(cfitsio_url "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz")
set(cfitsio_url "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.49.tar.gz")
#------------------------------------------------------------------------------
# HEALPix
#------------------------------------------------------------------------------
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.tar.gz/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.zip/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")#"https://sourceforge.net/projects/healpix/files/latest/download")
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.tar.gz/download")
set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.70/Healpix_3.70_2020Jul23.tar.gz/download")
#set(healpix_md5 "ed7c9a3d7593577628ed1286fa7a9250")
#set(healpix_md5 "540b243406596205a7a82434d99af41e")
#set(healpix_md5 "9b51b2fc919f4e70076d296826eebee0")
set(healpix_md5 "bdcc2a4b1ede3ed5a07be57e4aec01d2")
# this command is for healpix 3.50 and below
#set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#------------------------------------------------------------------------------
# CAMB
#------------------------------------------------------------------------------
# CAMB with custom CMake support
set(camb_git_url "https://github.com/maksymbrl/CAMB.git")
set(camb_git_tag "36cc4e546427b9b30359473be628d506b09eb658")
#------------------------------------------------------------------------------
# Doxygen
#------------------------------------------------------------------------------
# there is some weird errors appearing for doxygen v1.8.17 and above, so will stick with this one
set(doxygen_url "https://github.com/doxygen/doxygen/archive/Release_1_8_16.tar.gz")#"https://github.com/doxygen/doxygen/archive/Release_1_8_18.tar.gz")#"https://github.com/doxygen/doxygen.git")
# compile doxygen with cmake itself, so no need to specify configure options here
set(flex_url "http://sourceforge.net/projects/flex/files/flex-2.5.39.tar.gz/download")#"https://sourceforge.net/projects/flex/files/flex-2.6.0.tar.xz/download")
set(bison_url "http://ftp.gnu.org/gnu/bison/bison-3.6.tar.gz")#"http://ftp.gnu.org/gnu/bison/bison-3.6.2.tar.gz")
#------------------------------------------------------------------------------

# Include all project configuration files
foreach(project ${projects})
	include("${project}")
endforeach()
# TODO: use target_include_directories instead (for each subproject if necessary). 
include_directories(${CMAKE_INSTALL_PREFIX}/include)
message(STATUS "Finished looking for packages.")
message(STATUS "---------------------------------------------------------------")
