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
# Description: This script determines the location of FFTW on the host system.
# If it fails to do so, it will download, compile and install FFTW from source.
#================================================================================

message(STATUS "---------------------------------------------------------------")
# TODO: make it so components will matter because now it install everything because 
# I gave the command to add appropriate configure suboptions to configure command
#if(NOT (FFTW_FORCE_COMPILE OR ALL_FORCE_COMPILE))
if(USE_SYSTEM_FFTW AND USE_SYSTEM_LIBS)
	find_package(FFTW 
		COMPONENTS 
		DOUBLE 
		DOUBLE_THREADS
		FLOAT 
		#FLOAT_MPI 
		FLOAT_OPENMP
		FLOAT_THREADS
		)
endif()
# Is TRUE if one of the components were missing
if(NOT FFTW_FOUND)
	# Configure command to compile FFTW from source
	# DOUBLE (default one)
	list(APPEND fftw_double_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${MPI_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
		"${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw/configure" 
		"--prefix=<INSTALL_DIR>")
	# FLOAT
	list(APPEND fftw_float_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${MPI_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
		"${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw_float/configure" 
		"--prefix=<INSTALL_DIR>")
	# First, we determine which component were found, so we can link them
	# others will be compiled from source
	# double component is default one, so no configuration option required
	if(NOT FFTW_DOUBLE_FOUND)
		message(STATUS "Missing component - DOUBLE - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_DOUBLE_LIB: ${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_DOUBLE_THREADS_FOUND)
		message(STATUS "Missing component - DOUBLE_THREADS - will be compiled from source")	
		list(APPEND fftw_double_configure_command "--enable-threads")
	else()
		message(STATUS "Found FFTW_DOUBLE_THREADS_LIB: ${FFTW_DOUBLE_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		message(STATUS "Missing component - FLOAT - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-float")
	else()
		message(STATUS "Found FFTW_FLOAT_LIB: ${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_THREADS_FOUND)
		message(STATUS "Missing component - FLOAT_THREADS - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-threads")
	else()
		message(STATUS "Found FFTW_FLOAT_THREADS_LIB: ${FFTW_FLOAT_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		message(STATUS "Missing component - FLOAT_OPENMP - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-openmp")
	else()
		message(STATUS "Found FFTW_FLOAT_OPENMP_LIB: ${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_MPI_FOUND)
		message(STATUS "Missing component - FLOAT_MPI - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-mpi")
	else()
		message(STATUS "Found FFTW_FLOAT_MPI_LIB: ${FFTW_FLOAT_MPI_LIB}")
	endif()
	#------------------------------------------------------------------------------
	# Getting FFTW from source
	#------------------------------------------------------------------------------
	# Splitting external project add into 3 steps:
	# 1. To download the project
	# 2. To compile with single and double precision - requires by GNU compilers
	ExternalProject_Add(fftw
		DEPENDS required_libraries
		URL "${fftw_url}"
		URL_MD5 "${fftw_md5}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/fftw"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		# Ommiting Configuration, build and install steps
		CONFIGURE_COMMAND ""
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_float" 
		)

	#------------------------------------------------------------------------------
	# FFTW Float precision
	#------------------------------------------------------------------------------
	ExternalProject_Add(fftw_float
		DEPENDS fftw #${project}_copy_step	
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_float"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_float"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON 
		LOG_BUILD ON 
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${${project}_float_configure_command}"
		)
	
	#------------------------------------------------------------------------------
	# FFTW Double precision
	#------------------------------------------------------------------------------
	ExternalProject_Add(fftw_double
		# specifying that this project depends on the previous one
		DEPENDS fftw #${project}_copy_step
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON 
		LOG_BUILD ON 
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${${project}_double_configure_command}"
		)

	# adding fftw3, fftw3_threads, fftw3_mpi and fftws3_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
	if(NOT FFTW_DOUBLE_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_DOUBLE_THREADS_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3_threads${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_DOUBLE_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_THREADS_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_MPI_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_mpi${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_MPI_LIB}")
	endif()
else()
	# adding empty targets in case FFTW was found on the system
	add_custom_target(fftw ALL "")
	add_custom_target(fftw_double ALL "")
	add_custom_target(fftw_float ALL "")
	set(FFTW3_LIBRARIES
		${FFTW_DOUBLE_LIB}
		${FFTW_DOUBLE_THREADS_LIB}
		${FFTW_FLOAT_LIB}
		${FFTW_FLOAT_OPENMP_LIB}
		${FFTW_FLOAT_THREADS_LIB}
		#${FFTW_FLOAT_MPI_LIB}
		)
endif()
