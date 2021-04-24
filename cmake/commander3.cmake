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
# Description: This script compiles/installs Commander3 on host system and it also 
# links Commander3 to external dependencies (such as HDF5, CFitsio, HEALPix etc.)
#================================================================================

message(STATUS "---------------------------------------------------------------")

# To resolve conflict with the flags, we need to compile *.cpp
# into a stand alone library and then just link it to rest of 
# the files. However, with CMake it is possible to compile the
# file into object, which will not create *.a files. We go with
# this approach.
add_library(comm_system_backend 
	STATIC 
	${COMMANDER3_SOURCE_DIR}/comm_system_backend.cpp
	)
target_compile_options(comm_system_backend
	PRIVATE
	# setting flags depending on configuration
	#"$<$<CONFIG:Release>:${COMMANDER3_CXX_COMPILER_FLAGS_RELEASE}>"
	#"$<$<CONFIG:Debug>:${COMMANDER3_CXX_COMPILER_FLAGS_DEBUG}>"
	#"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_CXX_COMPILER_FLAGS_RELWITHDEBINFO}>"
	#"$<$<CONFIG:MinSizeRel>:${COMMANDER3_CXX_COMPILER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_CXX_COMPILER_FLAGS}
	)
# installing comm_system_backend as a library
install(TARGETS comm_system_backend ARCHIVE DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#target_link_options(${comm_system_backend}
#	PUBLIC
#	"$<$<CONFIG:Release>:${COMMANDER3_CXX_LINKER_FLAGS_RELEASE}>"
#	"$<$<CONFIG:Debug>:${COMMANDER3_CXX_LINKER_FLAGS_DEBUG}>"
#	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_CXX_LINKER_FLAGS_RELWITHDEBINFO}>"
#	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_CXX_LINKER_FLAGS_MINSIZEREL}>"
#	# setting other compiler dependent flags
#	${COMMANDER3_CXX_LINKER_FLAGS}
#	)

# TODO: add all sources manually instead of this command, as 
# there seems to be a problem with tempita language
#file(GLOB_RECURSE sources *.f90 *.cpp *.f)
set(sources
	${COMMANDER3_SOURCE_DIR}/commander.f90
	${COMMANDER3_SOURCE_DIR}/ARS_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_fft_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_N_QUcov_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_4D_map_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_int_0D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_N_rms_mod.f90
	# TOD processing modules
	${COMMANDER3_SOURCE_DIR}/comm_tod_noise_psd_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_driver_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_mapmaking_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_LFI_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_gain_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_noise_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_orbdipole_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_pointing_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_WMAP_mod.f90
  ${COMMANDER3_SOURCE_DIR}/comm_tod_SPIDER_mod.f90
  ${COMMANDER3_SOURCE_DIR}/comm_tod_LB_mod.f90
  ${COMMANDER3_SOURCE_DIR}/comm_tod_jump_mod.f90
  ${COMMANDER3_SOURCE_DIR}/comm_tod_driver_mod.f90
	# TOD simulations module (and submodules)
	${COMMANDER3_SOURCE_DIR}/comm_tod_simulations_mod.f90
	#
	${COMMANDER3_SOURCE_DIR}/comm_F_int_1D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_output_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_B_bl_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_int_2D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_param_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_beam_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_int_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_physdust_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_B_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_line_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_powlaw_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_bp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_ptsrc_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_zodi_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_bp_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_freefree_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_shared_arr_mod.f90
	${COMMANDER3_SOURCE_DIR}/d1mach.f
	${COMMANDER3_SOURCE_DIR}/comm_chisq_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_gain_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_shared_output_mod.f90
	${COMMANDER3_SOURCE_DIR}/drc3jj.f
	${COMMANDER3_SOURCE_DIR}/comm_Cl_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_hdf_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_signal_mod.f90
	${COMMANDER3_SOURCE_DIR}/hashtbl_4Dmap.f90
	${COMMANDER3_SOURCE_DIR}/comm_cmb_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_spindust2_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/hashtbl.f90
	${COMMANDER3_SOURCE_DIR}/comm_cmb_relquad_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_huffman_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_spindust_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/InvSamp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_line_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_status_mod.f90
	${COMMANDER3_SOURCE_DIR}/locate_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_conviqt_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_map_mod.f90
	${COMMANDER3_SOURCE_DIR}/math_tools.f90
	${COMMANDER3_SOURCE_DIR}/comm_cr_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_MBB_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_system_mod.f90
	${COMMANDER3_SOURCE_DIR}/powell_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_cr_precond_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_md_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_task_mod.f90
	${COMMANDER3_SOURCE_DIR}/sharp.f90
	${COMMANDER3_SOURCE_DIR}/comm_cr_utils.f90
	# MPI (sub)modules
	${COMMANDER3_SOURCE_DIR}/comm_mpi_mod.f90
	# 
	${COMMANDER3_SOURCE_DIR}/comm_template_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/sort_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_data_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_N_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_bandpass_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_defs.f90
	${COMMANDER3_SOURCE_DIR}/comm_noise_mod.f90
	${COMMANDER3_SOURCE_DIR}/spline_1D_mod.f90
	${COMMANDER3_SOURCE_DIR}/spline_2D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_diffuse_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_nonlin_mod.f90
        ${COMMANDER3_SOURCE_DIR}/comm_tod_adc_mod.f90
	#progressbar_mod.f90
	)

# manually setting Fortran compiler flags for different compilers
#if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
	# To prevent error with length and turn type mismatch into a warnings
	# -fallow-argument-mismatch <= is for 10x compilers
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -ffree-line-length-none -fallow-argument-mismatch -std=legacy")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wl -ffree-line-length-none -fallow-argument-mismatch")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fallow-argument-mismatch")
	# -Wno-argument-mismatch <= is th flag for 9x and older compilers
	#	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-argument-mismatch")
	#	message(STATUS "CMAKE_Fortran_FLAGS are: ${CMAKE_Fortran_FLAGS}")
	#elseif("${CMAKE_Fortran_COMPILER}" MATCHES "Intel")
	#	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -qopt-matmul")
	#endif()
# setting executable name
set(commander3 commander3)
add_executable(${commander3} "")
# make sure that commander executable will be built last
add_dependencies(${commander3} ${projects} fftw_float fftw_double)
target_sources(${commander3}
	PUBLIC	
	${sources}
	)
set_property(TARGET ${commander3} PROPERTY ENABLE_EXPORTS TRUE)
set_property(TARGET ${commander3} PROPERTY LINKER_LANGUAGE Fortran)
# adding compiler flags to commander3 target
target_compile_options(${commander3}
	PUBLIC
	# setting flags depending on configuration
	"$<$<CONFIG:Release>:${COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE}>"
	"$<$<CONFIG:Debug>:${COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG}>"
	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO}>"
	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_Fortran_COMPILER_FLAGS}
	)
# adding linker flags to commander3 target
target_link_options(${commander3}
	PUBLIC
	"$<$<CONFIG:Release>:${COMMANDER3_Fortran_LINKER_FLAGS_RELEASE}>"
	"$<$<CONFIG:Debug>:${COMMANDER3_Fortran_LINKER_FLAGS_DEBUG}>"
	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO}>"
	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_Fortran_LINKER_FLAGS}
	)

#message("cmake dl libs are ${CMAKE_DL_LIBS}")
# LINKING ORDER IN LIBRARIES IS IMPORTANT!
# Order is:
# MPI => OpenMP => Blas => LAPACK => HEALPix => 
# CFITSIO => cURL => libm => dl => HDF5 => ZLib
# => FFTW => comm_system_backend
target_link_libraries(${commander3} 
	PRIVATE
	# linking MPI
	MPI::MPI_Fortran
	# linking OpenMP
	OpenMP::OpenMP_Fortran
	# including MKL
	#-qopt-matmul
	${BLAS_LINKER_FLAGS} 
	${BLAS_LIBRARIES}
	${LAPACK_LINKER_FLAGS} 
	${LAPACK_LIBRARIES}
	# including sharp2
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libsharp2.a"
	#"${out_lib_dir}/libsharp2.a"
	#${SHARP2_LIBRARIES}
	# including healpix
	${HEALPIX_LIBRARIES}
	# including cfitsio
	${CFITSIO_LIBRARIES}
	# to avoid error error dlclose@@GLIBC_2.2.5', so 
	# we need to link math library
	#-lm
	${LIBM_LIBRARY}
	# and -ldl (dl library)
	${CMAKE_DL_LIBS}
	# including hdf5 - first fortran and then general
	${HDF5_Fortran_LIBRARIES}
	# hdf5 requires zlib (?), otherwise will get some stupid error
	#"-lz"
	#"/usr/lib64/libz.so"
	#ZLIB::ZLIB
	#zlib_lib
	#${ZLIB_LIBRARIES}
	# adding curl
	#${CURL_LIBRARIES}
	# these are sort of curl dependencies
	#-lcrypto 
	#-lssl
	#CURL::libcurl
	# Including FFTW3
	${FFTW3_LIBRARIES}
	# Linking commander *.cpp file(s)
	comm_system_backend
	)

# installing commander into appropriate folder
#install(TARGETS ${commander3} RUNTIME DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${commander3} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
