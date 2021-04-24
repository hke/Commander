!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_tod_WMAP_mod
  !   Module which contains all the WMAP time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
   use comm_tod_mod
   use comm_param_mod
   use comm_map_mod
   use comm_conviqt_mod
   use pix_tools
   use healpix_types
   use comm_huffman_mod
   use comm_hdf_mod
   use comm_fft_mod
   use spline_1D_mod
   use comm_4D_map_mod
   use comm_zodi_mod
   use comm_tod_mapmaking_mod
   use comm_tod_pointing_mod
   use comm_tod_gain_mod
   use comm_tod_bandpass_mod
   use comm_tod_orbdipole_mod
   use comm_tod_driver_mod
   use comm_utils
   implicit none

   private
   public comm_WMAP_tod

   type, extends(comm_tod) :: comm_WMAP_tod
      character(len=20), allocatable, dimension(:) :: labels ! names of fields
   contains
      procedure     :: process_tod           => process_WMAP_tod
   end type comm_WMAP_tod

   interface comm_WMAP_tod
      procedure constructor
   end interface comm_WMAP_tod

contains



   !**************************************************
   !             Constructor
   !**************************************************
   function constructor(cpar, id_abs, info, tod_type)
      !
      ! Constructor function that gathers all the instrument parameters in a pointer
      ! and constructs the objects
      !
      ! Arguments:
      ! ----------
      ! cpar:     derived type
      !           Object containing parameters from the parameterfile.
      ! id_abs:   integer
      !           The index of the current band within the parameters, related to cpar
      ! info:     map_info structure
      !           Information about the maps for this band, e.g.,
      !           how the maps are distributed in memory
      ! tod_type: string
      !           Instrument specific tod type
      !
      ! Returns
      ! ----------
      ! constructor: pointer
      !              Pointer that contains all instrument data

      implicit none
      type(comm_params),      intent(in) :: cpar
      integer(i4b),           intent(in) :: id_abs
      class(comm_mapinfo),    target     :: info
      character(len=128),     intent(in) :: tod_type
      class(comm_WMAP_tod),   pointer    :: constructor

      integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam
      logical(lgt) :: pol_beam

      ! Initialize common parameters
      allocate (constructor)

      ! Set up noise PSD type and priors
      constructor%freq            = cpar%ds_label(id_abs)
      constructor%n_xi            = 3
      constructor%noise_psd_model = 'oof'
      allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
      allocate(constructor%xi_n_P_rms(constructor%n_xi))
    
      constructor%xi_n_P_rms      = [-1.0, 0.1, 0.2]   ! [sigma0, fknee, alpha]; sigma0 is not used
      if (.true.) then
         constructor%xi_n_nu_fit     = [0.0, 0.200]    ! More than max(2*fknee_DPC)
         constructor%xi_n_P_uni(2,:) = [0.001, 0.250]  ! fknee
         constructor%xi_n_P_uni(3,:) = [-3.0, -0.4]    ! alpha
      else
         write(*,*) 'Invalid WMAP frequency label = ', trim(constructor%freq)
         stop
      end if


      call constructor%tod_constructor(cpar, id_abs, info, tod_type)

      ! Set up WMAP specific parameters
      constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
      constructor%nhorn           = 2
      constructor%ndiode          = 1
      constructor%n_xi            = 3
      constructor%compressed_tod  = .true.
      constructor%correct_sl      = .false.
      constructor%orb_4pi_beam    = .false.
      constructor%symm_flags      = .false.
      constructor%chisq_threshold = 400.d0 ! 9.d0
      constructor%nmaps           = info%nmaps
      constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",")
      constructor%verbosity       = cpar%verbosity
      constructor%x_im            = 0d0

      ! Iniitialize TOD labels
      allocate (constructor%labels(7))
      constructor%labels(1) = 'map'
      constructor%labels(2) = 'res'
      constructor%labels(3) = 'ncorr'
      constructor%labels(4) = 'bpcorr'
      constructor%labels(5) = 'orb'
      constructor%labels(6) = 'sl'
      constructor%labels(7) = 'zodi'

      ! Initialize beams
      nside_beam                  = 512
      nmaps_beam                  = 3
      pol_beam                    = .true.
      constructor%nside_beam      = nside_beam


      ! Get detector labels
      call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label)

      ! Define detector partners
      do i = 1, constructor%ndet
         if (mod(i,2) == 1) then
            constructor%partner(i) = i+1
         else
            constructor%partner(i) = i-1
         end if
         constructor%horn_id(i) = (i+1)/2
      end do

      ! Read the actual TOD
      call constructor%read_tod(constructor%label)

      ! Initialize bandpass mean and proposal matrix
      call constructor%initialize_bp_covar(trim(cpar%datadir)//cpar%ds_tod_bp_init(id_abs))

      ! Construct lookup tables
      call constructor%precompute_lookups()

      !load the instrument file
      call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

      ! Need precompute the main beam precomputation for both the A-horn and
      ! B-horn.
      ! Allocate sidelobe convolution data structures
      allocate(constructor%slconv(constructor%ndet), constructor%orb_dp)
      constructor%orb_dp => comm_orbdipole(constructor%mbeam)

      ! Initialize all baseline corrections to zero
      do i = 1, constructor%nscan
         constructor%scans(i)%d%baseline = 0.d0
      end do

   end function constructor

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
      !
      ! Routine that processes the WMAP time ordered data.
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
      ! Writes maps to disc in fits format
      !
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_WMAP_tod class
      !           Points to output of the constructor
      ! chaindir: string
      !           Directory for output files
      ! chain:    integer
      !           Index number of the chain being run
      ! iter:     integer
      !           Gibbs iteration number
      ! handle:   planck_rng derived type
      !           Healpix definition for random number generation
      !           so that the same sequence can be resumed later on from that same point
      ! map_in:   array
      !           Array of dimension (ndet,ndelta) with pointer to maps,
      !           with both access to maps and changing them.
      !           ndet is the number of detectors and
      !           ndelta is the number of bandpass deltas being considered
      ! delta:    array
      !           Array of bandpass corrections with dimensions (0:ndet,npar,ndelta)
      !           where ndet is number of detectors, npar is number of parameters
      !           and ndelta is the number of bandpass deltas being considered
      !
      ! Returns:
      ! ----------
      ! map_out: comm_map class
      !          Final output map after TOD processing combined for all detectors
      ! rms_out: comm_map class
      !          Final output rms map after TOD processing combined for all detectors
      implicit none
      class(comm_WMAP_tod),             intent(inout) :: self
      character(len=*),                    intent(in) :: chaindir
      integer(i4b),                        intent(in) :: chain, iter
      type(planck_rng),                 intent(inout) :: handle
      type(map_ptr), dimension(1:, 1:), intent(inout) :: map_in    ! (ndet,ndelta)
      real(dp),  dimension(0:, 1:, 1:), intent(inout) :: delta     ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map), intent(inout) :: map_out      ! Combined output map
      class(comm_map), intent(inout) :: rms_out      ! Combined output rms

      real(dp)     :: t1, t2
      integer(i4b) :: i, j, k, l, n
      integer(i4b) :: nside, npix, nmaps 
      integer(i4b) :: ierr, ndelta
      real(sp), allocatable, dimension(:, :)          :: s_buf
      real(sp), allocatable, dimension(:, :, :)       :: d_calib
      real(dp), allocatable, dimension(:, :)          :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)          :: M_diag
      real(dp), allocatable, dimension(:, :, :)       :: b_map, b_mono, sys_mono
      character(len=512) :: prefix, postfix
      character(len=2048) :: Sfilename

      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, output_scanlist
      type(comm_scandata) :: sd

      character(len=4)   :: ctext, myid_text
      character(len=6)   :: samptext, scantext
      character(len=512), allocatable, dimension(:) :: slist
      real(sp),       allocatable, dimension(:)     :: procmask, procmask2, sigma0
      real(sp),  allocatable, dimension(:, :, :, :) :: map_sky
      class(map_ptr),     allocatable, dimension(:) :: outmaps

      ! biconjugate gradient-stab parameters
      integer(i4b) :: num_cg_iters
      real(dp) ::  epsil(6)
      real(dp), allocatable, dimension(:, :, :) :: bicg_sol
      real(dp), allocatable, dimension(:)       :: map_full


      call int2string(iter, ctext)
      call update_status(status, "tod_start"//ctext)

      ! Toggle optional operations
      sample_rel_bandpass   = size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
      sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
      select_data           = self%first_call        ! only perform data selection the first time
      output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration

      ! Initialize local variables
      ndelta          = size(delta,3)
      self%n_bp_prop  = ndelta
      nside           = map_out%info%nside
      nmaps           = map_out%info%nmaps
      npix            = 12*nside**2
      self%output_n_maps = 3
      if (self%output_aux_maps > 0) then
         if (mod(iter-1,self%output_aux_maps) == 0) self%output_n_maps = 7
      end if

      call int2string(chain, ctext)
      call int2string(iter, samptext)
      call int2string(self%myid, myid_text)
      prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
      postfix = '_c' // ctext // '_k' // samptext // '.fits'

      ! Distribute maps
      allocate(map_sky(nmaps,self%nobs,0:self%ndet,ndelta))
      !call distribute_sky_maps(self, map_in, 1.e-3, map_sky) ! uK to mK
      call distribute_sky_maps(self, map_in, 1., map_sky) ! K to K?

      ! Distribute processing masks
      allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
      call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
      call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
      deallocate(m_buf)

      ! Allocate total map (for monopole sampling)
      allocate(map_full(0:npix-1))
      map_full = 0.d0


      ! Precompute far sidelobe Conviqt structures
      if (self%correct_sl) then
         do i = 1, self%ndet
            !TODO: figure out why this is rotated
            call map_in(i,1)%p%YtW()  ! Compute sky a_lms
            self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
                 & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
                 & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
         end do
      end if

      call update_status(status, "tod_init")

      !if (self%myid == 0) then
      !  write(*,*) 'Input map statistics, I/Q/U'
      !  write(*,*) minval(map_sky(1,:,:,1)), maxval(map_sky(1,:,:,1)), minval(abs(map_sky(1,:,:,1))), sum(map_sky(1,:,:,1)), sum(abs(map_sky(1,:,:,1)))
      !  write(*,*) minval(map_sky(2,:,:,1)), maxval(map_sky(2,:,:,1)), minval(abs(map_sky(2,:,:,1))), sum(map_sky(2,:,:,1)), sum(abs(map_sky(2,:,:,1)))
      !  write(*,*) minval(map_sky(3,:,:,1)), maxval(map_sky(3,:,:,1)), minval(abs(map_sky(3,:,:,1))), sum(map_sky(3,:,:,1)), sum(abs(map_sky(3,:,:,1)))
      !end if
      !------------------------------------
      ! Perform main sampling steps
      !------------------------------------
      call sample_baseline(self, handle, map_sky, procmask, procmask2)
      !do i = 1, self%nscan
      !  if (self%scanid(i) == 30) then
      !    call sd%init_differential(self, i, map_sky, procmask, procmask2, &
      !      & init_s_bp=.true.)
      !    write(*,*) "S_orb"
      !    write(*,*) sd%tod(1,1), sd%s_orb(1,1)
      !    write(*,*) sd%tod(1,2), sd%s_orb(1,2)
      !    write(*,*) sd%tod(1,3), sd%s_orb(1,3)
      !    write(*,*) sd%tod(1,4), sd%s_orb(1,4)
      !    write(*,*) 'baseline', self%scans(i)%d(1)%baseline
      !    do j = 1, 4
      !      write(*,*) 'j, sum(sd%s_sky(:,j), sum(sd%s_orb(:,j))', j, sum(sd%s_sky(:,j)), sum(sd%s_orb(:,j))
      !    end do
      !    call sd%dealloc
      !  end if
      !end do
      ! The baseline sampling and the orbital dipole template seem to be exactly
      ! the same. Something must be strange with the accumulation step.
      call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2)
      call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
      call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)
      call sample_calibration(self, 'imbal',  handle, map_sky, procmask, procmask2)



      ! Prepare intermediate data structures
      if (sample_abs_bandpass .or. sample_rel_bandpass) then
         allocate(chisq_S(self%ndet,size(delta,3)))
         chisq_S = 0.d0
      end if
      if (output_scanlist) then
         allocate(slist(self%nscan))
         slist   = ''
      end if

      allocate (M_diag(0:npix-1, nmaps+1))
      allocate ( b_map(0:npix-1, nmaps, self%output_n_maps))
      M_diag = 0d0
      b_map = 0d0

      ! Perform loop over scans
      if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
      do i = 1, self%nscan
         


         ! Skip scan if no accepted data
         if (.not. any(self%scans(i)%d%accept)) cycle
         call wall_time(t1)

         ! Prepare data
         if (sample_rel_bandpass) then
            call sd%init_differential(self, i, map_sky, procmask, procmask2, &
              & init_s_bp=.true., init_s_bp_prop=.true.)
         else if (sample_abs_bandpass) then
            call sd%init_differential(self, i, map_sky, procmask, procmask2, &
              & init_s_bp=.true., init_s_sky_prop=.true.)
         else
            call sd%init_differential(self, i, map_sky, procmask, procmask2, &
              & init_s_bp=.true.)
         end if
         allocate(s_buf(sd%ntod,sd%ndet))

         ! Calling Simulation Routine
         ! Not implemented for differential
         !if (self%enable_tod_simulations) then
         !   call simulate_tod(self, i, sd%s_tot, handle)
         !   call sd%dealloc
         !   cycle
         !end if

         ! Sample correlated noise
         call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, &
           & sd%pix(:,1,:), dospike=.false.)

         ! Compute noise spectrum parameters
         call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

         ! Compute chisquare
         do j = 1, sd%ndet
            if (.not. self%scans(i)%d(j)%accept) cycle
            call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), &
              & sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
         end do

         ! Select data
         ! if (select_data) call remove_bad_data(self, i, sd%flag)

         ! Compute chisquare for bandpass fit
         if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

         ! Compute binned map
         allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
         call compute_calibrated_data(self, i, sd, d_calib)
         if (.false. .and. i==1 .and. self%first_call) then
            call int2string(self%scanid(i), scantext)
            if (self%myid == 0 .and. self%verbosity > 0) write(*,*) 'Writing tod to txt'
            do k = 1, self%ndet
               open(78,file=trim(chaindir)//'/tod_'//trim(self%label(k))//'_pid'//scantext//'.dat', recl=1024)
               write(78,*) "# Sample   uncal_TOD (mK)  n_corr (mK) cal_TOD (mK)  skyA (mK)  skyB (mK)"// &
                    & " s_orbA (mK)  s_orbB (mK)  mask, baseline, flag"
               do j = 1, sd%ntod
                  write(78,*) j, sd%tod(j, k), sd%n_corr(j, k), d_calib(1,j,k), &
                   &  sd%s_sky(j,k), sd%s_orb(j,k), sd%mask(j, k), self%scans(i)%d(k)%baseline
               end do
               close(78)
            end do
         end if
         
         ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
         if (self%output_4D_map > 0) then
            if (mod(iter-1,self%output_4D_map) == 0) then
               allocate(sigma0(sd%ndet))
               do j = 1, sd%ndet
                  sigma0(j) = self%scans(i)%d(j)%N_psd%sigma0/self%scans(i)%d(j)%gain
               end do
               call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
                    & samptext, self%scanid(i), self%nside, self%npsi, &
                    & self%label, self%horn_id, real(self%polang*180/pi,sp), sigma0, &
                    & sd%pix(:,:,1), sd%psi(:,:,1)-1, d_calib(1,:,:), iand(sd%flag,self%flag0), &
                    & self%scans(i)%d(:)%accept)
               deallocate(sigma0)
            end if
         end if

         ! Bin TOD
         call bin_differential_TOD(self, d_calib, sd%pix(:,1,:),  &
           & sd%psi(:,1,:), sd%flag(:,1), self%x_im, procmask, b_map, M_diag, i, &
           & .false.)

         ! Update scan list
         call wall_time(t2)
         self%scans(i)%proctime   = self%scans(i)%proctime   + t2-t1
         self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
         if (output_scanlist) then
            write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
                 & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),&
                 & real(self%spinaxis(i,:),sp)
         end if

         ! Clean up
         call sd%dealloc
         deallocate(s_buf, d_calib)

      end do

      if (self%myid == 0) write(*,*) '   --> Finalizing maps, bp'

      ! Output latest scan list with new timing information
      if (output_scanlist) call self%output_scan_list(slist)


      call update_status(status, "Running allreduce on M_diag")
      call mpi_allreduce(mpi_in_place, M_diag, size(M_diag), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
      call update_status(status, "Running allreduce on b")
      call mpi_allreduce(mpi_in_place, b_map, size(b_map), &
           & MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)

      where (M_diag == 0d0)
         M_diag = 1d0
      end where

      ! Conjugate Gradient solution to (P^T Ninv P) m = P^T Ninv d, or Ax = b
      call update_status(status, "Allocating cg arrays")
      allocate (bicg_sol(0:npix-1, nmaps, self%output_n_maps))
      if (self%myid == 0) then 
         bicg_sol = 0.0d0
         epsil(1)   = 1d-10
         epsil(2:6) = 1d-6
         num_cg_iters = 0
         if (self%verbosity > 0) write(*,*) '  Running BiCG'
      end if

      ! Solve for maps
      call update_status(status, "Starting bicg-stab")
      do l=1, self%output_n_maps
         if (self%verbosity > 0 .and. self%myid == 0) then
           write(*,*) '    Solving for ', trim(adjustl(self%labels(l)))
         end if
         call run_bicgstab(self, handle, bicg_sol, npix, nmaps, num_cg_iters, &
                          & epsil(l), procmask, map_full, M_diag, b_map, l)
      end do
      if (self%verbosity > 0 .and. self%myid == 0) write(*,*) '  Finished BiCG'

      call mpi_bcast(bicg_sol, size(bicg_sol),  MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
      call mpi_bcast(num_cg_iters, 1,  MPI_INTEGER, 0, self%info%comm, ierr)
      allocate(outmaps(self%output_n_maps))
      do i = 1, self%output_n_maps
         outmaps(i)%p => comm_map(self%info)
      end do
      do k = 1, self%output_n_maps
         do j = 1, nmaps
            outmaps(k)%p%map(:, j) = bicg_sol(self%info%pix, j, k)
         end do
      end do


      map_out%map = outmaps(1)%p%map
      rms_out%map = M_diag(self%info%pix, 1:nmaps)**-0.5
      call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
      call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
      do n = 2, self%output_n_maps
        call outmaps(n)%p%writeFITS(trim(prefix)//trim(adjustl(self%labels(n)))//trim(postfix))
      end do

      ! Sample bandpass parameters
      if (sample_rel_bandpass .or. sample_abs_bandpass) then
         call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
         self%bp_delta = delta(:,:,1)
      end if

      ! Clean up temporary arrays
      deallocate(procmask, procmask2)
      deallocate(b_map, M_diag)
      deallocate(map_full)
      if (allocated(chisq_S)) deallocate (chisq_S)
      if (allocated(b_mono)) deallocate (b_mono)
      if (allocated(sys_mono)) deallocate (sys_mono)
      if (allocated(slist)) deallocate (slist)

      if (allocated(outmaps)) then
         do i = 1, self%output_n_maps
            call outmaps(i)%p%dealloc
         end do
         deallocate (outmaps)
      end if

      deallocate (map_sky, bicg_sol)

      if (self%correct_sl) then
         do i = 1, self%ndet
            call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
         end do
      end if

      call int2string(iter, ctext)
      call update_status(status, "tod_end"//ctext)

      ! Parameter to check if this is first time routine has been
      self%first_call = .false.

   end subroutine process_WMAP_tod



end module comm_tod_WMAP_mod
