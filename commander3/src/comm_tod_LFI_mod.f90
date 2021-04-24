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
module comm_tod_LFI_mod
  !   Module which contains all the LFI time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates
  !       all data needed for TOD processing
  !   process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
  !
  use comm_tod_mod
  use comm_param_mod
  use comm_map_mod
  use comm_conviqt_mod
  use pix_tools
  use healpix_types
  use comm_huffman_mod
  use comm_hdf_mod
  use comm_fft_mod
  use comm_shared_arr_mod
  use spline_1D_mod
  use comm_4D_map_mod
  use comm_tod_driver_mod
  use comm_utils
  use comm_tod_adc_mod
  implicit none

  private
  public comm_LFI_tod

  type, extends(comm_tod) :: comm_LFI_tod
     integer(i4b) :: nbin_spike
     real(dp),          allocatable, dimension(:)       :: mb_eff
     real(dp),          allocatable, dimension(:,:)     :: diode_weights
     type(adc_pointer), allocatable, dimension(:,:,:,:) :: adc_corrections ! ndet, n_diode, nadc_templates, (sky, load)
     real(dp),          allocatable, dimension(:,:)     :: spike_templates ! nbin, ndet
     real(dp),          allocatable, dimension(:,:)     :: spike_amplitude ! nscan, ndet
   contains
     procedure     :: process_tod             => process_LFI_tod
     procedure     :: diode2tod_inst          => diode2tod_LFI
     procedure     :: load_instrument_inst    => load_instrument_LFI
     procedure     :: dumpToHDF_inst          => dumpToHDF_LFI
     procedure     :: construct_corrtemp_inst => construct_corrtemp_LFI
     procedure     :: filter_reference_load
  end type comm_LFI_tod

  interface comm_LFI_tod
     procedure constructor
  end interface comm_LFI_tod

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
    !           Information about the maps for this band, like how the maps are distributed in memory
    ! tod_type: string
    !           Instrument specific tod type
    !
    ! Returns
    ! ----------
    ! constructor: pointer
    !              Pointer that contains all instrument data

    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id_abs
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type
    class(comm_LFI_tod),     pointer    :: constructor

    integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam, ierr
    logical(lgt) :: pol_beam

    ! Allocate object
    allocate(constructor)

    ! Set up noise PSD type and priors
    constructor%freq            = cpar%ds_label(id_abs)
    constructor%n_xi            = 3
    constructor%noise_psd_model = 'oof'
    allocate(constructor%xi_n_P_uni(constructor%n_xi,2))
    allocate(constructor%xi_n_P_rms(constructor%n_xi))
    
    constructor%xi_n_P_rms      = [-1.d0, 0.1d0, 0.2d0] ! [sigma0, fknee, alpha]; sigma0 is not used
    if (trim(constructor%freq) == '030') then
       constructor%xi_n_nu_fit     = [0.d0, 1.225d0]    ! More than max(7*fknee_DPC)
       constructor%xi_n_P_uni(2,:) = [0.010d0, 0.45d0]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-2.5d0, -0.4d0]   ! alpha
    else if (trim(constructor%freq) == '044') then
       constructor%xi_n_nu_fit     = [0.d0, 1.00d0]    ! More than max(2*fknee_DPC)
       constructor%xi_n_P_uni(2,:) = [0.002d0, 0.40d0]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-2.5d0, -0.4d0]   ! alpha
    else if (trim(constructor%freq) == '070') then
       constructor%xi_n_nu_fit     = [0.d0, 0.140d0]    ! More than max(2*fknee_DPC)
       constructor%xi_n_P_uni(2,:) = [0.001d0, 0.25d0]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-3.0d0, -0.4d0]   ! alpha
    else
       write(*,*) 'Invalid LFI frequency label = ', trim(constructor%freq)
       stop
    end if

    ! Initialize instrument-specific parameters
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%nhorn           = 1
    constructor%ndiode          = 4
    constructor%compressed_tod  = .true.
    constructor%correct_sl      = .true.
    constructor%orb_4pi_beam    = .true.
    constructor%symm_flags      = .true.
    constructor%chisq_threshold = 20.d6 ! 9.d0
    constructor%nmaps           = info%nmaps
    constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",")

    nside_beam                  = 512
    nmaps_beam                  = 3
    pol_beam                    = .true.
    constructor%nside_beam      = nside_beam

    ! Initialize common parameters
    call constructor%tod_constructor(cpar, id_abs, info, tod_type)

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

    ! Define diode labels
    do i = 1, constructor%ndet
       if (index(constructor%label(i), 'M') /= 0) then
          constructor%diode_names(i,1) = 'sky00'
          constructor%diode_names(i,2) = 'sky01'
          constructor%diode_names(i,3) = 'ref00'
          constructor%diode_names(i,4) = 'ref01'
       else
          constructor%diode_names(i,1) = 'sky10'
          constructor%diode_names(i,2) = 'sky11'
          constructor%diode_names(i,3) = 'ref10'
          constructor%diode_names(i,4) = 'ref11'
       end if
    end do

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    if (trim(constructor%freq) == '030') then
       constructor%polang = -[-3.428, -3.428, 2.643, 2.643]*pi/180.
    else if (trim(constructor%freq) == '044') then
       constructor%polang = -[-2.180, -2.180,  7.976, 7.976, -4.024, -4.024]*pi/180.
    else if (trim(constructor%freq) == '070') then
       constructor%polang = -[ 0.543, 0.543,  1.366, 1.366,  -1.811, -1.811, -1.045, -1.045,  -2.152, -2.152,  -0.960, -0.960]*pi/180.
    end if

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(cpar%datadir)//'/'//cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    call constructor%precompute_lookups()

    ! allocate LFI specific instrument file data
    constructor%nbin_spike      = nint(constructor%samprate*sqrt(3.d0))
    allocate(constructor%mb_eff(constructor%ndet))
    allocate(constructor%diode_weights(constructor%ndet, 2))
    allocate(constructor%spike_templates(constructor%nbin_spike, constructor%ndet))
    allocate(constructor%spike_amplitude(constructor%nscan,constructor%ndet))
    allocate(constructor%adc_corrections(constructor%ndet, 2, 2, 2))

    ! Load the instrument file
    call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)
    constructor%spike_amplitude = 0.d0

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
  subroutine process_LFI_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    !
    ! Routine that processes the LFI time ordered data.
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
    ! Writes maps to disc in fits format
    !
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_LFI_tod class
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
    class(comm_LFI_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms

    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
    logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, output_scanlist
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext
    character(len=512)  :: prefix, postfix, prefix4D, filename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp), allocatable, dimension(:)       :: procmask, procmask2, sigma0
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

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
    call distribute_sky_maps(self, map_in, 1.e-6, map_sky) ! uK to K

    ! Distribute processing masks
    allocate(m_buf(0:npix-1,nmaps), procmask(0:npix-1), procmask2(0:npix-1))
    call self%procmask%bcast_fullsky_map(m_buf);  procmask  = m_buf(:,1)
    call self%procmask2%bcast_fullsky_map(m_buf); procmask2 = m_buf(:,1)
    deallocate(m_buf)

    ! Precompute far sidelobe Conviqt structures
    if (self%correct_sl) then
       if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
       do i = 1, self%ndet
          !TODO: figure out why this is rotated
          call map_in(i,1)%p%YtW()  ! Compute sky a_lms
          self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, self%slbeam(i)%p%info%nside, &
               & 100, 3, 100, self%slbeam(i)%p, map_in(i,1)%p, 2)
       end do
    end if
!    write(*,*) 'qqq', self%myid
!    if (.true. .or. self%myid == 78) write(*,*) 'a', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
!!$    call mpi_finalize(ierr)
!!$    stop

    call update_status(status, "tod_init")

    !------------------------------------
    ! Perform main sampling steps
    !------------------------------------

    ! Sample 1Hz spikes
    call sample_1Hz_spikes(self, handle, map_sky, procmask, procmask2)

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
    call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2)
    call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
    call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)

    ! Prepare intermediate data structures
    call binmap%init(self, .true., sample_rel_bandpass)
    if (sample_abs_bandpass .or. sample_rel_bandpass) then
       allocate(chisq_S(self%ndet,size(delta,3)))
       chisq_S = 0.d0
    end if
    if (output_scanlist) then
       allocate(slist(self%nscan))
       slist   = ''
    end if

    ! Perform loop over scans
    if (self%myid == 0) write(*,*) '   --> Sampling ncorr, xi_n, maps'
    do i = 1, self%nscan
       
       ! Skip scan if no accepted data
       if (.not. any(self%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (sample_rel_bandpass) then
!          if (.true. .or. self%myid == 78) write(*,*) 'b', self%myid, self%correct_sl, self%ndet, self%slconv(1)%p%psires
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_bp_prop=.true.)
       else if (sample_abs_bandpass) then
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true., init_s_sky_prop=.true.)
       else
          call sd%init_singlehorn(self, i, map_sky, procmask, procmask2, init_s_bp=.true.)
       end if
       allocate(s_buf(sd%ntod,sd%ndet))

       ! Calling Simulation Routine
       if (self%enable_tod_simulations) then
          call simulate_tod(self, i, sd%s_tot, handle)
          call sd%dealloc
          cycle
       end if

       ! Sample correlated noise
       call sample_n_corr(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr, sd%pix(:,:,1), dospike=.true.)

       ! Compute noise spectrum parameters
       call sample_noise_psd(self, sd%tod, handle, i, sd%mask, sd%s_tot, sd%n_corr)

       ! Compute chisquare
       do j = 1, sd%ndet
          if (.not. self%scans(i)%d(j)%accept) cycle
          call self%compute_chisq(i, j, sd%mask(:,j), sd%s_sky(:,j), sd%s_sl(:,j) + sd%s_orb(:,j), sd%n_corr(:,j), sd%tod(:,j))
       end do

       ! Select data
       if (select_data) call remove_bad_data(self, i, sd%flag)

       ! Compute chisquare for bandpass fit
       if (sample_abs_bandpass) call compute_chisq_abs_bp(self, i, sd, chisq_S)

       ! Compute binned map
       allocate(d_calib(self%output_n_maps,sd%ntod, sd%ndet))
       call compute_calibrated_data(self, i, sd, d_calib)
       
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
       call bin_TOD(self, i, sd%pix(:,:,1), sd%psi(:,:,1), sd%flag, d_calib, binmap)

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

    ! Solve for maps
    call syncronize_binmap(binmap, self)
    if (sample_rel_bandpass) then
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d6, chisq_S=chisq_S, mask=procmask2)
    else
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d6)
    end if
    map_out%map = binmap%outmaps(1)%p%map

    ! Sample bandpass parameters
    if (sample_rel_bandpass .or. sample_abs_bandpass) then
       call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
       self%bp_delta = delta(:,:,1)
    end if
   
    ! Output maps to disk
    call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
    call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))
    if (self%output_n_maps > 1) call binmap%outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
    if (self%output_n_maps > 2) call binmap%outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
    if (self%output_n_maps > 3) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
    if (self%output_n_maps > 4) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
    if (self%output_n_maps > 5) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    if (self%output_n_maps > 6) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

    ! Clean up
    call binmap%dealloc()
    if (allocated(slist)) deallocate(slist)
    deallocate(map_sky, procmask, procmask2)
    if (self%correct_sl) then
       do i = 1, self%ndet
          call self%slconv(i)%p%dealloc(); deallocate(self%slconv(i)%p)
       end do
    end if

    ! Parameter to check if this is first time routine has been
    self%first_call = .false.

    call update_status(status, "tod_end"//ctext)

  end subroutine process_LFI_tod
  
  
  subroutine load_instrument_LFI(self, instfile, band)
    !
    ! Reads the LFI specific fields from the instrument file
    ! Implements comm_tod_mod::load_instrument_inst
    !
    ! Arguments:
    !
    ! self : comm_LFI_tod
    !    the LFI tod object (this class)
    ! file : hdf_file
    !    the open file handle for the instrument file
    ! band : int
    !    the index of the current detector
    ! 
    ! Returns : None
    implicit none
    class(comm_LFI_tod),                 intent(inout) :: self
    type(hdf_file),                      intent(in)    :: instfile
    integer(i4b),                        intent(in)    :: band

    integer(i4b) :: i, j
    integer(i4b) :: ext(2)
    real(dp) :: weight
    character(len=2) :: diode_name
    real(dp), dimension(:,:), allocatable :: adc_buffer

    ! Read in mainbeam_eff
    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'mbeam_eff', self%mb_eff(band))

    ! read in the diode weights
    call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'diodeWeight', weight)
    self%diode_weights(band,1:1) = weight
    self%diode_weights(band,2:2) = 1.d0 - weight

    do i=0, 1
      if(index(self%label(band), 'M') /= 0) then
        if(i == 0) then
          diode_name = '00'
        else
          diode_name = '01'
        end if
      else
        if(i == 0) then
          diode_name = '10'
        else
          diode_name = '11'
        end if
      end if
!      write(self%diode_names(band,i+1), 'sky'//diode_name)
!      write(self%diode_names(band,i+3), 'ref'//diode_name)
      ! read in adc correction templates
      call get_size_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'adc91-'//diode_name, ext)
      allocate(adc_buffer(ext(1), ext(2)))
      call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'adc91-'//diode_name, adc_buffer)      
      self%adc_corrections(band, i+1, 1, 1)%p => comm_adc(adc_buffer(:,1), adc_buffer(:,2)) !adc correction for first half, sky
      self%adc_corrections(band, i+1, 1, 2)%p => comm_adc(adc_buffer(:,3), adc_buffer(:,4)) !adc correction for first half, load
      deallocate(adc_buffer)

      call get_size_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'adc953-'//diode_name, ext)
      allocate(adc_buffer(ext(1), ext(2)))
      call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'adc953-'//diode_name, adc_buffer)
      self%adc_corrections(band, i+1, 2, 1)%p => comm_adc(adc_buffer(:,1), adc_buffer(:,2)) !adc correction for second half, sky
      self%adc_corrections(band, i+1, 2, 2)%p => comm_adc(adc_buffer(:,3), adc_buffer(:,4)) !adc corrections for second half, load
      deallocate(adc_buffer)

      if (index(self%label(band), '44') /= 0) then ! read spike templates
         !call read_hdf(instfile, trim(adjustl(self%label(band)))//'/'//'spikes-'//diode_name, self%spike_templates(band, i+1, :, :))
         self%spike_templates = 0.d0
      end if     
 
    end do

  end subroutine load_instrument_LFI
  
  subroutine diode2tod_LFI(self, scan, tod)
    ! 
    ! Generates detector-coadded TOD from low-level diode data
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! scan:     int
    !           Scan ID number
    !
    ! Returns
    ! ----------
    ! tod:      ntod x ndet sp array
    !           Output detector TOD generated from raw diode data
    !
    implicit none
    class(comm_LFI_tod),                 intent(inout) :: self
    integer(i4b),                        intent(in)    :: scan
    real(sp),          dimension(:,:),   intent(out)   :: tod

    integer(i4b) :: i,j,k,half,horn
    real(sp), allocatable, dimension(:,:) :: diode_data, corrected_data, filtered_data

    allocate(diode_data(self%scans(scan)%ntod, self%ndiode))
    allocate(corrected_data(self%scans(scan)%ntod, self%ndiode))
    allocate(filtered_data(self%scans(scan)%ntod, self%ndiode))


    !determine which of the two adc templates we should use
    half = 1
    if (self%scanid(scan) >= 25822) half = 2 !first scan of day 953


    do i=1, self%ndet

        ! Decompress diode TOD for current scan
        call self%decompress_diodes(scan, i, diode_data)

        ! Apply ADC corrections

        do j=1, self%ndiode
          horn=1
          if(index('ref', self%diode_names(i,j)) /= 0) horn=2

          !call self%adc_corrections(i, j, half, horn)%p%adc_correct(diode_data(:,j), corrected_data(:,j))

          !do k = 1, 10
          !   write(*,*) diode_data(k,j), corrected_data(k,j)
          !end do
          !stop

          corrected_data(:,j) = diode_data(:,j)
        end do

        ! Wiener-filter load data         
        call self%filter_reference_load(corrected_data, filtered_data, self%diode_weights(i,:))
        
        ! Compute output differenced TOD

        !w1(sky00 - ref00) + w2(sky01 - ref01)
        tod(:,i) = self%diode_weights(i,1) * (corrected_data(:,1) - corrected_data(:,3)) + self%diode_weights(i,2)*( corrected_data(:,2) - corrected_data(:,4))
        
!!$        open(58,file='comm3_L2fromL1_27M_1670.dat', recl=1024)
!!$        do j = 1, size(tod,1)
!!$           write(58,*) tod(j,1), diode_data(j,:), corrected_data(j,:)
!!$        end do
!!$        close(58)
!!$        stop
        
    end do

    deallocate(diode_data, corrected_data)

  end subroutine diode2tod_LFI

  subroutine filter_reference_load(self, data_in, data_out, weights)
    ! 
    ! Wiener filters the reference load data to minimize noise
    !
    ! Arguments:
    ! ----------
    ! 
    ! self:     comm_tod_LFI object
    !           TOD processing class
    ! data_in:  float array (ntod, ndiode)
    !           input diode timestreams
    !
    ! Returns:
    ! --------
    !
    ! data_out: float array (ntod, ndiode)
    !           filtered output timestreams
    ! weights:  float array (2)
    !           output weights for radiometer 1,2
    implicit none
    class(comm_LFI_tod),               intent(in)   :: self
    real(sp), dimension(:,:),          intent(in)   :: data_in
    real(sp), dimension(:,:),          intent(out)  :: data_out
    real(dp), dimension(:),            intent(inout):: weights    

    integer(i4b) :: i, j, nfft, n
    integer*8    :: plan_fwd, plan_back

    real(dp),     allocatable, dimension(:) :: dt_sky, dt_ref
    real(dp),     allocatable, dimension(:) :: filter
    complex(dpc), allocatable, dimension(:) :: dv_sky, dv_ref

    data_out = data_in

    nfft = size(data_in(:,1))

    allocate(dt_sky(nfft), dt_ref(nfft), dv_sky(0:nfft-1), dv_ref(0:nfft-1), filter(nfft))
    
    call sfftw_plan_dft_r2c_1d(plan_fwd, nfft, dt_ref, dv_ref, fftw_estimate + fftw_unaligned)

    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv_ref, dt_ref, fftw_estimate + fftw_unaligned)

    do i = 1, self%ndiode/2

      ! Check if data is all zeros
      dt_ref = data_in(:, 2*i -1)
      dt_sky = data_in(:, 2*i)
      if(all(dt_ref == 0) .or. all(dt_sky == 0)) return

      ! FFT of ref signal
      call dfftw_execute_dft_r2c(plan_fwd, dt_ref, dv_ref)

      ! FFT of sky signal
      call dfftw_execute_dft_r2c(plan_fwd, dt_sky, dv_sky)     

      if(self%myid == 1) write(*,*) dt_sky(1:100), dv_sky(1:100)
      if(self%myid == 1) write(*,*) dt_sky(nfft-100:nfft-1), dv_sky(nfft-100:nfft-1)
  
 
      ! Compute cross correlation
      filter = 0.5*(dv_sky*conjg(dv_ref) + dv_ref*conjg(dv_sky))/sqrt(dv_sky*conjg(dv_sky) * dv_ref*conjg(dv_ref))

      filter(1) = 1

      do j = 2, size(filter) -1
        if (filter(j+1) > filter(j)) then
          filter(j) = 0.5*(filter(j-1) + filter(j+1))
        end if

        write(*,*) filter(j)

      end do

      stop

      ! Filter ref with cross correlation transfer function
      dv_ref = dv_ref * filter

      ! IFFT ref signal
      call dfftw_execute_dft_c2r(plan_back, dv_ref, dt_ref)

      data_out(:, 2*i-1) = dt_ref/nfft

    end do

  end subroutine filter_reference_load

  subroutine dumpToHDF_LFI(self, chainfile, path)
    ! 
    ! Writes instrument-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_tod)
    !           TOD object
    ! chainfile: derived type (hdf_file)
    !           Already open HDF file handle to existing chainfile
    ! path:   string
    !           HDF path to current dataset, e.g., "000001/tod/030"
    !
    ! Returns
    ! ----------
    ! None
    !
    implicit none
    class(comm_LFI_tod),                 intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:,:) :: amp, amp_tot

    allocate(amp(self%nscan_tot,self%ndet), amp_tot(self%nscan_tot,self%ndet))
    amp = 0.d0
    amp(self%scanid,:) = self%spike_amplitude
    call mpi_reduce(amp, amp_tot, size(amp), MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%info%comm, ierr)

    if (self%myid == 0) then
       call write_hdf(chainfile, trim(adjustl(path))//'1Hz_temp', self%spike_templates)
       call write_hdf(chainfile, trim(adjustl(path))//'1Hz_ampl', amp_tot)
    end if

    deallocate(amp, amp_tot)

  end subroutine dumpToHDF_LFI

  subroutine sample_1Hz_spikes(tod, handle, map_sky, procmask, procmask2)
    !   Sample LFI specific 1Hz spikes shapes and amplitudes
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    !             so that the same sequence can be resumed later on from that same point
    !   map_sky:
    implicit none
    class(comm_LFI_tod),                          intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j, k, b, ierr, nbin
    real(dp)     :: dt, t_tot, t
    real(dp)     :: t1, t2
    character(len=6) :: scantext
    real(dp), allocatable, dimension(:)     :: nval
    real(sp), allocatable, dimension(:)     :: res
    real(dp), allocatable, dimension(:,:)   :: s_sum
    real(dp), allocatable, dimension(:,:,:) :: s_bin
    type(comm_scandata) :: sd

    if (tod%myid == 0) write(*,*) '   --> Sampling 1Hz spikes'

    dt    = 1.d0/tod%samprate   ! Sample time
    t_tot = 1.d0                 ! Time range in sec
    nbin  = tod%nbin_spike      ! Number of bins 

    allocate(s_bin(0:nbin-1,tod%ndet,tod%nscan), s_sum(0:nbin-1,tod%ndet), nval(0:nbin-1))

    ! Compute template per scan
    s_bin = 0.d0
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       tod%apply_inst_corr = .false. ! Disable 1Hz correction for just this call
       call sd%init_singlehorn(tod, i, map_sky, procmask, procmask2)
       tod%apply_inst_corr = .true.  ! Enable 1Hz correction again

       allocate(res(tod%scans(i)%ntod))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle

          res = sd%tod(:,j)/tod%scans(i)%d(j)%gain - (sd%s_sky(:,j) + &
               & sd%s_sl(:,j) + sd%s_orb(:,j))

          nval = 0.d0
          do k = 1, tod%scans(i)%ntod
             if (sd%mask(k,j) == 0.) cycle
             t = modulo(tod%scans(i)%t0(2)/65536.d0 + (k-1)*dt,t_tot)    ! OBT is stored in units of 2**-16 = 1/65536 sec
             b = min(int(t*nbin),nbin-1)
             s_bin(b,j,i) = s_bin(b,j,i)  + res(k)
             nval(b)      = nval(b)       + 1.d0
          end do
          s_bin(:,j,i) = s_bin(:,j,i) / nval
          s_bin(:,j,i) = s_bin(:,j,i) - mean(s_bin(1:nbin/3,j,i))
       end do

!!$       if (trim(tod%freq) == '070') then 
!!$          call int2string(tod%scanid(i),scantext)
!!$          open(58,file='temp_1Hz_22S_PID'//scantext//'.dat')
!!$          do k = 0, nbin-1
!!$             write(58,*) s_bin(k,10,i)
!!$          end do
!!$          close(58)
!!$       end if
!!$
!!$       if (trim(tod%freq) == '044') then 
!!$          call int2string(tod%scanid(i),scantext)
!!$          open(58,file='temp_1Hz_26S_PID'//scantext//'.dat')
!!$          do k = 0, nbin-1
!!$             write(58,*) s_bin(k,6,i)
!!$          end do
!!$          close(58)
!!$       end if

       ! Clean up
        call sd%dealloc
       deallocate(res)
    end do

    ! Compute smoothed templates
    s_sum = 0.d0
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       do j = 1, tod%ndet
          s_sum(:,j) = s_sum(:,j) + s_bin(:,j,i)
       end do
    end do
    call mpi_allreduce(mpi_in_place, s_sum,  size(s_sum),  &
         & MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)

    ! Normalize to maximum of unity
    do j = 1, tod%ndet
       !s_sum(:,j) = s_sum(:,j) - median(s_sum(:,j)) 
       s_sum(:,j) = s_sum(:,j) / maxval(abs(s_sum(:,j))) 
       tod%spike_templates(:,j) = s_sum(:,j) 
    end do

    ! Compute amplitudes per scan and detector
    tod%spike_amplitude = 0.
    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       do j = 1, tod%ndet
          tod%spike_amplitude(i,j) = sum(s_sum(:,j)*s_bin(:,j,i))/sum(s_sum(:,j)**2)
       end do
    end do

    ! Clean up
    deallocate(s_bin, s_sum, nval)

  end subroutine sample_1Hz_spikes

  subroutine construct_corrtemp_LFI(self, scan, pix, psi, s)
    !  Construct an LFI instrument-specific correction template; for now contains 1Hz template only
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  pix: int
    !       index for pixel
    !  psi: int
    !       integer label for polarization angle
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_LFI_tod),                   intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t_tot, t

    dt    = 1.d0/self%samprate   ! Sample time
    t_tot = 1.d0                ! Time range in sec
    nbin  = self%nbin_spike      ! Number of bins 

    do j = 1, self%ndet
       if (.not. self%scans(scan)%d(j)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t = modulo(self%scans(scan)%t0(2)/65536.d0 + (k-1)*dt,t_tot)    ! OBT is stored in units of 2**-16 = 1/65536 sec
          b = min(int(t*nbin),nbin-1)
          s(k,j) = self%spike_amplitude(scan,j) * self%spike_templates(b,j)
       end do
    end do

  end subroutine construct_corrtemp_LFI



end module comm_tod_LFI_mod
