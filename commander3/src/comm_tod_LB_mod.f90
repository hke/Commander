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
module comm_tod_LB_mod
  !   Module which contains all the LiteBIRD time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_LB_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
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
  implicit none

  private
  public comm_LB_tod

  type, extends(comm_tod) :: comm_LB_tod
   contains
     procedure     :: process_tod          => process_LB_tod
  end type comm_LB_tod

  interface comm_LB_tod
     procedure constructor
  end interface comm_LB_tod

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
    type(comm_params),       intent(in) :: cpar          !comm_param structure, list of all the input parameters
    integer(i4b),            intent(in) :: id_abs        !index of the current band within the parameters 
    class(comm_mapinfo),     target     :: info
    character(len=128),      intent(in) :: tod_type      !
    class(comm_LB_tod),      pointer    :: constructor

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
    
    constructor%xi_n_P_rms      = [-1.0, 0.1, 0.2] ! [sigma0, fknee, alpha]; sigma0 is not used
    if (.true.) then
       constructor%xi_n_nu_fit     = [0.,    0.200] ! More than max(2*fknee_default)
       constructor%xi_n_P_uni(2,:) = [0.001, 0.45]  ! fknee
       constructor%xi_n_P_uni(3,:) = [-2.5, -0.4]   ! alpha
    else
       write(*,*) 'Invalid LiteBIRD frequency label = ', trim(constructor%freq)
       stop
    end if

    ! Initialize common parameters
    call constructor%tod_constructor(cpar, id_abs, info, tod_type)

    ! Initialize instrument-specific parameters
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%nhorn           = 1
    constructor%ndiode          = 1
    constructor%compressed_tod  = .false.
    constructor%correct_sl      = .false.
    constructor%orb_4pi_beam    = .false.
    constructor%symm_flags      = .true.
    constructor%chisq_threshold = 100000000000.d0 !20.d0 ! 9.d0
    constructor%nmaps           = info%nmaps
    constructor%ndet            = num_tokens(cpar%ds_tod_dets(id_abs), ",")

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
    call constructor%initialize_bp_covar(trim(cpar%datadir)//'/'//cpar%ds_tod_bp_init(id_abs))

    ! Construct lookup tables
    call constructor%precompute_lookups()

    ! Load the instrument file
    call constructor%load_instrument_file(nside_beam, nmaps_beam, pol_beam, cpar%comm_chain)

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
  subroutine process_LB_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    ! 
    ! Routine that processes the LiteBIRD time ordered data. 
    ! Samples absolute and relative bandpass, gain and correlated noise in time domain, 
    ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms. 
    ! Writes maps to disc in fits format
    ! 
    ! Arguments:
    ! ----------
    ! self:     pointer of comm_LB_tod class
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
    class(comm_LB_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms

    real(dp)            :: t1, t2
    integer(i4b)        :: i, j, k, l, ierr, ndelta, nside, npix, nmaps
    logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, sample_gain, output_scanlist
    type(comm_binmap)   :: binmap
    type(comm_scandata) :: sd
    character(len=4)    :: ctext, myid_text
    character(len=6)    :: samptext, scantext
    character(len=512)  :: prefix, postfix, prefix4D, filename
    character(len=512), allocatable, dimension(:) :: slist
    real(sp), allocatable, dimension(:)       :: procmask, procmask2
    real(sp), allocatable, dimension(:,:)     :: s_buf
    real(sp), allocatable, dimension(:,:,:)   :: d_calib
    real(sp), allocatable, dimension(:,:,:,:) :: map_sky
    real(dp), allocatable, dimension(:,:)     :: chisq_S, m_buf

    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    ! Toggle optional operations
    sample_rel_bandpass   = .false. !size(delta,3) > 1      ! Sample relative bandpasses if more than one proposal sky
    sample_abs_bandpass   = .false.                ! don't sample absolute bandpasses
    select_data           = self%first_call        ! only perform data selection the first time
    output_scanlist       = mod(iter-1,10) == 0    ! only output scanlist every 10th iteration
    sample_gain           = .false.                ! Gain sampling, LB TOD sims have perfect gain

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
    call distribute_sky_maps(self, map_in, 1.e0, map_sky) ! uK to K

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

    ! Sample gain components in separate TOD loops; marginal with respect to n_corr
     if (sample_gain) then
       ! 'abscal': the global constant gain factor
       call sample_calibration(self, 'abscal', handle, map_sky, procmask, procmask2)
       ! 'relcal': the gain factor that is constant in time but varying between detectors
       call sample_calibration(self, 'relcal', handle, map_sky, procmask, procmask2)
       ! 'deltaG': the time-variable and detector-variable gain
       call sample_calibration(self, 'deltaG', handle, map_sky, procmask, procmask2)
    end if

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
!!$       if (self%output_4D_map > 0) then
!!$          if (mod(iter-1,self%output_4D_map) == 0) then
!!$             call output_4D_maps_hdf(trim(chaindir) // '/tod_4D_chain'//ctext//'_proc' // myid_text // '.h5', &
!!$                  & samptext, self%scanid(i), self%nside, self%npsi, &
!!$                  & self%label, self%horn_id, real(self%polang*180/pi,sp), &
!!$                  & real(self%scans(i)%d%sigma0/self%scans(i)%d%gain,sp), &
!!$                  & sd%pix(:,:,1), sd%psi(:,:,1)-1, d_calib(1,:,:), iand(sd%flag,self%flag0), &
!!$                  & self%scans(i)%d(:)%accept)
!!$          end if
!!$       end if

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
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d0, chisq_S=chisq_S, mask=procmask2)
    else
       call finalize_binned_map(self, binmap, handle, rms_out, 1.d0)
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
    !if (self%output_n_maps > 3) call binmap%outmaps(8)%p%writeFITS(trim(prefix)//'hitmap'//trim(postfix))
    if (self%output_n_maps > 4) call binmap%outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
    if (self%output_n_maps > 5) call binmap%outmaps(5)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
    if (self%output_n_maps > 6) call binmap%outmaps(6)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
    if (self%output_n_maps > 7) call binmap%outmaps(7)%p%writeFITS(trim(prefix)//'zodi'//trim(postfix))

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

  end subroutine process_LB_tod   


end module comm_tod_LB_mod
