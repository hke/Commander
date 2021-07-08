module comm_tod_SPIDER_mod
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
  use comm_zodi_mod
  use comm_tod_mapmaking_mod
  use comm_tod_pointing_mod
  use comm_tod_gain_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  use comm_tod_jump_mod
  implicit none

  private
  public comm_SPIDER_tod

  integer(i4b), parameter :: N_test      = 19
  integer(i4b), parameter :: samp_N      = 1
  integer(i4b), parameter :: prep_G      = 15
  integer(i4b), parameter :: samp_G      = 2
  integer(i4b), parameter :: prep_acal   = 3
  integer(i4b), parameter :: samp_acal   = 4
  integer(i4b), parameter :: prep_rcal   = 18
  integer(i4b), parameter :: samp_rcal   = 19
  integer(i4b), parameter :: prep_relbp  = 5
  integer(i4b), parameter :: prep_absbp  = 16
  integer(i4b), parameter :: samp_bp     = 11
  integer(i4b), parameter :: samp_sl     = 6
  integer(i4b), parameter :: samp_N_par  = 7
  integer(i4b), parameter :: sel_data    = 8
  integer(i4b), parameter :: bin_map     = 9
  integer(i4b), parameter :: calc_chisq  = 10
  integer(i4b), parameter :: output_slist = 12 
  integer(i4b), parameter :: samp_mono   = 13
  integer(i4b), parameter :: sub_sl      = 14
  integer(i4b), parameter :: sub_zodi    = 17
  logical(lgt), dimension(N_test) :: do_oper


  type, extends(comm_tod) :: comm_SPIDER_tod
     !class(orbdipole_pointer), allocatable :: orb_dp !orbital dipole calculator
   contains
     procedure     :: process_tod        => process_SPIDER_tod
     procedure     :: read_tod_inst      => read_tod_inst_SPIDER
     procedure     :: read_scan_inst     => read_scan_inst_SPIDER
     procedure     :: initHDF_inst       => initHDF_SPIDER
     procedure     :: dumpToHDF_inst     => dumpToHDF_SPIDER
  end type comm_SPIDER_tod

  interface comm_SPIDER_tod
     procedure constructor
  end interface comm_SPIDER_tod

contains

  !**************************************************
  !             Constructor
  !**************************************************
  function constructor(cpar, id_abs, info)
    implicit none
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id_abs
    class(comm_mapinfo),     target     :: info
    class(comm_SPIDER_tod),     pointer    :: constructor

    integer(i4b) :: i, j, k, l, nside_beam, lmax_beam, nmaps_beam, ndelta, np_vec, ierr
    real(dp)     :: f_fill, f_fill_lim(3), theta, phi
    character(len=512) :: datadir
    logical(lgt) :: pol_beam
    type(hdf_file) :: h5_file
    integer(i4b), allocatable, dimension(:) :: pix
    real(dp), dimension(3) :: v

    ! Set up SPIDER specific parameters
    allocate(constructor)
    constructor%myid          = cpar%myid_chain
    constructor%comm          = cpar%comm_chain
    constructor%numprocs      = cpar%numprocs_chain
    constructor%myid_shared   = cpar%myid_shared
    constructor%comm_shared   = cpar%comm_shared
    constructor%myid_inter    = cpar%myid_inter
    constructor%comm_inter    = cpar%comm_inter
    constructor%info          => info
    constructor%output_n_maps = 3
    constructor%init_from_HDF = cpar%ds_tod_initHDF(id_abs)
    constructor%freq          = cpar%ds_label(id_abs)
    constructor%operation     = cpar%operation
    constructor%outdir        = cpar%outdir
    constructor%first_call    = .true.
    constructor%first_scan    = cpar%ds_tod_scanrange(id_abs,1)
    constructor%last_scan     = cpar%ds_tod_scanrange(id_abs,2)
    constructor%flag0         = cpar%ds_tod_flag(id_abs)
    constructor%orb_abscal    = cpar%ds_tod_orb_abscal(id_abs)
    constructor%nscan_tot     = cpar%ds_tod_tot_numscan(id_abs)
    constructor%output_4D_map = cpar%output_4D_map_nth_iter
    constructor%subtract_zodi = cpar%include_TOD_zodi
    constructor%central_freq  = cpar%ds_nu_c(id_abs)
    constructor%samprate_lowres = 1.d0  ! Lowres samprate in Hz
    constructor%compressed_tod = .false.

    call mpi_comm_size(cpar%comm_shared, constructor%numprocs_shared, ierr)

    if (constructor%first_scan > constructor%last_scan) then
       write(*,*) 'Error: First scan larger than last scan for ', trim(constructor%freq)
       call mpi_finalize(ierr)
       stop
    end if

    datadir = trim(cpar%datadir)//'/' 
    constructor%filelist    = trim(datadir)//trim(cpar%ds_tod_filelist(id_abs))
    constructor%procmaskf1  = trim(datadir)//trim(cpar%ds_tod_procmask1(id_abs))
    constructor%procmaskf2  = trim(datadir)//trim(cpar%ds_tod_procmask2(id_abs))
    constructor%instfile    = trim(datadir)//trim(cpar%ds_tod_instfile(id_abs))

    call constructor%get_scan_ids(constructor%filelist)

    constructor%procmask => comm_map(constructor%info, constructor%procmaskf1)
    constructor%procmask2 => comm_map(constructor%info, constructor%procmaskf2)
    do i = 0, constructor%info%np-1
       if (any(constructor%procmask%map(i,:) < 0.5d0)) then
          constructor%procmask%map(i,:) = 0.d0
       else
          constructor%procmask%map(i,:) = 1.d0
       end if
       if (any(constructor%procmask2%map(i,:) < 0.5d0)) then
          constructor%procmask2%map(i,:) = 0.d0
       else
          constructor%procmask2%map(i,:) = 1.d0
       end if
    end do


    constructor%nmaps    = info%nmaps
   !  constructor%ndet     = num_tokens(cpar%ds_tod_dets(id_abs), ",") ! Harald
    constructor%ndet     = count_lines(cpar%ds_tod_dets(id_abs),datadir)
    constructor%nhorn    = 1
    allocate(constructor%stokes(constructor%nmaps))
    allocate(constructor%w(constructor%nmaps,constructor%nhorn,constructor%ndet))
    allocate(constructor%label(constructor%ndet))
    allocate(constructor%partner(constructor%ndet))
    allocate(constructor%horn_id(constructor%ndet))
    constructor%stokes = [1,2,3]
    constructor%w      = 1.d0
    !call get_tokens(cpar%ds_tod_dets(id_abs), ",", constructor%label) ! Harald
    call get_det_labels(cpar%ds_tod_dets(id_abs),datadir,constructor%label)
    do i = 1, constructor%ndet
       if (mod(i,2) == 1) then
          constructor%partner(i) = i+1
       else
          constructor%partner(i) = i-1
       end if
       constructor%horn_id(i) = (i+1)/2
    end do


    if (trim(cpar%ds_bpmodel(id_abs)) == 'additive_shift') then
       ndelta = 1
    else if (trim(cpar%ds_bpmodel(id_abs)) == 'powlaw_tilt') then
       ndelta = 1
    else
       write(*,*) 'Unknown bandpass model'
       stop
    end if
    allocate(constructor%bp_delta(0:constructor%ndet,ndelta))
    constructor%bp_delta = 0.d0

    ! Read the actual TOD
    call constructor%read_tod(constructor%label)

    ! Initialize bandpass mean and proposal matrix
    call constructor%initialize_bp_covar(trim(datadir)//cpar%ds_tod_bp_init(id_abs))

    ! Initialize beams
    allocate(constructor%fwhm(constructor%ndet))
    allocate(constructor%elip(constructor%ndet))
    allocate(constructor%psi_ell(constructor%ndet))
    allocate(constructor%mb_eff(constructor%ndet))
    allocate(constructor%nu_c(constructor%ndet))
    
    call open_hdf_file(constructor%instfile, h5_file, 'r')
    nside_beam = 128
    nmaps_beam = 3
    pol_beam   = .true.
    ! call read_hdf(h5_file, trim(adjustl(constructor%label(1)))//'/'//'sllmax', lmax_beam)
    ! constructor%slinfo => comm_mapinfo(cpar%comm_chain, nside_beam, lmax_beam, &
    !      & nmaps_beam, pol_beam)
    ! allocate(constructor%slbeam(constructor%ndet), constructor%slconv(constructor%ndet))
    ! allocate(constructor%orb_dp_s(constructor%ndet, 9))
    ! do i = 1, constructor%ndet
    !    constructor%slbeam(i)%p => comm_map(constructor%slinfo, h5_file, .true., .true., &
    !         & trim(constructor%label(i)))
    !    !Perform integrals for orbital dipole sidelobe correction term
    !    call constructor%slbeam(i)%p%Y()
    !    constructor%orb_dp_s(i, :) = 0.d0
    !    do j = 0, constructor%slbeam(i)%p%info%np-1
    !      call pix2vec_ring(constructor%nside, constructor%slbeam(i)%p%info%pix(j+1),v)
    !      !x
    !      constructor%orb_dp_s(i, 1) = constructor%orb_dp_s(i, 1) + constructor%slbeam(i)%p%map(j, 1) * v(1)
    !      !y 
    !      constructor%orb_dp_s(i, 2) = constructor%orb_dp_s(i, 2) + constructor%slbeam(i)%p%map(j, 1) * v(2)
    !      !z 
    !      constructor%orb_dp_s(i, 3) = constructor%orb_dp_s(i, 3) + constructor%slbeam(i)%p%map(j, 1) * v(3)
    !      !x^2 
    !      constructor%orb_dp_s(i, 4) = constructor%orb_dp_s(i, 4) + constructor%slbeam(i)%p%map(j, 1) * v(1) * v(1)
    !      !2xy 
    !      constructor%orb_dp_s(i, 5) = constructor%orb_dp_s(i, 5) + 2.d0 * constructor%slbeam(i)%p%map(j, 1) * v(1) * v(2)
    !      !2xz 
    !      constructor%orb_dp_s(i, 6) = constructor%orb_dp_s(i, 6) + 2.d0 * constructor%slbeam(i)%p%map(j, 1) * v(1) * v(3)
    !      !y^2 
    !      constructor%orb_dp_s(i, 7) = constructor%orb_dp_s(i, 7) + constructor%slbeam(i)%p%map(j, 1) * v(2) * v(2)
    !      !2yz 
    !      constructor%orb_dp_s(i, 8) = constructor%orb_dp_s(i, 8) + 2.d0 * constructor%slbeam(i)%p%map(j, 1) * v(2) * v(3)
    !      !z^2 
    !      constructor%orb_dp_s(i, 9) = constructor%orb_dp_s(i, 9) + constructor%slbeam(i)%p%map(j, 1) * v(3) * v(3)
    !      !if(constructor%myid == 0) then
    !      !   write(*,*) v(1), v(2), v(3), constructor%slbeam(i)%p%map(j, 1)
    !      !end if
    !    end do 

    !    constructor%orb_dp_s = constructor%orb_dp_s*4*pi/real(constructor%slbeam(i)%p%info%npix)

    ! end do
   
    !call mpi_allreduce(MPI_IN_PLACE, constructor%orb_dp_s, size(constructor%orb_dp_s), MPI_DOUBLE_PRECISION, MPI_SUM, constructor%info%comm, ierr)

!!$    if (constructor%myid == 0) then 
!!$      do i = 1, 9
!!$        write(*,*) constructor%orb_dp_s(1, i)
!!$      end do
!!$    end if
 
    do i = 1, constructor%ndet
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'fwhm', constructor%fwhm(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'elip', constructor%elip(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'psi_ell', constructor%psi_ell(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'mbeam_eff', constructor%mb_eff(i))
       call read_hdf(h5_file, trim(adjustl(constructor%label(i)))//'/'//'centFreq', constructor%nu_c(i))
    end do

    constructor%mb_eff = 1.d0 
    constructor%mb_eff = constructor%mb_eff / mean(constructor%mb_eff)
    constructor%nu_c   = constructor%nu_c * 1d9

    call close_hdf_file(h5_file)
 
    ! Lastly, create a vector pointing table for fast look-up for orbital dipole
    np_vec = 12*constructor%nside**2 !npix
    allocate(constructor%pix2vec(3,0:np_vec-1))
    do i = 0, np_vec-1
       call pix2vec_ring(constructor%nside, i, constructor%pix2vec(:,i))
    end do

    ! Construct observed pixel array
    allocate(constructor%pix2ind(0:12*constructor%nside**2-1))
    constructor%pix2ind = -1
    do i = 1, constructor%nscan
       allocate(pix(constructor%scans(i)%ntod))          
       do j = 1, constructor%ndet
          do l = 1, constructor%nhorn
             call huffman_decode(constructor%scans(i)%hkey, &
                  constructor%scans(i)%d(j)%pix(l)%p, pix)
             constructor%pix2ind(pix(1)) = 1
             do k = 2, constructor%scans(i)%ntod
                pix(k)  = pix(k-1)  + pix(k)
                constructor%pix2ind(pix(k)) = 1
             end do
          end do
       end do
       deallocate(pix)
    end do
    constructor%nobs = count(constructor%pix2ind == 1) 
    allocate(constructor%ind2pix(constructor%nobs))
    allocate(constructor%ind2sl(constructor%nobs))
    allocate(constructor%ind2ang(2,constructor%nobs))
    j = 1
    do i = 0, 12*constructor%nside**2-1
       if (constructor%pix2ind(i) == 1) then
          constructor%ind2pix(j) = i 
          constructor%pix2ind(i) = j
          call pix2ang_ring(constructor%nside, i, theta, phi)
          call ang2pix_ring(nside_beam, theta, phi, constructor%ind2sl(j))
          constructor%ind2ang(:,j) = [theta,phi]
          j = j+1
       end if
    end do
    f_fill = constructor%nobs/(12.*constructor%nside**2)
    call mpi_reduce(f_fill, f_fill_lim(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, constructor%info%comm, ierr)
    call mpi_reduce(f_fill, f_fill_lim(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, constructor%info%comm, ierr)
    call mpi_reduce(f_fill, f_fill_lim(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, constructor%info%comm, ierr)
    if (constructor%myid == 0) then
       write(*,*) '  Min/mean/max TOD-map f_sky = ', real(100*f_fill_lim(1),sp), real(100*f_fill_lim(3)/constructor%info%nprocs,sp), real(100*f_fill_lim(2),sp)
    end if

  end function constructor

  !**************************************************
  !             Driver routine
  !**************************************************
  subroutine process_SPIDER_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
    implicit none
    class(comm_SPIDER_tod),                      intent(inout) :: self
    character(len=*),                         intent(in)    :: chaindir
    integer(i4b),                             intent(in)    :: chain, iter
    type(planck_rng),                         intent(inout) :: handle
    type(map_ptr),       dimension(1:,1:),    intent(inout) :: map_in       ! (ndet,ndelta)
    real(dp),            dimension(0:,1:,1:), intent(inout) :: delta        ! (0:ndet,npar,ndelta) BP corrections
    class(comm_map),                          intent(inout) :: map_out      ! Combined output map
    class(comm_map),                          intent(inout) :: rms_out      ! Combined output rms

    integer(i4b) :: i, j, k, l, start_chunk, end_chunk, chunk_size, ntod, ndet
    integer(i4b) :: nside, npix, nmaps, naccept, nhorn, ntot, ext(2), nscan_tot
    integer(i4b) :: ierr, main_iter, n_main_iter, ndelta, ncol, n_A, nout=1
    real(dp)     :: t1, t2, t3, t4, t5, t6, t7, t8, chisq_threshold
    real(dp)     :: t_tot(22)
    real(sp)     :: inv_gain
    real(sp),     allocatable, dimension(:,:)     :: n_corr, s_sl, s_sky, s_orb, mask,mask2, s_bp
    real(sp),     allocatable, dimension(:,:)     :: s_mono, s_buf, s_tot, s_zodi
    real(sp),     allocatable, dimension(:,:)     :: s_invN, s_lowres
    real(sp),     allocatable, dimension(:,:,:)   :: s_sky_prop, s_bp_prop
    real(sp),     allocatable, dimension(:,:,:)   :: d_calib
    real(dp),     allocatable, dimension(:)       :: A_abscal, b_abscal
    real(dp),     allocatable, dimension(:,:)     :: chisq_S, m_buf
    real(dp),     allocatable, dimension(:,:)     :: A_map, dipole_mod
    real(dp),     allocatable, dimension(:,:,:)   :: b_map, b_mono, sys_mono
    integer(i4b), allocatable, dimension(:,:)     :: flag
    integer(i4b), allocatable, dimension(:,:,:)   :: pix, psi
    logical(lgt)       :: correct_sl
    character(len=512) :: prefix, postfix, prefix4D, filename
    character(len=2048) :: Sfilename
    character(len=4)   :: ctext
    character(len=6)   :: samptext, scantext
    character(len=512), allocatable, dimension(:) :: slist
    type(shared_1d_int) :: sprocmask, sprocmask2
    real(sp),           allocatable, dimension(:,:,:,:) :: map_sky
    type(shared_2d_dp) :: sA_map
    type(shared_3d_dp) :: sb_map, sb_mono
    class(comm_map), pointer :: condmap 
    class(map_ptr), allocatable, dimension(:) :: outmaps

    real(sp),     allocatable, dimension(:,:)     :: tod_gapfill, s_jump
    integer(i4b), allocatable, dimension(:,:)     :: jumps, offset_range
    real(sp),     allocatable, dimension(:)       :: offset_level
    character(len=100)                            :: dir_name, run_label
    integer(i4b)                                  :: i_det
    integer(i4b), allocatable, dimension(:,:)     :: jumpflag_range
    real(dp)  :: time_start, time_end, time_tot, chisq_threshold_negative
    character(len=2)                              :: it_text

    call int2string(iter, it_text)

    if (iter > 1) self%first_call = .false.
    call int2string(iter, ctext)
    call update_status(status, "tod_start"//ctext)

    t_tot   = 0.d0
    call wall_time(t5)


    ! Set up full-sky map structures
    call wall_time(t1)
    correct_sl      = .false.
    chisq_threshold = 0.d0 !1E10 ! 0.d0 !100.d0 !7.d0
    chisq_threshold_negative = -10.d0 !-1E10 !-10.d0 !50.d0
    n_main_iter     = 4
   !  chisq_threshold = 9.d0  !3000.d0
    !this ^ should be 7.d0, is currently 2000 to debug sidelobes
    ndet            = self%ndet
    ndelta          = size(delta,3)
    nside           = map_out%info%nside
    nmaps           = map_out%info%nmaps
    npix            = 12*nside**2
    nscan_tot       = self%nscan_tot
    chunk_size      = npix/self%numprocs_shared
    if (chunk_size*self%numprocs_shared /= npix) chunk_size = chunk_size+1
    allocate(A_abscal(self%ndet), b_abscal(self%ndet))
    allocate(map_sky(nmaps,self%nobs,0:ndet,ndelta)) 
    allocate(chisq_S(ndet,ndelta))
    allocate(slist(self%nscan))
    allocate(dipole_mod(nscan_tot, ndet))
    slist   = ''

    call int2string(chain, ctext)
    call int2string(iter, samptext)
    prefix = trim(chaindir) // '/tod_' // trim(self%freq) // '_'
    postfix = '_c' // ctext // '_k' // samptext // '.fits'

    if (self%myid == 0) then
       do k = 1, ndelta
          write(*,*) 'delta in =', real(delta(:,1,k),sp)
       end do
    end if

    if (.false.) then
       do i = 1, self%ndet
          filename = trim(chaindir) // '/BP_fg_' // trim(self%label(i)) // '_v1.fits'
          call map_in(i,1)%p%writeFITS(filename)
       end do
       deallocate(A_abscal, chisq_S, slist)
       return
    end if

    ! Distribute fullsky maps
    allocate(m_buf(0:npix-1,nmaps))
    do j = 1, ndelta
       do i = 1, self%ndet
          map_in(i,j)%p%map = 1d-6 * map_in(i,j)%p%map ! uK to K
          call map_in(i,j)%p%bcast_fullsky_map(m_buf)
          do k = 1, self%nobs
             map_sky(:,k,i,j) = m_buf(self%ind2pix(k),:) ! uK to K
          end do
       end do
       do k = 1, self%nobs
          do l = 1, nmaps
             map_sky(l,k,0,j) = sum(map_sky(l,k,1:ndet,j))/ndet
          end do
       end do
    end do
    deallocate(m_buf)

    ! Set up shared processing mask
    call init_shared_1d_int(self%myid_shared, self%comm_shared, &
         & self%myid_inter, self%comm_inter, [npix], sprocmask)
    call sync_shared_1d_int_map(sprocmask, self%procmask%info%pix, &
         & nint(self%procmask%map(:,1)))
    call init_shared_1d_int(self%myid_shared, self%comm_shared, &
         & self%myid_inter, self%comm_inter, [npix], sprocmask2)
    call sync_shared_1d_int_map(sprocmask2, self%procmask2%info%pix, &
         & nint(self%procmask2%map(:,1)))
    call wall_time(t2); t_tot(9) = t2-t1


    ! Compute far sidelobe Conviqt structures
    call wall_time(t1)
    if (self%myid == 0) write(*,*) 'Precomputing sidelobe convolved sky'
    do i = 1, self%ndet
       if (.not. correct_sl) exit

       !TODO: figure out why this is rotated
       call map_in(i,1)%p%YtW()  ! Compute sky a_lms
       self%slconv(i)%p => comm_conviqt(self%myid_shared, self%comm_shared, &
            & self%myid_inter, self%comm_inter, 128, 100, 3, 100, &
            & self%slbeam(i)%p, map_in(i,1)%p, 2)
    end do
    call wall_time(t2); t_tot(13) = t2-t1

    call update_status(status, "tod_init")

    ! Compute output map and rms
    call wall_time(t3)
    do_oper             = .true.
    do main_iter = 1, n_main_iter

       call wall_time(t7)
       call update_status(status, "tod_istart")

       if (self%myid == 0) write(*,*) '  Performing main iteration = ', main_iter
          
       ! Select operations for current iteration
       do_oper(samp_acal)    = .false. !(main_iter == n_main_iter-3) .and. .not. self%first_call
       do_oper(samp_rcal)    = .false. !(main_iter == n_main_iter-2) .and. .not. self%first_call
       do_oper(samp_G)       = .false. !(main_iter == n_main_iter-1) .and. .not. self%first_call
       do_oper(samp_N)       = (main_iter >= n_main_iter-0)
       do_oper(samp_N_par)   = do_oper(samp_N)
       do_oper(prep_relbp)   = ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call .and. mod(iter,2) == 0
       do_oper(prep_absbp)   = .false. !ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call .and. mod(iter,2) == 1
       do_oper(samp_bp)      = ndelta > 1 .and. (main_iter == n_main_iter-0) .and. .not. self%first_call
       do_oper(samp_mono)    = .false.  !do_oper(bin_map)             !.and. .not. self%first_call
       do_oper(bin_map)      = (main_iter == n_main_iter  )
       do_oper(sel_data)     = (main_iter == n_main_iter  ) .and. (iter==4) !(iter==2)
       do_oper(calc_chisq)   = (main_iter == n_main_iter  ) 
       do_oper(sub_sl)       = correct_sl
       do_oper(sub_zodi)     = self%subtract_zodi
       do_oper(output_slist) = mod(iter, 1) == 0

       dipole_mod = 0

       ! Perform pre-loop operations
       if (do_oper(bin_map) .or. do_oper(prep_relbp)) then
          if (do_oper(prep_relbp)) then
             ncol = nmaps + ndet - 1
             n_A  = nmaps*(nmaps+1)/2 + 4*(ndet-1)
             nout = self%output_n_maps + ndelta - 1
          else
             ncol = nmaps
             n_A  = nmaps*(nmaps+1)/2
             nout = self%output_n_maps
          end if
          allocate(outmaps(nout)) 
          do i = 1, nout
             outmaps(i)%p => comm_map(map_out%info)
          end do
          if (allocated(A_map)) deallocate(A_map, b_map)
          allocate(A_map(n_A,self%nobs), b_map(nout,ncol,self%nobs))
          A_map = 0.d0; b_map = 0.d0       
          if (sA_map%init) call dealloc_shared_2d_dp(sA_map)
          call init_shared_2d_dp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [n_A,npix], sA_map)
          call mpi_win_fence(0, sA_map%win, ierr)
          if (sA_map%myid_shared == 0) sA_map%a = 0.d0
          call mpi_win_fence(0, sA_map%win, ierr)
          if (sb_map%init) call dealloc_shared_3d_dp(sb_map)
          call init_shared_3d_dp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [nout,ncol,npix], sb_map)
          call mpi_win_fence(0, sb_map%win, ierr)
          if (sb_map%myid_shared == 0) sb_map%a = 0.d0
          call mpi_win_fence(0, sb_map%win, ierr)
       end if

       if (do_oper(samp_mono)) then
          allocate(b_mono(nmaps,self%nobs,ndet), sys_mono(nmaps,nmaps+ndet,0:self%info%np-1))
          b_mono = 0.d0
          call init_shared_3d_dp(self%myid_shared, self%comm_shared, &
               & self%myid_inter, self%comm_inter, [nmaps,npix,ndet], sb_mono)
          call mpi_win_fence(0, sb_mono%win, ierr)
          if (sb_map%myid_shared == 0) sb_mono%a = 0.d0
          call mpi_win_fence(0, sb_mono%win, ierr)
       end if

       if (do_oper(sel_data)) then
          naccept = 0; ntot    = 0
       end if

       if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
          A_abscal = 0.d0; b_abscal = 0.d0
       end if

       if (do_oper(prep_absbp)) then
          chisq_S = 0.d0
       end if

       call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

       ! Perform main analysis loop
       do i = 1, self%nscan

          !call update_status(status, "tod_loop1")
          call wall_time(t7)

          if (.not. any(self%scans(i)%d%accept)) cycle

          ! Short-cuts to local variables
          call wall_time(t1)
          ntod = self%scans(i)%ntod
          ndet = self%ndet
          
          ! Set up local data structure for current scan
          allocate(n_corr(ntod, ndet))                 ! Correlated noise in V
          allocate(s_sl(ntod, ndet))                   ! Sidelobe in uKcm
          allocate(s_sky(ntod, ndet))                  ! Sky signal in uKcmb
          allocate(s_sky_prop(ntod, ndet,2:ndelta))    ! Sky signal in uKcmb
          allocate(s_bp(ntod, ndet))                   ! Signal minus mean
          allocate(s_bp_prop(ntod, ndet, 2:ndelta))    ! Signal minus mean
          allocate(s_orb(ntod, ndet))                  ! Orbital dipole in uKcmb
          if (do_oper(samp_mono)) allocate(s_mono(ntod, ndet))                 ! Monopole correction in uKcmb
          allocate(s_buf(ntod, ndet))                  ! Buffer
          allocate(s_tot(ntod, ndet))                  ! Sum of all sky compnents
          if (do_oper(sub_zodi)) allocate(s_zodi(ntod, ndet))                 ! Zodical light
          allocate(mask(ntod, ndet))                   ! Processing mask in time
          allocate(mask2(ntod, ndet))                  ! Processing mask in time
          allocate(pix(ntod, ndet, nhorn))                    ! Decompressed pointing
          allocate(psi(ntod, ndet, nhorn))                    ! Decompressed pol angle
          allocate(flag(ntod, ndet))                   ! Decompressed flags

          ! Harald's stuff
          allocate(tod_gapfill(ntod,ndet))
          allocate(s_jump(ntod,ndet))
          allocate(jumps(ntod,ndet))
          
          ! Initializing arrays to zero
          call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1
          
          ! --------------------
          ! Analyze current scan
          ! --------------------
         
          ! Decompress pointing, psi and flags for current scan
          call wall_time(t1)
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             call self%decompress_pointing_and_flags(i, j, pix(:,j,:), &
                  & psi(:,j,:), flag(:,j))
          end do
         !  call self%symmetrize_flags(flag)
          !call validate_psi(self%scanid(i), psi)
          call wall_time(t2); t_tot(11) = t_tot(11) + t2-t1
          !call update_status(status, "tod_decomp")
          ! Construct sky signal template
          call wall_time(t1)
!!$ HKE: commented out this
!!$          if (do_oper(bin_map) .or. do_oper(prep_relbp)) then 
!!$             call project_sky(self, map_sky(:,:,:,1), pix(:,:,1), psi(:,:,1), flag, &
!!$                  & sprocmask%a, i, s_sky, mask, s_bp=s_bp)  
!!$          else 
!!$             call project_sky(self, map_sky(:,:,:,1), pix(:,:,1), psi(:,:,1), flag, &
!!$                  & sprocmask%a, i, s_sky, mask)
!!$          end if
!!$          if (do_oper(prep_relbp)) then
!!$             do j = 2, ndelta
!!$                call project_sky(self, map_sky(:,:,:,j), pix(:,:,1), psi(:,:,1), flag, &
!!$                     & sprocmask2%a, i, s_sky_prop(:,:,j), mask2, s_bp=s_bp_prop(:,:,j))  
!!$             end do
!!$          else if (do_oper(prep_absbp)) then
!!$             do j = 2, ndelta
!!$                call project_sky(self, map_sky(:,:,:,j), pix(:,:,1), psi(:,:,1), flag, &
!!$                     & sprocmask2%a, i, s_sky_prop(:,:,j), mask2)  
!!$             end do
!!$          end if
          if (main_iter == 1 .and. self%first_call) then
             do j = 1, ndet
                if (all(mask(:,j) == 0)) self%scans(i)%d(j)%accept = .false.
                if (self%scans(i)%d(j)%N_psd%sigma0 <= 0.d0) self%scans(i)%d(j)%accept = .false.
             end do
          end if
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             if (self%scans(i)%d(j)%N_psd%sigma0 <= 0) write(*,*) main_iter, self%scanid(i), j, self%scans(i)%d(j)%N_psd%sigma0
          end do
          call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1
          !call update_status(status, "tod_project")
          






          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          ! Harald
          dir_name = 'chains11_test/'
          run_label = 'cut'
          do j=1, ndet
            call wall_time(time_start)
            if ((sum(flag(:,j)) > 0.8*ntod) .or. (.not. self%scans(i)%d(j)%accept)) then
               self%scans(i)%d(j)%accept = .false.
               cycle
            end if
            
            ! Retrieve offsets from previous run, if they exist
            if (allocated(self%scans(i)%d(j)%offset_range)) then
               call expand_offset_list(              &
                  & self%scans(i)%d(j)%offset_range, &
                  & self%scans(i)%d(j)%offset_level, &
                  & s_jump(:,j))
            else
               s_jump(:,j) = 0
            end if

            ! Retrieve jump flags from previous run, if they exist
            if (allocated(self%scans(i)%d(j)%jumpflag_range)) then
               call add_jumpflags(                      &
                  & self%scans(i)%d(j)%jumpflag_range,  &
                  & flag(:,j))
            end if

            ! Scanning for jumps
            if (main_iter==1) then
               call jump_scan(                                         &
                  & self%scans(i)%d(j)%tod - s_sky(:,j) - s_jump(:,j), &
                  & flag(:,j),                                         & 
                  & jumps(:,j),                                        &
                  & offset_range,                                      &
                  & offset_level,                                      &
                  & handle,                                            &
                  & jumpflag_range)

               ! Add offsets to persistent list
               if (.not. allocated(self%scans(i)%d(j)%offset_range)) then
                  allocate(self%scans(i)%d(j)%offset_range(size(offset_level),2))
                  allocate(self%scans(i)%d(j)%offset_level(size(offset_level)))
               
                  self%scans(i)%d(j)%offset_range = offset_range
                  self%scans(i)%d(j)%offset_level = offset_level
               else
                  call update_offset_list(              &
                     & offset_range,                    &
                     & offset_level,                    &
                     & self%scans(i)%d(j)%offset_range, &
                     & self%scans(i)%d(j)%offset_level)
               end if

               ! Add jump flags to persistent list
               if (allocated(jumpflag_range)) then
                  if (.not. allocated(self%scans(i)%d(j)%jumpflag_range)) then
                     allocate(self%scans(i)%d(j)%jumpflag_range(size(jumpflag_range)/2,2))
                     self%scans(i)%d(j)%jumpflag_range = jumpflag_range
                  else
                     call update_jumpflag(jumpflag_range, self%scans(i)%d(j)%jumpflag_range)
                  end if
               end if
            
               call expand_offset_list(               &
                  & self%scans(i)%d(j)%offset_range,  &
                  & self%scans(i)%d(j)%offset_level,  & 
                  & s_jump(:,j))
            end if


               
            call gap_fill_linear(                      &
               & self%scans(i)%d(j)%tod - s_jump(:,j), &
               & flag(:,j),                            &
               & tod_gapfill(:,j),                     &
               & handle,                               &
               & .true.)


            ! call wall_time(time_end)
            ! time_tot = time_end - time_start
            ! write(*,*) "time", time_tot

            ! Output tods
            ! if (main_iter==n_main_iter) then
            if (.false.) then
               call tod2file(                                                           &
                  & trim(dir_name)//'raw_tod_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & self%scans(i)%d(j)%tod)
                  
               call tod2file(                                                          &
                  & trim(dir_name)//'s_jump_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & s_jump(:,j))

               call tod2file(                                                               &
                  & trim(dir_name)//'tod_gapfill_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & tod_gapfill(:,j))
 
               call tod2file(                                                             &
                  & trim(dir_name)//'jump_flag_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & jumps(:,j))

               call tod2file(                                                        &
                  & trim(dir_name)//'flag_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & flag(:,j))  

               call tod2file(                                                         &
                  & trim(dir_name)//'s_sky_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & s_sky(:,j))

               call tod2file(                                                         &
                  & trim(dir_name)//'jumps_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & jumps(:,j))
            end if

            if (allocated(jumpflag_range)) deallocate(jumpflag_range)
            if (allocated(offset_level))   deallocate(offset_level)
            if (allocated(offset_range))   deallocate(offset_range)
            !if (main_iter==n_main_iter) then
            !   deallocate(self%scans(i)%d(j)%offset_range)
            !   deallocate(self%scans(i)%d(j)%offset_level)
            !   if (allocated(self%scans(i)%d(j)%jumpflag_range)) deallocate(self%scans(i)%d(j)%jumpflag_range)
            !end if
          end do































          ! Construct orbital dipole template
          call wall_time(t1)
         !  call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix, psi, s_orb) ! Harald
          s_orb = 0 ! Harald
          call wall_time(t2); t_tot(2) = t_tot(2) + t2-t1
          !call update_status(status, "tod_orb")

          ! Construct zodical light template
          if (do_oper(sub_zodi)) then
             call compute_zodi_template(self%nside, pix(:,:,1), self%scans(i)%satpos, [30.d9, 30.d9, 30.d9, 30.d9], s_zodi)
          end if

          ! Construct sidelobe template 
          call wall_time(t1)
          if (do_oper(sub_sl)) then
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                call self%construct_sl_template(self%slconv(j)%p, &
                     & pix(:,j,1), psi(:,j,1), s_sl(:,j), self%polang(j))
                s_sl(:,j) = 2 * s_sl(:,j) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
             end do
          else
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_sl(:,j) = 0.
             end do
          end if
          call wall_time(t2); t_tot(12) = t_tot(12) + t2-t1

          ! Construct monopole correction template 
          call wall_time(t1)
          if (do_oper(samp_mono)) then
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                !s_mono(:,j) = self%mono(j)
                s_mono(:,j) = 0.d0 ! Disabled for now
             end do
          end if
          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             s_tot(:,j) = s_sky(:,j) + s_sl(:,j) + s_orb(:,j)
             if (do_oper(samp_mono)) s_tot(:,j) = s_tot(:,j) + s_mono(:,j)
          end do
          call wall_time(t2); t_tot(1) = t_tot(1) + t2-t1

          ! Precompute filtered signal for calibration
          if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. do_oper(samp_acal)) then
             call self%downsample_tod(s_orb(:,1), ext)
             allocate(s_invN(ext(1):ext(2), ndet))      ! s * invN
             allocate(s_lowres(ext(1):ext(2), ndet))      ! s * invN
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. .not. self%orb_abscal) then
                   s_buf(:,j) = s_tot(:,j)
                   call fill_all_masked(s_buf(:,j), mask(:,j), ntod, trim(self%operation)=='sample', real(self%scans(i)%d(j)%N_psd%sigma0, sp), handle, self%scans(i)%chunk_num)
                   call self%downsample_tod(s_buf(:,j), ext, &
                        & s_lowres(:,j))!, mask(:,j))
                else
                   call self%downsample_tod(s_orb(:,j), ext, &
                        & s_lowres(:,j))!, mask(:,j))
                end if
             end do
             s_invN = s_lowres
             call multiply_inv_N(self, i, s_invN,   sampfreq=self%samprate_lowres, pow=0.5d0)
             call multiply_inv_N(self, i, s_lowres, sampfreq=self%samprate_lowres, pow=0.5d0)
          end if

          ! Prepare for absolute calibration
          if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle

                if (do_oper(samp_acal)) then
                   if (self%orb_abscal) then
                      s_buf(:, j) = real(self%gain0(0),sp) * (s_tot(:, j) - s_orb(:, j)) + &
                           & real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                   else
                      s_buf(:, j) = real(self%gain0(j) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                   end if
                else
!                   s_buf(:,j) =  s_tot(:,j) - s_orb(:,j) !s_sky(:,j) + s_sl(:,j) + s_mono(:,j)
                   s_buf(:,j) = real(self%gain0(0) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                end if
             end do
             call accumulate_abscal(self, i, mask, s_buf, s_invN, A_abscal, b_abscal, handle, .false.)

             if (.false.) then
                call int2string(self%scanid(i), scantext)
                open(78,file='tod_'//trim(self%label(j))//'_pid'//scantext//'_k'//samptext//'.dat', recl=1024)
                write(78,*) "# Sample     Data (V)    Res (V)    s_sub (K)   s_orb (K)   mask"
                do k = 1, ntod, 60
                   write(78,*) k, mean(1.d0*self%scans(i)%d(j)%tod(k:k+59)), mean(1.d0*self%scans(i)%d(j)%tod(k:k+59) - &
                        & self%scans(i)%d(j)%gain*s_buf(k:k+59,j)), mean(1.d0*s_orb(k:k+59,j)),  mean(1.d0*s_buf(k:k+59,j)),  minval(mask(k:k+59,j))
                end do
                close(78)
             end if
             call wall_time(t2); t_tot(14) = t_tot(14) + t2-t1
          end if

          ! Fit gain 
          if (do_oper(samp_G)) then
             call wall_time(t1)
             call calculate_gain_mean_std_per_scan(self, i, s_invN, mask, s_tot, handle)
             call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
          end if


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Harald
          if (self%first_call) then
            do j = 1, ndet
               self%scans(i)%d(j)%N_psd%sigma0 = 0.0018
            end do
          end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ! Fit correlated noise
          if (do_oper(samp_N)) then
            call wall_time(t1)
            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               if (do_oper(samp_mono)) then
                  s_buf(:,j) = s_tot(:,j)-s_mono(:,j)
               else
                  s_buf(:,j) = s_tot(:,j)
               end if
            end do
            !call sample_n_corr(self, handle, i, mask, s_buf, n_corr, pix(:,:,1), tod_gapfill)
            call wall_time(t2); t_tot(3) = t_tot(3) + t2-t1
          else
            n_corr = 0.
          end if



          ! Use n_corr to find jumps (Harald)
          if (do_oper(samp_N) .and. (iter==1) .and. .true.) then
            do j = 1, ndet
               ! call tod2file(trim(dir_name)//'res_test_'//trim(it_text)//'.txt', tod_gapfill(:,j)-s_sky(:,j))
               if (.not. self%scans(i)%d(j)%accept) cycle
               call jump_scan_stage2( &
                  & n_corr(:,j),    &
                  & flag(:,j),      &
                  & jumps(:,j),     &
                  & offset_range,   &
                  & offset_level,   & 
                  & handle,         & 
                  & jumpflag_range, &
                  & it_text,        &
                  & dir_name)

               call expand_offset_list( &
                  & offset_range,       &
                  & offset_level,       &
                  & s_jump(:,j))

               ! call tod2file(trim(dir_name)//'s_jump_test_'//trim(it_text)//'.txt', s_jump(:,j))
               ! call tod2file(trim(dir_name)//'flag_after_test_'//trim(it_text)//'.txt', flag(:,j))

               ! Add offsets to persistent list
               if (.not. allocated(self%scans(i)%d(j)%offset_range)) then
                  allocate(self%scans(i)%d(j)%offset_range(size(offset_level),2))
                  allocate(self%scans(i)%d(j)%offset_level(size(offset_level)))
               
                  self%scans(i)%d(j)%offset_range = offset_range
                  self%scans(i)%d(j)%offset_level = offset_level
               else
                  call update_offset_list(              &
                     & offset_range,                    &
                     & offset_level,                    &
                     & self%scans(i)%d(j)%offset_range, &
                     & self%scans(i)%d(j)%offset_level)
               end if

               call expand_offset_list(              &
                  & self%scans(i)%d(j)%offset_range, &
                  & self%scans(i)%d(j)%offset_level, &
                  & s_jump(:,j))

               ! Add jump flags to persistent list
               if (allocated(jumpflag_range)) then
                  if (.not. allocated(self%scans(i)%d(j)%jumpflag_range)) then
                     allocate(self%scans(i)%d(j)%jumpflag_range(size(jumpflag_range)/2,2))
                     self%scans(i)%d(j)%jumpflag_range = jumpflag_range
                  else
                     call update_jumpflag( &
                        & jumpflag_range,  &
                        & self%scans(i)%d(j)%jumpflag_range)
                  end if
               end if

               
               if (allocated(jumpflag_range)) deallocate(jumpflag_range)
               if (allocated(offset_level))   deallocate(offset_level)
               if (allocated(offset_range))   deallocate(offset_range)
               
            end do



         end if
         
         
         ! Compute noise spectrum
         if (do_oper(samp_N_par)) then
            call wall_time(t1)
            !call sample_noise_psd(self, handle, i, mask, s_tot, n_corr, tod_gapfill)
            call wall_time(t2); t_tot(6) = t_tot(6) + t2-t1
         end if


         
         ! Harald
         ! if (main_iter==n_main_iter) then
         if (.false.) then
            do j=1, ndet
               call tod2file(                                                          &
                  & trim(dir_name)//'n_corr_'//trim(self%scans(i)%d(j)%label)//'_'//trim(run_label)//'_'//trim(it_text)//'.txt', &
                  & n_corr(:,j))
            end do
         end if

         ! Compute chisquare
         if (do_oper(calc_chisq)) then
            call wall_time(t1)
            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               s_buf(:,j) =  s_sl(:,j) + s_orb(:,j) 
                if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
                call self%compute_chisq(i, j, 1.0-flag(:,j), s_sky(:,j), &
                     & s_buf(:,j), n_corr(:,j), s_jump(:,j))

                call write2file(trim(chaindir)//'/chisq_'//trim(run_label)//'.txt',  iter, self%scans(i)%d(j)%chisq) 
               !  call write2file(trim(dir_name)//'sigma0_'//trim(run_label)//'.txt', iter, self%scans(i)%d(j)%sigma0) 
               !  call write2file(trim(dir_name)//'fknee_'//trim(run_label)//'.txt',  iter, self%scans(i)%d(j)%fknee) 
               !  call write2file(trim(dir_name)//'alpha_'//trim(run_label)//'.txt',  iter, self%scans(i)%d(j)%alpha) 

                
             end do
             call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
          end if
          
          ! Select data
          if (do_oper(sel_data)) then
             call wall_time(t1)
             do j = 1, ndet
                ntot= ntot + 1
                if (.not. self%scans(i)%d(j)%accept) cycle
               !  if (count(iand(flag(:,j),self%flag0) .ne. 0) > 0.1*ntod) then  ! Harald
                if (sum(flag(:,j)) > 0.8*ntod) then
                   self%scans(i)%d(j)%accept = .false.
                else if ((self%scans(i)%d(j)%chisq) > chisq_threshold .or. &
                & isNaN(self%scans(i)%d(j)%chisq) .or. (self%scans(i)%d(j)%chisq) < chisq_threshold_negative) then
                   write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
                        & self%scanid(i), j, ', chisq = ', &
                        & self%scans(i)%d(j)%chisq
                   self%scans(i)%d(j)%accept = .false.
                   cycle
                !else
                   !write(*,fmt='(a,i8,i5,a,f12.1)') 'Accept scan, det = ', &
                   !   & self%scanid(i), j, ', chisq = ', &
                   !   & self%scans(i)%d(j)%chisq
                   !naccept = naccept + 1
                end if
             end do
            !  if (any(.not. self%scans(i)%d%accept)) self%scans(i)%d%accept = .false. ! Harald
            !  do j = 1, ndet
            !     if (.not. self%scans(i)%d(j)%accept) then
            !        self%scans(i)%d(self%partner(j))%accept = .false.
            !     end if
            !  end do
             call wall_time(t2); t_tot(15) = t_tot(15) + t2-t1
          end if

          ! Compute chisquare for bandpass fit
          if (do_oper(prep_absbp)) then
             call wall_time(t1)
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                s_buf(:,j) =  s_sl(:,j) + s_orb(:,j)
                if (do_oper(samp_mono)) s_buf(:,j) =  s_buf(:,j) + s_mono(:,j)
                call self%compute_chisq(i, j, mask2(:,j), s_sky(:,j), &
                     & s_buf(:,j), n_corr(:,j), absbp=.true.)
                chisq_S(j,1) = chisq_S(j,1) + self%scans(i)%d(j)%chisq_prop
                do k = 2, ndelta
                   call self%compute_chisq(i, j, mask2(:,j), s_sky_prop(:,j,k), &
                        & s_buf(:,j), n_corr(:,j), absbp=.true.)
                   chisq_S(j,k) = chisq_S(j,k) + self%scans(i)%d(j)%chisq_prop
                end do
             end do
             call wall_time(t2); t_tot(7) = t_tot(7) + t2-t1
          end if


          ! Compute binned map 
          if (do_oper(bin_map)) then
             call wall_time(t1)
             allocate(d_calib(nout,ntod, ndet)) 
             do j = 1, ndet
                if (.not. self%scans(i)%d(j)%accept) cycle
                inv_gain = 1.0 / real(self%scans(i)%d(j)%gain,sp)
                d_calib(1,:,j) = (self%scans(i)%d(j)%tod - n_corr(:,j) - s_jump(:,j)) * &
                     & inv_gain - s_tot(:,j) + s_sky(:,j) - s_bp(:,j)
                if (do_oper(bin_map) .and. nout > 1) d_calib(2,:,j) = d_calib(1,:,j) - s_sky(:,j) + s_bp(:,j) ! Residual
                if (do_oper(bin_map) .and. nout > 2) d_calib(3,:,j) = (n_corr(:,j) - sum(n_corr(:,j)/ntod)) * inv_gain
                if (do_oper(bin_map) .and. nout > 3) d_calib(4,:,j) = s_bp(:,j)
                if (do_oper(bin_map) .and. do_oper(samp_mono) .and. nout > 4) d_calib(5,:,j) = s_mono(:,j)
                if (do_oper(bin_map) .and. nout > 5) d_calib(6,:,j) = s_orb(:,j)
                if (do_oper(bin_map) .and. nout > 6) d_calib(7,:,j) = s_sl(:,j)
                if (do_oper(prep_relbp)) then
                   do k = 2, ndelta
                      d_calib(self%output_n_maps+k-1,:,j) = d_calib(1,:,j) + s_bp(:,j) - s_bp_prop(:,j,k)
                   end do
                end if
             end do
            

             
             if (do_oper(bin_map) .and. self%output_4D_map > 0 .and. mod(iter,self%output_4D_map) == 0) then

                ! Output 4D map; note that psi is zero-base in 4D maps, and one-base in Commander
!!$                call int2string(self%scanid(i), scantext)
!!$                prefix4D = "!"//trim(prefix) // '4D_pid' // scantext
!!$                call output_4D_maps(prefix4D, postfix, self%scanid(i), self%nside, self%npsi, &
!!$                     & self%label, self%horn_id, real(self%polang*180/pi,sp), &
!!$                     & real(self%scans(i)%d%N_psdsigma0/self%scans(i)%d%gain,sp), &
!!$                     & pix(:,:,1), psi(:,:,1)-1, d_calib(1,:,:), iand(flag,self%flag0), &
!!$                     & self%scans(i)%d(:)%accept)
             end if


             call wall_time(t2); t_tot(5) = t_tot(5) + t2-t1

             if (.false. .and. self%scanid(i) == 1000 .and. do_oper(bin_map) ) then
                call int2string(self%scanid(i), scantext)
                do k = 1, self%ndet
                   open(78,file='tod_'//trim(self%label(k))//'_pid'//scantext//'.dat', recl=1024)
                   write(78,*) "# Sample     Data (V)     Mask    cal_TOD (K)   res (K)"// &
                        & "   n_corr (K)   s_corr (K)   s_mono (K)   s_orb  (K)   s_sl (K)"
                   do j = 1, ntod
                      write(78,*) j, self%scans(i)%d(k)%tod(j), mask(j,1), d_calib(:,j,1)
                   end do
                   close(78)
                end do
             end if


             call wall_time(t1)

!!$ HKE: commented out this
!!$             if (do_oper(samp_mono)) then
!!$                call bin_TOD(self, d_calib, pix(:,:,1), &
!!$                     & psi(:,:,1), flag, A_map, b_map, i, do_oper(prep_relbp), b_mono=b_mono)
!!$             else
!!$                call bin_TOD(self, d_calib, pix(:,:,1), &
!!$                     & psi(:,:,1), flag, A_map, b_map, i, do_oper(prep_relbp))
!!$             end if
             deallocate(d_calib)
             call wall_time(t2); t_tot(8) = t_tot(8) + t2-t1
          end if

          do j = 1, ndet
             if (.not. self%scans(i)%d(j)%accept) cycle
             dipole_mod(self%scanid(i), j) = masked_variance(s_sky(:, j), mask(:, j))
          end do


          ! Clean up
          call wall_time(t1)
          deallocate(n_corr, s_sl, s_sky, s_orb, s_tot, s_buf)
          deallocate(s_bp, s_sky_prop, s_bp_prop)
          deallocate(mask, mask2, pix, psi, flag)
          deallocate(tod_gapfill, s_jump, jumps)
          if (allocated(s_lowres)) deallocate(s_lowres)
          if (allocated(s_invN)) deallocate(s_invN)
          if (do_oper(sub_zodi)) deallocate(s_zodi)
          if (do_oper(samp_mono)) deallocate(s_mono)
          call wall_time(t2); t_tot(18) = t_tot(18) + t2-t1

          call wall_time(t8); t_tot(19) = t_tot(19) + t8-t7

          if (do_oper(output_slist)) then
             self%scans(i)%proctime   = self%scans(i)%proctime   + t8-t7
             self%scans(i)%n_proctime = self%scans(i)%n_proctime + 1
             if (main_iter == n_main_iter) then
                write(slist(i),*) self%scanid(i), '"',trim(self%hdfname(i)), &
                     & '"', real(self%scans(i)%proctime/self%scans(i)%n_proctime,sp),real(self%spinaxis(i,:),sp)
             end if
          end if
          !call update_status(status, "tod_loop2")

       end do

       call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)

!       if (self%myid == 0) then
!          write(*, *) "CHECKPOINT"
!          write(*, *) dipole_mod
!       end if

       if (do_oper(samp_acal)) then
          call wall_time(t1)
          call sample_abscal_from_orbital(self, handle, A_abscal, b_abscal)
          call wall_time(t2); t_tot(16) = t_tot(16) + t2-t1
       end if

       if (do_oper(samp_rcal)) then
          call wall_time(t1)
          call sample_relcal(self, handle, A_abscal, b_abscal)
          call wall_time(t2); t_tot(16) = t_tot(16) + t2-t1
       end if

       if (do_oper(samp_G)) then
          call wall_time(t1)
          call sample_smooth_gain(self, handle, dipole_mod)
          call wall_time(t2); t_tot(4) = t_tot(4) + t2-t1
       end if

       call wall_time(t7)

    end do
    call wall_time(t4)

    ! Output latest scan list with new timing information
    if (do_oper(output_slist)) then
       call update_status(status, "scanlist1")
       call wall_time(t1)
       call self%output_scan_list(slist)
       call wall_time(t2); t_tot(20) = t_tot(20) + t2-t1
       call update_status(status, "scanlist2")
    end if

    ! Solve combined map, summed over all pixels 
    call wall_time(t1)
    call update_status(status, "shared1")
    if (sA_map%init) then
       do i = 0, self%numprocs_shared-1
          !call update_status(status, "tod_share_loop")
          start_chunk = mod(self%myid_shared+i,self%numprocs_shared)*chunk_size
          end_chunk   = min(start_chunk+chunk_size-1,npix-1)
          do while (start_chunk < npix)
             if (self%pix2ind(start_chunk) /= -1) exit
             start_chunk = start_chunk+1
          end do
          do while (end_chunk >= start_chunk)
             if (self%pix2ind(end_chunk) /= -1) exit
             end_chunk = end_chunk-1
          end do
          if (start_chunk < npix)       start_chunk = self%pix2ind(start_chunk)
          if (end_chunk >= start_chunk) end_chunk   = self%pix2ind(end_chunk)
          
          call mpi_win_fence(0, sA_map%win, ierr)
          call mpi_win_fence(0, sb_map%win, ierr)
          do j = start_chunk, end_chunk
             sA_map%a(:,self%ind2pix(j)+1) = sA_map%a(:,self%ind2pix(j)+1) + &
                  & A_map(:,j)
             sb_map%a(:,:,self%ind2pix(j)+1) = sb_map%a(:,:,self%ind2pix(j)+1) + &
                  & b_map(:,:,j)
          end do
          if (do_oper(samp_mono)) then
             call mpi_win_fence(0, sb_mono%win, ierr)
             do j = start_chunk, end_chunk
                sb_mono%a(:,self%ind2pix(j)+1,:) = sb_mono%a(:,self%ind2pix(j)+1,:) + &
                     & b_mono(:,j,:)
             end do
          end if
       end do
       call mpi_win_fence(0, sA_map%win, ierr)
       call mpi_win_fence(0, sb_map%win, ierr)
       if (do_oper(samp_mono)) call mpi_win_fence(0, sb_mono%win, ierr)

       call update_status(status, "shared2")

       call update_status(status, "finalize1")
       Sfilename = trim(prefix) // 'Smap'// trim(postfix)
!!$ HKE: commented out this
!!$       if (do_oper(samp_mono)) then
!!$          if (do_oper(prep_relbp)) then
!!$             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps, sb_mono=sb_mono, sys_mono=sys_mono, chisq_S=chisq_S, Sfile=Sfilename, mask=sprocmask2%a)
!!$          else
!!$             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps, sb_mono=sb_mono, sys_mono=sys_mono)
!!$          end if
!!$       else
!!$          !condmap => comm_map(self%info)
!!$          if (do_oper(prep_relbp)) then
!!$             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps, chisq_S=chisq_S, Sfile=Sfilename, mask=sprocmask2%a)
!!$          else
!!$             call finalize_binned_map(self, handle, sA_map, sb_map, rms_out, outmaps=outmaps)
!!$          end if
!!$          !call condmap%writeFITS("cond.fits")
!!$          !call condmap%dealloc()
!!$       end if

       if (do_oper(samp_bp)) then
          call wall_time(t1)
          call sample_bp(self, iter, delta, map_sky, handle, chisq_S)
          call wall_time(t2); t_tot(17) = t_tot(17) + t2-t1
       end if

       call update_status(status, "finalize2")
       map_out%map = outmaps(1)%p%map
       
       ! Sample monopole coefficients
!!$ HKE: Commented out this; no support for monopole sampling anymore
!!$       if (do_oper(samp_mono)) then
!!$          call sample_mono(self, handle, sys_mono, outmaps(2)%p, rms_out, &
!!$               & self%procmask)
!!$       end if

       ! Update bandpass parameters
       self%bp_delta = delta(:,:,1)

       ! Output maps to disk
       call map_out%writeFITS(trim(prefix)//'map'//trim(postfix))
       call rms_out%writeFITS(trim(prefix)//'rms'//trim(postfix))

       if (self%output_n_maps > 1) call outmaps(2)%p%writeFITS(trim(prefix)//'res'//trim(postfix))
       if (self%output_n_maps > 2) call outmaps(3)%p%writeFITS(trim(prefix)//'ncorr'//trim(postfix))
       if (self%output_n_maps > 3) call outmaps(4)%p%writeFITS(trim(prefix)//'bpcorr'//trim(postfix))
       if (self%output_n_maps > 4) call outmaps(5)%p%writeFITS(trim(prefix)//'mono'//trim(postfix))
       if (self%output_n_maps > 5) call outmaps(6)%p%writeFITS(trim(prefix)//'orb'//trim(postfix))
       if (self%output_n_maps > 6) call outmaps(7)%p%writeFITS(trim(prefix)//'sl'//trim(postfix))
       call update_status(status, "finalize3")

    end if



   !  if (self%first_call) then ! Harald
     if (do_oper(sel_data)) then 
       call mpi_reduce(ntot,    i, 1, MPI_INTEGER, MPI_SUM, &
            & self%numprocs/2, self%info%comm, ierr)
       ntot = i
       call mpi_reduce(naccept, i, 1, MPI_INTEGER, MPI_SUM, &
            & self%numprocs/2, self%info%comm, ierr)
       naccept = i
    end if
    call wall_time(t2); t_tot(10) = t_tot(10) + t2-t1

    call wall_time(t6)
    if (self%myid == self%numprocs/2) then
       write(*,*) '  Time dist sky   = ', t_tot(9)
       write(*,*) '  Time sl precomp = ', t_tot(13)
       write(*,*) '  Time decompress = ', t_tot(11)
       write(*,*) '  Time alloc      = ', t_tot(18)
       write(*,*) '  Time project    = ', t_tot(1)
       write(*,*) '  Time orbital    = ', t_tot(2)
       write(*,*) '  Time sl interp  = ', t_tot(12)
       write(*,*) '  Time ncorr      = ', t_tot(3)
       write(*,*) '  Time gain       = ', t_tot(4)
       write(*,*) '  Time absgain    = ', t_tot(14)
       write(*,*) '  Time sel data   = ', t_tot(15)
       write(*,*) '  Time clean      = ', t_tot(5)
       write(*,*) '  Time noise      = ', t_tot(6)
       write(*,*) '  Time samp abs   = ', t_tot(16)
       write(*,*) '  Time samp bp    = ', t_tot(17)
       write(*,*) '  Time chisq      = ', t_tot(7)
       write(*,*) '  Time bin        = ', t_tot(8)
       write(*,*) '  Time scanlist   = ', t_tot(20)
       write(*,*) '  Time final      = ', t_tot(10)
      !  if (self%first_call) then ! Harald
       if (do_oper(sel_data)) then
          write(*,*) '  Time total      = ', t6-t5, &
               & ', accept rate = ', real(naccept,sp) / ntot 
          write(*,*) 'naccept', naccept
          write(*,*) 'real(naccept,sp)', real(naccept,sp)
          write(*,*) 'ntot', ntot
       else
          write(*,*) '  Time total      = ', t6-t5, sum(t_tot(1:18))
       end if      
    end if


    ! Clean up temporary arrays
    deallocate(A_abscal, b_abscal, chisq_S)
    if (allocated(A_map)) deallocate(A_map, b_map)
    if (sA_map%init)  call dealloc_shared_2d_dp(sA_map)
    if (sb_map%init)  call dealloc_shared_3d_dp(sb_map)
    if (sb_mono%init) call dealloc_shared_3d_dp(sb_mono)
    if (allocated(b_mono)) deallocate(b_mono)
    if (allocated(sys_mono)) deallocate(sys_mono)
    if (allocated(slist)) deallocate(slist)

    if (allocated(outmaps)) then
       do i = 1, nout
          call outmaps(i)%p%dealloc
       end do
       deallocate(outmaps)
    end if

    if (sprocmask%init)  call dealloc_shared_1d_int(sprocmask)
    if (sprocmask2%init) call dealloc_shared_1d_int(sprocmask2)
    deallocate(map_sky)

    if (correct_sl) then
       do i = 1, self%ndet
          call self%slconv(i)%p%dealloc()
       end do
    end if

    call int2string(iter, ctext)
    call update_status(status, "tod_end"//ctext)

    ! Parameter to check if this is first time routine has been 
    self%first_call = .false.

  end subroutine process_SPIDER_tod


  integer(i4b) function count_lines(filename,datadir)
   implicit none
   character(len=*), intent(in) :: filename, datadir

   character(len=500)           :: detlist_file
   integer(i4b)                 :: unit,io_error,counter
   logical                      :: counting
   character(len=8)             :: line

   unit = 20
   detlist_file = trim(datadir)//trim(filename)

   open(unit,file=trim(detlist_file),status='old',action='read',iostat=io_error)
   if (io_error == 0) then
      ! Do nothing
   else
      stop 'Could not open file.'
   end if

   counting = .true.
   counter = 0
   do while(counting)
      read(unit,'(a)',end=1) line
      if ((line(1:1) == '#') .or. (line(1:1) == '')) then
         cycle
      else
         counter = counter + 1
      end if
   end do
1  close(unit)

   count_lines = counter
  end function count_lines

  subroutine get_det_labels(filename,datadir,labels)
   implicit none
   character(len=*), intent(in)    :: filename, datadir
   character(len=*), intent(inout) :: labels(:)

   character(len=500)           :: detlist_file
   integer(i4b)                 :: unit,io_error,counter, ndet, i
   character(len=8)             :: line

   ndet = size(labels)
   unit = 20
   detlist_file = trim(datadir)//trim(filename)

   open(unit,file=trim(detlist_file),status='old',action='read',iostat=io_error)
   if (io_error == 0) then
      ! Do nothing
   else
      stop 'Could not open file.'
   end if

   do i=1, ndet
      read(unit,'(a)') line
      if ((line(1:1) == '#') .or. (line(1:1) == '')) then
         cycle
      else
         labels(i) = line
      end if
   end do
end subroutine get_det_labels


subroutine write2file(filename, iter, param)
   implicit none

   character(len=*), intent(in)         :: filename
   real(dp), intent(in)                 :: param
   integer(i4b), intent(in)             :: iter

   integer(i4b)                           :: unit, io_error
   logical                                :: existing

   unit = 22

   inquire(file=trim(filename),exist=existing)
   if (existing) then
      open(unit,file=trim(filename),status='old',position='append',action='write',iostat=io_error)
   else
      open(unit,file=trim(filename),status='replace',action='write',iostat=io_error)
   end if

   write(unit,*), iter, param

   close(unit)
 end subroutine write2file


  subroutine read_tod_inst_SPIDER(self, file)
    ! 
    ! Reads SPIDER-specific common fields from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_SPIDER_tod)
    !           SPIDER-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle; only root includes this
    !
    ! Returns
    ! ----------
    ! None, but updates self
    !
    implicit none
    class(comm_SPIDER_tod),              intent(inout)          :: self
    type(hdf_file),                      intent(in),   optional :: file
  end subroutine read_tod_inst_SPIDER
  
  subroutine read_scan_inst_SPIDER(self, file, slabel, detlabels, scan)
    ! 
    ! Reads SPIDER-specific scan information from TOD fileset
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_SPIDER_tod)
    !           SPIDER-specific TOD object
    ! file:     derived type (hdf_file)
    !           Already open HDF file handle
    ! slabel:   string
    !           Scan label, e.g., "000001/"
    ! detlabels: string (array)
    !           Array of detector labels, e.g., ["27M", "27S"]
    ! scan:     derived class (comm_scan)
    !           Scan object
    !
    ! Returns
    ! ----------
    ! None, but updates scan object
    !
    implicit none
    class(comm_SPIDER_tod),              intent(in)    :: self
    type(hdf_file),                      intent(in)    :: file
    character(len=*),                    intent(in)    :: slabel
    character(len=*), dimension(:),      intent(in)    :: detlabels
    class(comm_scan),                    intent(inout) :: scan
  end subroutine read_scan_inst_SPIDER

  subroutine initHDF_SPIDER(self, chainfile, path)
    ! 
    ! Initializes SPIDER-specific TOD parameters from existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_SPIDER_tod)
    !           SPIDER-specific TOD object
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
    class(comm_SPIDER_tod),              intent(inout)  :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine initHDF_SPIDER
  
  subroutine dumpToHDF_SPIDER(self, chainfile, path)
    ! 
    ! Writes SPIDER-specific TOD parameters to existing chain file
    ! 
    ! Arguments:
    ! ----------
    ! self:     derived class (comm_SPIDER_tod)
    !           SPIDER-specific TOD object
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
    class(comm_SPIDER_tod),              intent(in)     :: self
    type(hdf_file),                      intent(in)     :: chainfile
    character(len=*),                    intent(in)     :: path
  end subroutine dumpToHDF_SPIDER



end module comm_tod_SPIDER_mod
