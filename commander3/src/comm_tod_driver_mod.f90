module comm_tod_driver_mod
  use comm_tod_mod
  use comm_tod_gain_mod
  use comm_tod_noise_mod
  use comm_tod_pointing_mod
  use comm_tod_bandpass_mod
  use comm_tod_orbdipole_mod
  use comm_tod_simulations_mod
  use comm_tod_mapmaking_mod
  use comm_zodi_mod
  use comm_shared_arr_mod
  use omp_lib
  implicit none


  ! Class for uncompressed data for a given scan
  type :: comm_scandata
     integer(i4b) :: ntod, ndet, nhorn, ndelta
     real(sp),     allocatable, dimension(:,:)     :: tod        ! Raw data
     real(sp),     allocatable, dimension(:,:)     :: n_corr     ! Correlated noise in V
     real(sp),     allocatable, dimension(:,:)     :: s_sl       ! Sidelobe correction
     real(sp),     allocatable, dimension(:,:)     :: s_sky      ! Stationary sky signal
     real(sp),     allocatable, dimension(:,:,:)   :: s_sky_prop ! Stationary sky signal proposal for bandpass sampling
     real(sp),     allocatable, dimension(:,:)     :: s_orb      ! Orbital dipole
     real(sp),     allocatable, dimension(:,:)     :: s_mono     ! Detector monopole correction 
     real(sp),     allocatable, dimension(:,:)     :: s_bp       ! Bandpass correction
     real(sp),     allocatable, dimension(:,:,:)   :: s_bp_prop  ! Bandpass correction proposal     
     real(sp),     allocatable, dimension(:,:)     :: s_zodi     ! Zodiacal light
     real(sp),     allocatable, dimension(:,:)     :: s_inst     ! Instrument-specific correction template
     real(sp),     allocatable, dimension(:,:)     :: s_tot      ! Total signal
     real(sp),     allocatable, dimension(:,:)     :: mask       ! TOD mask (flags + main processing mask)
     real(sp),     allocatable, dimension(:,:)     :: mask2      ! Small TOD mask, for bandpass sampling
     integer(i4b), allocatable, dimension(:,:,:)   :: pix        ! Discretized pointing 
     integer(i4b), allocatable, dimension(:,:,:)   :: psi        ! Discretized polarization angle
     integer(i4b), allocatable, dimension(:,:)     :: flag       ! Quality flags

     real(sp),     allocatable, dimension(:,:)     :: s_totA     ! Total signal, horn A (differential only)
     real(sp),     allocatable, dimension(:,:)     :: s_totB     ! Total signal, horn B (differential only)
   contains
     procedure  :: init_singlehorn   => init_scan_data_singlehorn
     procedure  :: init_differential => init_scan_data_differential
     procedure  :: dealloc           => dealloc_scan_data
  end type comm_scandata


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Scan data routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_scan_data_singlehorn(self, tod, scan, map_sky, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop)
    implicit none
    class(comm_scandata),                      intent(inout)          :: self    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_

    if (tod%nhorn /= 1) then
       write(*,*) 'Error: init_scan_data_singlehorn only applicable for 1-horn experiments'
       stop
    end if
    
    !if (.true. .or. tod%myid == 78) write(*,*) 'c', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    init_s_bp_ = .false.; if (present(init_s_bp)) init_s_bp_ = init_s_bp
    init_s_sky_prop_ = .false.; if (present(init_s_sky_prop)) init_s_sky_prop_ = init_s_sky_prop
    init_s_bp_prop_ = .false.
    if (present(init_s_bp_prop)) then
       init_s_bp_prop_  = init_s_bp_prop
       init_s_sky_prop_ = init_s_bp_prop
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c1', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    self%ntod   = tod%scans(scan)%ntod
    self%ndet   = tod%ndet
    self%nhorn  = tod%nhorn
    self%ndelta = 0; if (present(init_s_sky_prop)) self%ndelta = size(map_sky,4)

    ! Allocate data structures
    allocate(self%tod(self%ntod, self%ndet))
    allocate(self%n_corr(self%ntod, self%ndet))
    allocate(self%s_sl(self%ntod, self%ndet))
    allocate(self%s_sky(self%ntod, self%ndet))
    allocate(self%s_bp(self%ntod, self%ndet))
    allocate(self%s_orb(self%ntod, self%ndet))
    allocate(self%s_tot(self%ntod, self%ndet))
    allocate(self%mask(self%ntod, self%ndet))
    allocate(self%pix(self%ntod, self%ndet, self%nhorn))
    allocate(self%psi(self%ntod, self%ndet, self%nhorn))
    allocate(self%flag(self%ntod, self%ndet))
    if (init_s_sky_prop_)    allocate(self%s_sky_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_bp_prop_)     allocate(self%s_bp_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_sky_prop_)    allocate(self%mask2(self%ntod, self%ndet))
    if (tod%sample_mono)     allocate(self%s_mono(self%ntod, self%ndet))
    if (tod%subtract_zodi)   allocate(self%s_zodi(self%ntod, self%ndet))
    if (tod%apply_inst_corr) allocate(self%s_inst(self%ntod, self%ndet))

    !if (.true. .or. tod%myid == 78) write(*,*) 'c2', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Decompress pointing, psi and flags for current scan
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       call tod%decompress_pointing_and_flags(scan, j, self%pix(:,j,:), &
            & self%psi(:,j,:), self%flag(:,j))
    end do
    !if (tod%myid == 78) write(*,*) 'c3', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    if (tod%symm_flags) call tod%symmetrize_flags(self%flag)
    !if (.true. .or. tod%myid == 78) write(*,*) 'c4', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    
    ! Prepare TOD
    if (tod%ndiode == 1) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          if (tod%compressed_tod) then
             call tod%decompress_tod(scan, j, self%tod(:,j))
          else
             self%tod(:,j) = tod%scans(scan)%d(j)%tod
          end if
       end do
    else
       call tod%diode2tod_inst(scan, self%tod)
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c5', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Construct sky signal template
    if (init_s_bp_) then
       call project_sky(tod, map_sky(:,:,:,1), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
            & procmask, scan, self%s_sky, self%mask, s_bp=self%s_bp)
    else
       call project_sky(tod, map_sky(:,:,:,1), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
            & procmask, scan, self%s_sky, self%mask)
    end if
    !if (tod%myid == 78) write(*,*) 'c6', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    if (init_s_bp_prop_) then
       do j = 2, size(map_sky,4)
          !if (.true. .or. tod%myid == 78) write(*,*) 'c61', j, tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires, size(map_sky,4)
          call project_sky(tod, map_sky(:,:,:,j), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
               & procmask2, scan, self%s_sky_prop(:,:,j), self%mask2, s_bp=self%s_bp_prop(:,:,j))
       end do
    else if (init_s_sky_prop_) then
       do j = 2, size(map_sky,4)
          !if (.true. .or. tod%myid == 78) write(*,*) 'c62', j, tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
          call project_sky(tod, map_sky(:,:,:,j), self%pix(:,:,1), self%psi(:,:,1), self%flag, &
               & procmask2, scan, self%s_sky_prop(:,:,j), self%mask2)
       end do
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c71', tod%myid, tod%correct_sl
    !if (.true. .or. tod%myid == 78) write(*,*) 'c72', tod%myid, tod%ndet
    !if (.true. .or. tod%myid == 78) write(*,*) 'c73', tod%myid, tod%slconv(1)%p%psires

    ! Perform sanity tests
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(self%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%N_psd%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    !if (.true. .or. tod%myid == 78) write(*,*) 'c8', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    
    ! Construct orbital dipole template
    call tod%construct_dipole_template(scan, self%pix(:,:,1), self%psi(:,:,1), .true., self%s_orb)
    !if (.true. .or. tod%myid == 78) write(*,*) 'c9', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Construct zodical light template
    if (tod%subtract_zodi) then
       call compute_zodi_template(tod%nside, self%pix(:,:,1), tod%scans(scan)%satpos, tod%nu_c, self%s_zodi)
    end if
    !if (.true. .or. tod%myid == 78) write(*,*) 'c10', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires

    ! Construct sidelobe template
    !if (.true. .or. tod%myid == 78) write(*,*) 'd', tod%myid, tod%correct_sl, tod%ndet, tod%slconv(1)%p%psires
    if (tod%correct_sl) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !if (.true. .or. tod%myid == 78) write(*,*) 'e', tod%myid, j, tod%slconv(j)%p%psires, tod%slconv(j)%p%psisteps
          call tod%construct_sl_template(tod%slconv(j)%p, &
               & self%pix(:,j,1), self%psi(:,j,1), self%s_sl(:,j), tod%polang(j))
          self%s_sl(:,j) = 2.d0 * self%s_sl(:,j) ! Scaling by a factor of 2, by comparison with LevelS. Should be understood
       end do
    else
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          self%s_sl(:,j) = 0.
       end do
    end if

    ! Construct monopole correction template
    if (tod%sample_mono) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !self%s_mono(:,j) = tod%mono(j)
          self%s_mono(:,j) = 0.d0 ! Disabled for now
       end do
    end if

    ! Generate and apply instrument-specific correction template
    if (tod%apply_inst_corr) then
       call tod%construct_corrtemp_inst(scan, self%pix(:,:,1), self%psi(:,:,1), self%s_inst)
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          self%tod(:,j) = self%tod(:,j) - self%s_inst(:,j)
       end do
    end if

    ! Construct total sky signal
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       self%s_tot(:,j) = self%s_sky(:,j) + self%s_sl(:,j) + self%s_orb(:,j)
       if (tod%sample_mono) self%s_tot(:,j) = self%s_tot(:,j) + self%s_mono(:,j)
    end do

  end subroutine init_scan_data_singlehorn


  subroutine init_scan_data_differential(self, tod, scan, map_sky, procmask, procmask2, &
       & init_s_bp, init_s_bp_prop, init_s_sky_prop)
    implicit none
    class(comm_scandata),                      intent(inout)          :: self    
    class(comm_tod),                           intent(inout)          :: tod
    integer(i4b),                              intent(in)             :: scan
    real(sp),          dimension(1:,1:,0:,1:), intent(in)             :: map_sky
    real(sp),          dimension(0:),          intent(in)             :: procmask
    real(sp),          dimension(0:),          intent(in)             :: procmask2
    logical(lgt),                              intent(in),   optional :: init_s_bp
    logical(lgt),                              intent(in),   optional :: init_s_bp_prop
    logical(lgt),                              intent(in),   optional :: init_s_sky_prop

    integer(i4b) :: j, k, ndelta
    logical(lgt) :: init_s_bp_, init_s_bp_prop_, init_s_sky_prop_
    real(sp),     allocatable, dimension(:,:)     :: s_bufA, s_bufB, s_buf2A, s_buf2B      ! Buffer

    if (tod%nhorn /= 2) then
       write(*,*) 'Error: init_scan_data_differential only applicable for 2-horn experiments'
       stop
    end if

    init_s_bp_ = .false.; if (present(init_s_bp)) init_s_bp_ = init_s_bp
    init_s_sky_prop_ = .false.; if (present(init_s_sky_prop)) init_s_sky_prop_ = init_s_sky_prop
    init_s_bp_prop_ = .false.
    if (present(init_s_bp_prop)) then
       init_s_bp_prop_ = init_s_bp_prop
       init_s_sky_prop_ = init_s_sky_prop
    end if

    self%ntod   = tod%scans(scan)%ntod
    self%ndet   = tod%ndet
    self%nhorn  = tod%nhorn
    self%ndelta = 0; if (present(init_s_sky_prop)) self%ndelta = size(map_sky,4)

    ! Allocate data structures
    allocate(self%tod(self%ntod, self%ndet))
    allocate(self%n_corr(self%ntod, self%ndet))
    allocate(self%s_sl(self%ntod, self%ndet))
    allocate(self%s_sky(self%ntod, self%ndet))
    allocate(self%s_bp(self%ntod, self%ndet))
    allocate(self%s_orb(self%ntod, self%ndet))
    allocate(self%s_tot(self%ntod, self%ndet))
    allocate(self%s_totA(self%ntod, self%ndet))
    allocate(self%s_totB(self%ntod, self%ndet))
    allocate(self%mask(self%ntod, self%ndet))
    allocate(self%pix(self%ntod, 1, self%nhorn))
    allocate(self%psi(self%ntod, 1, self%nhorn))
    allocate(self%flag(self%ntod, 1))
    if (init_s_sky_prop_)   allocate(self%s_sky_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_bp_prop_)    allocate(self%s_bp_prop(self%ntod, self%ndet, 2:self%ndelta))
    if (init_s_sky_prop_)   allocate(self%mask2(self%ntod, self%ndet))
    if (tod%sample_mono)    allocate(self%s_mono(self%ntod, self%ndet))
    if (tod%subtract_zodi)  allocate(self%s_zodi(self%ntod, self%ndet))
    self%s_tot  = 0.
    self%s_totA = 0.
    self%s_totB = 0.

    allocate(s_bufA(self%ntod, self%ndet))
    allocate(s_bufB(self%ntod, self%ndet))
    allocate(s_buf2A(self%ntod, self%ndet))
    allocate(s_buf2B(self%ntod, self%ndet))

    ! Decompress pointing, psi and flags for current scan
    ! Only called for one detector, det=1, since the pointing and polarization
    ! angles are the same for all detectors
    call tod%decompress_pointing_and_flags(scan, 1, self%pix(:,1,:), &
            & self%psi(:,1,:), self%flag(:,1))
    
    ! Prepare TOD
    if (tod%ndiode == 1) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          if (tod%compressed_tod) then
             call tod%decompress_tod(scan, j, self%tod(:,j))
          else
             self%tod(:,j) = tod%scans(scan)%d(j)%tod
          end if
       end do
    else
       call tod%diode2tod_inst(scan, self%tod)
    end if

    ! Construct sky signal template
    if (init_s_bp_) then
       call project_sky_differential(tod, map_sky(:,:,:,1), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, self%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
    else
       call project_sky_differential(tod, map_sky(:,:,:,1), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
            & procmask, scan, s_bufA, s_bufB, self%mask)
    end if
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       self%s_sky(:,j)  = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
       self%s_totA(:,j) = self%s_totA(:,j) + s_bufA(:,j)
       self%s_totB(:,j) = self%s_totB(:,j) + s_bufB(:,j)
       self%s_tot(:,j)  = self%s_tot(:,j)  + self%s_sky(:,j)
       if (init_s_bp_) self%s_bp(:,j)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
    end do

    ! Set up (optional) bandpass sampling quantities (s_sky_prop, mask2 and bp_prop)
    if (init_s_bp_prop_) then
       do k = 2, size(map_sky,4)
          call project_sky_differential(tod, map_sky(:,:,:,k), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
               & procmask, scan, s_bufA, s_bufB, self%mask, s_bpA=s_buf2A, s_bpB=s_buf2B)
          do j = 1, self%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             self%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
             if (init_s_bp_) self%s_bp_prop(:,j,k)  = (1.+tod%x_im(j))*s_buf2A(:,j) - (1.-tod%x_im(j))*s_buf2B(:,j)
          end do
       end do
    else if (init_s_sky_prop_) then
       do k = 2, size(map_sky,4)
          call project_sky_differential(tod, map_sky(:,:,:,k), self%pix(:,1,:), self%psi(:,1,:), self%flag(:,1), &
               & procmask, scan, s_bufA, s_bufB, self%mask)
          do j = 1, self%ndet
             if (.not. tod%scans(scan)%d(j)%accept) cycle
             self%s_sky_prop(:,j,k) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          end do
       end do
    end if

    ! Perform sanity tests
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (all(self%mask(:,j) == 0)) tod%scans(scan)%d(j)%accept = .false.
       if (tod%scans(scan)%d(j)%N_psd%sigma0 <= 0.d0) tod%scans(scan)%d(j)%accept = .false.
    end do
    
    ! Construct orbital dipole template
    call tod%construct_dipole_template_diff(scan, self%pix(:,:,1), self%psi(:,:,1), .true., s_bufA, 1d3)
    call tod%construct_dipole_template_diff(scan, self%pix(:,:,2), self%psi(:,:,2), .true., s_bufB, 1d3)
    do j = 1, self%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       self%s_orb(:,j)  = (1.+tod%x_im(j))*s_bufA(:,j)  - (1.-tod%x_im(j))*s_bufB(:,j)
       self%s_totA(:,j) = self%s_totA(:,j) + s_bufA(:,j)
       self%s_totB(:,j) = self%s_totB(:,j) + s_bufB(:,j)
       self%s_tot(:,j)  = self%s_tot(:,j)  + self%s_orb(:,j)
    end do


    ! Construct zodical light template
    if (tod%subtract_zodi) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          call compute_zodi_template(tod%nside, self%pix(:,1:1,1), tod%scans(scan)%satpos, tod%nu_c(j:j), s_bufA)
          call compute_zodi_template(tod%nside, self%pix(:,1:1,2), tod%scans(scan)%satpos, tod%nu_c(j:j), s_bufB)
          self%s_zodi(:,j) = (1.+tod%x_im(j))*s_bufA(:,j) - (1.-tod%x_im(j))*s_bufB(:,j)
          self%s_tot(:,j)  = self%s_tot(:,j) + self%s_zodi(:,j)
          self%s_totA(:,j) = self%s_totA(:,j) + s_bufA(:,j)
          self%s_totB(:,j) = self%s_totB(:,j) + s_bufB(:,j)
       end do
    else
       self%s_zodi = 0.
    end if

    ! Construct sidelobe template
    if (tod%correct_sl) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          call tod%construct_sl_template(tod%slconv(1)%p, self%pix(:,1,1), self%psi(:,1,1), s_bufA(:,j), 0d0)
          call tod%construct_sl_template(tod%slconv(3)%p, self%pix(:,1,2), self%psi(:,1,2), s_bufB(:,j), 0d0)
          self%s_sl(:,j)  = 2.*((1d0+tod%x_im(j))*s_bufA(:,j) - (1d0-tod%x_im(j))*s_bufB(:,j))
          self%s_tot(:,j) = self%s_tot(:,j) + self%s_sl(:,j)
          self%s_totA(:,j) = self%s_totA(:,j) + 2.*s_bufA(:,j)
          self%s_totB(:,j) = self%s_totB(:,j) + 2.*s_bufB(:,j)
       end do
    else
       self%s_sl = 0.
    end if

    ! Construct monopole correction template
    if (tod%sample_mono) then
       do j = 1, self%ndet
          if (.not. tod%scans(scan)%d(j)%accept) cycle
          !self%s_mono(:,j) = tod%mono(j)
          self%s_mono(:,j) = 0.d0 ! Disabled for now
          self%s_tot(:,j)  = self%s_tot(:,j) + self%s_mono(:,j)
          self%s_totA(:,j) = self%s_totA(:,j) + self%s_mono(:,j)
          self%s_totB(:,j) = self%s_totB(:,j) + self%s_mono(:,j)
       end do
    else
       self%s_mono = 0.d0 
    end if

    ! Clean-up
    deallocate(s_bufA, s_bufB, s_buf2A, s_buf2B)

  end subroutine init_scan_data_differential


  subroutine dealloc_scan_data(self)
    implicit none
    class(comm_scandata), intent(inout) :: self    

    self%ntod = -1; self%ndet = -1; self%nhorn = -1

    ! Deallocate data structures
    deallocate(self%tod, self%n_corr, self%s_sl, self%s_sky)
    deallocate(self%s_orb, self%s_tot, self%mask)
    deallocate(self%pix, self%psi, self%flag)
    if (allocated(self%s_sky_prop))  deallocate(self%s_sky_prop)
    if (allocated(self%s_bp_prop))   deallocate(self%s_bp_prop)
    if (allocated(self%s_bp))        deallocate(self%s_bp)
    if (allocated(self%s_mono))      deallocate(self%s_mono)
    if (allocated(self%mask2))       deallocate(self%mask2)
    if (allocated(self%s_zodi))      deallocate(self%s_zodi)
    if (allocated(self%s_totA))      deallocate(self%s_totA)
    if (allocated(self%s_totB))      deallocate(self%s_totB)
    if (allocated(self%s_inst))      deallocate(self%s_inst)

  end subroutine dealloc_scan_data


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Sampling drivers etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sample_calibration(tod, mode, handle, map_sky, procmask, procmask2)
    !   Sample calibration modes
    !   Supported modes = {abscal, relcal, deltaG, imbal}
    !
    !   Subroutine that implements gain and horn imbalance sampling, following
    !   the formalism of Gjerlow et al. 2021
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   mode:     character array
    !             specifies sampling mode. currently supports absolute calibration,
    !             relative calibration, time-variable calibration, and horn
    !             imbalance.
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    !             so that the same sequence can be resumed later on from that same point
    !   map_sky:
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    character(len=*),                             intent(in)    :: mode
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j, ext(2), ierr
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:)   :: A, b
    real(sp), allocatable, dimension(:,:) :: s_invN, mask_lowres, s_buf
    real(dp), allocatable, dimension(:,:) :: dipole_mod
    type(comm_scandata) :: sd

    if (tod%myid == 0) write(*,*) '   --> Sampling calibration, mode = ', trim(mode)

    if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
       allocate(A(tod%ndet), b(tod%ndet))
       A = 0.d0; b = 0.d0
    else if (trim(mode) == 'deltaG') then
       allocate(dipole_mod(tod%nscan_tot, tod%ndet))
       dipole_mod = 0.d0
    else
       write(*,*) 'Unsupported sampling mode!'
       stop
    end if

    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (tod%nhorn == 1) then
          call sd%init_singlehorn(tod, i, map_sky, procmask, procmask2)
       else
          call sd%init_differential(tod, i, map_sky, procmask, procmask2)
       end if

       ! Set up filtered calibration signal, conditional contribution and mask
       call tod%downsample_tod(sd%s_orb(:,1), ext)
       allocate(s_invN(ext(1):ext(2), tod%ndet))      ! s * invN
       allocate(s_buf(sd%ntod, sd%ndet))
       allocate(mask_lowres(ext(1):ext(2), tod%ndet))
       do j = 1, tod%ndet
          if (.not. tod%scans(i)%d(j)%accept) cycle
          call tod%downsample_tod(sd%mask(:,j), ext, mask_lowres(:,j), threshold=0.9)
          if (trim(mode) == 'abscal' .and. tod%orb_abscal) then
             ! Calibrator = orbital dipole only
             call tod%downsample_tod(sd%s_orb(:,j), ext, s_invN(:,j))
          else if (trim(mode) == 'imbal') then
             ! Calibrator = common mode signal
             s_buf(:,j) = sd%s_totA(:,j) + sd%s_totB(:,j)
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., &
               & real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invN(:,j))
          else
             ! Calibrator = total signal
             s_buf(:,j) = sd%s_tot(:,j)
             call fill_all_masked(s_buf(:,j), sd%mask(:,j), sd%ntod, .false., real(tod%scans(i)%d(j)%N_psd%sigma0, sp), handle, tod%scans(i)%chunk_num)
             call tod%downsample_tod(s_buf(:,j), ext, s_invN(:,j))
          end if
       end do
       !do j = 1, tod%ndet
       !      if (tod%scanid(i) == 30) then
       !        write(*,*) 'abscaltest1', j, 'sum(s_invN(:,j))', sum(s_invN(:,j))
       !      end if
       !end do
       call multiply_inv_N(tod, i, s_invN, sampfreq=tod%samprate_lowres, pow=0.5d0)
       !do j = 1, tod%ndet
       !      if (tod%scanid(i) == 30) then
       !        write(*,*) 'abscaltest2', j, 'sum(s_invN(:,j))', sum(s_invN(:,j))
       !      end if
       !end do

       if (trim(mode) == 'abscal' .or. trim(mode) == 'relcal' .or. trim(mode) == 'imbal') then
          ! Constant gain terms; accumulate contribution from this scan
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             !if (trim(mode) == 'abscal') then
             !  write(*,*) 'test', tod%scanid(i), j, tod%gain0(j), tod%scans(i)%d(j)%dgain, &
             !             & sum(abs(1.d0*sd%s_tot(:, j))), sum(1.d0*abs(sd%s_orb(:,j)))
             !end if
             !if (tod%scanid(i)==30) write(*,*) 'j, tod%gain0(j), tod%scans(i)%d(j)%dgain,'//&
             !                &' sum(abs(1.d0*sd%s_tot(:, j))), sum(1.d0*abs(sd%s_orb(:,j)))'
             !if (tod%scanid(i)==30) write(*,*) j, tod%gain0(j), tod%scans(i)%d(j)%dgain, &
             !                       & sum(abs(1.d0*sd%s_tot(:, j))), sum(1.d0*abs(sd%s_orb(:,j)))
             if (trim(mode) == 'abscal' .and. tod%orb_abscal) then
                s_buf(:,j) = real(tod%gain0(0),sp) * (sd%s_tot(:,j) - sd%s_orb(:,j)) + &
                     & real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j) + &
                     & tod%scans(i)%d(j)%baseline
             else if (trim(mode) == 'abscal' .and. .not. tod%orb_abscal) then
                s_buf(:,j) = real(tod%gain0(j) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)! + &
                     !& tod%scans(i)%d(j)%baseline
             else if (trim(mode) == 'relcal') then
                s_buf(:,j) = real(tod%gain0(0) + tod%scans(i)%d(j)%dgain,sp) * sd%s_tot(:,j)! + &
                     !& tod%scans(i)%d(j)%baseline
             else if (trim(mode) == 'imbal') then
                s_buf(:,j) = tod%scans(i)%d(j)%gain * (sd%s_totA(:,j) - sd%s_totB(:,j))! + &
                     !& tod%scans(i)%d(j)%baseline
             end if
          end do
          if (tod%compressed_tod) then
            call accumulate_abscal(tod, i, sd%mask, s_buf, s_invN, s_invN, A, b, handle, &
              & out=.true., mask_lowres=mask_lowres, tod_arr=sd%tod)
          else
            call accumulate_abscal(tod, i, sd%mask, s_buf, s_invN, s_invN, A, b, handle, &
              & out=.true., mask_lowres=mask_lowres)
          end if
       else
          ! Time-variable gain terms
          if (tod%compressed_tod) then
            call calculate_gain_mean_std_per_scan(tod, i, s_invN, sd%mask, s_invN, sd%s_tot, &
              & handle, mask_lowres=mask_lowres, tod_arr=sd%tod)
          else
            call calculate_gain_mean_std_per_scan(tod, i, s_invN, sd%mask, s_invN, sd%s_tot, &
              & handle, mask_lowres=mask_lowres)
          end if
          do j = 1, tod%ndet
             if (.not. tod%scans(i)%d(j)%accept) cycle
             dipole_mod(tod%scanid(i),j) = masked_variance(sd%s_sky(:,j), sd%mask(:,j))
          end do
       end if

       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call sd%dealloc
       deallocate(s_invN, s_buf, mask_lowres)
    end do

    ! Perform sampling operations
    if (trim(mode) == 'abscal') then
       call sample_abscal_from_orbital(tod, handle, A, b)
    else if (trim(mode) == 'relcal') then
       call sample_relcal(tod, handle, A, b)
    else if (trim(mode) == 'deltaG') then
       call mpi_allreduce(mpi_in_place, dipole_mod, size(dipole_mod), MPI_DOUBLE_PRECISION, MPI_SUM, tod%info%comm, ierr)
       call sample_smooth_gain(tod, handle, dipole_mod)
    else if (trim(mode) == 'imbal') then
       call sample_imbal_cal(tod, handle, A, b)
    end if

    ! Clean up
    if (allocated(A))          deallocate(A)
    if (allocated(b))          deallocate(b)
    if (allocated(dipole_mod)) deallocate(dipole_mod)

  end subroutine sample_calibration


  ! Sample baseline
  subroutine sample_baseline(tod, handle, map_sky, procmask, procmask2)
    implicit none
    class(comm_tod),                              intent(inout) :: tod
    type(planck_rng),                             intent(inout) :: handle
    real(sp),            dimension(0:,1:,1:,1:),  intent(in)    :: map_sky
    real(sp),            dimension(0:),           intent(in)    :: procmask, procmask2

    integer(i4b) :: i, j
    real(dp)     :: t1, t2
    type(comm_scandata) :: sd

    if (tod%myid == 0) write(*,*) '   --> Sampling baseline'

    do i = 1, tod%nscan
       if (.not. any(tod%scans(i)%d%accept)) cycle
       call wall_time(t1)

       ! Prepare data
       if (tod%nhorn == 1) then
          call sd%init_singlehorn(tod, i, map_sky, procmask, procmask2)
       else
          call sd%init_differential(tod, i, map_sky, procmask, procmask2)
       end if

       do j = 1, tod%ndet
          tod%scans(i)%d(j)%baseline =sum((sd%tod(:,j) - tod%scans(i)%d(j)%gain*sd%s_tot(:,j)) &
            & *sd%mask(:,j))/sum(sd%mask(:,j))
          if (trim(tod%operation) == 'sample') then
            tod%scans(i)%d(j)%baseline = tod%scans(i)%d(j)%baseline &
             &  + rand_gauss(handle)/sqrt(sum(sd%mask(:,j)*tod%scans(i)%d(j)%N_psd%sigma0**2))
          end if
       end do

       ! Clean up
       call wall_time(t2)
       tod%scans(i)%proctime   = tod%scans(i)%proctime   + t2-t1
       tod%scans(i)%n_proctime = tod%scans(i)%n_proctime + 1
       call sd%dealloc
    end do
    !do j = 1, tod%ndet
    !  if (tod%myid == 0) then
    !    call sd%init_differential(tod, 1, map_sky, procmask, procmask2)
    !    write(*,*) 'Detector',j
    !    write(*,*) tod%scans(1)%d(j)%baseline
    !    write(*,*) sum(sd%tod(:,j))/size(sd%tod(:,j))
    !    write(*,*) sum(sd%tod(:,j) - tod%scans(1)%d(j)%baseline)/size(sd%tod(:,j))
    !    call sd%dealloc
    !  end if
    !end do

  end subroutine sample_baseline

  subroutine remove_bad_data(tod, scan, flag)
    !   Perform data selection on TOD object
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information. Bad data are removed by 
    !             setting scan%det%accept = .false.
    !   scan:     int (scalar)
    !             Local scan ID for the current core 
    !   flag:     int (ntod x ndet array)
    !             Array with data quality flags
    implicit none
    class(comm_tod),                   intent(inout) :: tod
    integer(i4b),    dimension(1:,1:), intent(in)    :: flag
    integer(i4b),                      intent(in)    :: scan

    integer(i4b) :: j, ntod, ndet
    
    ntod = size(flag,1)
    ndet = size(flag,2)
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (count(iand(flag(:,j),tod%flag0) .ne. 0) > 0.1*ntod) then    ! Discard scans with less than 10% good data
          tod%scans(scan)%d(j)%accept = .false.
       else if (abs(tod%scans(scan)%d(j)%chisq) > tod%chisq_threshold .or. &  ! Discard scans with high chisq or NaNs
            & isNaN(tod%scans(scan)%d(j)%chisq)) then
          write(*,fmt='(a,i8,i5,a,f12.1)') 'Reject scan, det = ', &
               & tod%scanid(scan), j, ', chisq = ', tod%scans(scan)%d(j)%chisq
          tod%scans(scan)%d(j)%accept = .false.
       end if
    end do
       !if (any(.not. tod%scans(scan)%d%accept)) tod%scans(scan)%d%accept = .false. ! Do we actually want this..?
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) tod%scans(scan)%d(tod%partner(j))%accept = .false.
    end do

  end subroutine remove_bad_data

  subroutine compute_chisq_abs_bp(tod, scan, sd, chisq)
    implicit none
    class(comm_tod),                       intent(inout) :: tod
    integer(i4b),                          intent(in)    :: scan
    type(comm_scandata),                   intent(in)    :: sd
    real(dp),            dimension(:,:),   intent(inout) :: chisq

    integer(i4b) :: j, k
    real(sp), allocatable, dimension(:,:) :: s_buf

    allocate(s_buf(sd%ntod,sd%ndet))
    do j = 1, tod%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       s_buf(:,j) =  sd%s_sl(:,j) + sd%s_orb(:,j)
       call tod%compute_chisq(scan, j, sd%mask2(:,j), sd%s_sky(:,j), &
            & s_buf(:,j), sd%n_corr(:,j), sd%tod(:,j), absbp=.true.)
       chisq(j,1) = chisq(j,1) + tod%scans(scan)%d(j)%chisq_prop
       do k = 2, size(sd%s_sky_prop)
          call tod%compute_chisq(scan, j, sd%mask2(:,j), sd%s_sky_prop(:,j,k), &
               & s_buf(:,j), sd%n_corr(:,j), sd%tod(:,j), absbp=.true.)
          chisq(j,k) = chisq(j,k) + tod%scans(scan)%d(j)%chisq_prop
       end do
    end do
    deallocate(s_buf)

  end subroutine compute_chisq_abs_bp

  subroutine compute_calibrated_data(tod, scan, sd, d_calib)
    !
    !  gets calibrated timestreams
    !
    !  Arguments:
    !  ----------
    !  tod: comm_tod object
    !
    !  scan: integer
    !     integer label for scan
    !  sd:  comm_scandata object
    !
    !  Returns:
    !  --------
    !  d_calib: real(sp) array
    !     nout x ndet x ntod array of calibrated timestreams
    !
    !  d_calib(1,:,:) - best estimate of calibrated data, with all known
    !    calibrations applied
    !  d_calib(2,:,:) - calibrated TOD with expected sky signal subtracted,
    !    i.e., residual
    !  d_calib(3,:,:) - correlated noise, mean subtracted, in temperature
    !    units
    !  d_calib(4,:,:) - bandpass difference contribution
    !  d_calib(5,:,:) - orbital dipole
    !  d_calib(6,:,:) - sidelobe
    !  d_calib(7,:,:) - zodiacal light emission
    !
    implicit none
    class(comm_tod),                       intent(in)   :: tod
    integer(i4b),                          intent(in)   :: scan
    type(comm_scandata),                   intent(in)   :: sd
    real(sp),            dimension(:,:,:), intent(out)  :: d_calib

    integer(i4b) :: j, nout
    real(dp)     :: inv_gain

    nout = size(d_calib,1)
    do j = 1, sd%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       inv_gain = 1.0 / tod%scans(scan)%d(j)%gain
       if (tod%compressed_tod) then
        d_calib(1,:,j) = (sd%tod(:,j) - tod%scans(scan)%d(j)%baseline- sd%n_corr(:,j)) &
          & * inv_gain - sd%s_tot(:,j) + sd%s_sky(:,j) - sd%s_bp(:,j)
       else
        d_calib(1,:,j) = (tod%scans(scan)%d(j)%tod - tod%scans(scan)%d(j)%baseline- sd%n_corr(:,j)) &
          & * inv_gain - sd%s_tot(:,j) + sd%s_sky(:,j) - sd%s_bp(:,j)
       end if
       if (nout > 1) d_calib(2,:,j) = d_calib(1,:,j) - sd%s_sky(:,j) + sd%s_bp(:,j)              ! residual
       if (nout > 2) d_calib(3,:,j) = (sd%n_corr(:,j) - sum(sd%n_corr(:,j)/sd%ntod)) * inv_gain  ! ncorr
       if (nout > 3) d_calib(4,:,j) = sd%s_bp(:,j)                                               ! bandpass
       if (nout > 4) d_calib(5,:,j) = sd%s_orb(:,j)                                              ! orbital dipole
       if (nout > 5) d_calib(6,:,j) = sd%s_sl(:,j)                                               ! sidelobes
       if (nout > 6) then
          if (allocated(sd%s_zodi)) then
             d_calib(7,:,j) = sd%s_zodi(:,j)                                                     ! zodi
          else
             d_calib(7,:,j) = 0.
          end if
       end if
    end do

  end subroutine compute_calibrated_data

  subroutine distribute_sky_maps(tod, map_in, scale, map_out)
    implicit none
    class(comm_tod),                       intent(in)     :: tod
    type(map_ptr), dimension(1:,1:),       intent(inout)  :: map_in       ! (ndet,ndelta)    
    real(sp),                              intent(in)     :: scale
    real(sp),      dimension(1:,1:,0:,1:), intent(out)    :: map_out

    integer(i4b) :: i, j, k, l, npix, nmaps
    real(dp),     allocatable, dimension(:,:) :: m_buf

    npix  = map_in(1,1)%p%info%npix
    nmaps = map_in(1,1)%p%info%nmaps
    allocate(m_buf(0:npix-1,nmaps))
    do j = 1, size(map_in,2)
       do i = 1, size(map_in,1)
          map_in(i,j)%p%map = scale * map_in(i,j)%p%map ! unit conversion
          call map_in(i,j)%p%bcast_fullsky_map(m_buf)
          do k = 1, tod%nobs
             map_out(:,k,i,j) = m_buf(tod%ind2pix(k),:)
          end do
       end do
       do k = 1, tod%nobs
          do l = 1, tod%nmaps
             map_out(l,k,0,j) = sum(map_out(l,k,1:tod%ndet,j))/tod%ndet
          end do
       end do
    end do
    deallocate(m_buf)

  end subroutine distribute_sky_maps

  ! ************************************************
  !
  !> @brief Commander3 native simulation module. It
  !! simulates correlated noise and then rewrites
  !! the original timestreams inside the files.
  !
  !> @author Maksym Brilenkov
  !
  !> @param[in]
  !> @param[out]
  !
  ! ************************************************
  subroutine simulate_tod(self, scan_id, s_tot, handle)
    implicit none
    class(comm_tod), intent(inout) :: self
    ! Parameter file variables
    !type(comm_params),                     intent(in)    :: cpar
    ! Other input/output variables
    real(sp), allocatable, dimension(:,:), intent(in)    :: s_tot   !< total sky signal
    integer(i4b),                          intent(in)    :: scan_id !< current PID
    type(planck_rng),                      intent(inout) :: handle
    ! Simulation variables
    real(sp), allocatable, dimension(:,:) :: tod_per_detector !< simulated tods per detector
    real(sp)                              :: gain   !< detector's gain value
    real(sp)                              :: sigma0
    real(sp) :: N_c
    real(sp) :: samprate
    real(sp) :: fft_norm
    integer(i4b)                          :: ntod !< total amount of ODs
    integer(i4b)                          :: ndet !< total amount of detectors
    ! HDF5 variables
    character(len=512) :: mystring, mysubstring !< dummy values for string manipulation
    integer(i4b)       :: myindex     !< dummy value for string manipulation
    character(len=512) :: currentHDFFile !< hdf5 file which stores simulation output
    character(len=6)   :: pidLabel
    character(len=3)   :: detectorLabel
    type(hdf_file)     :: hdf5_file   !< hdf5 file to work with
    integer(i4b)       :: hdf5_error  !< hdf5 error status
    integer(HID_T)     :: hdf5_file_id !< File identifier
    integer(HID_T)     :: dset_id     !< Dataset identifier
    integer(HSIZE_T), dimension(1) :: dims
    ! Other variables
    integer(i4b)                          :: i, j, k !< loop variables
    integer(i4b)       :: mpi_err !< MPI error status
    integer(i4b)       :: nomp !< Number of threads available
    integer(i4b)       :: omp_err !< OpenMP error status
    integer(i4b) :: omp_get_max_threads
    integer(i4b) :: n, nfft
    integer*8    :: plan_back
    real(sp) :: nu
    real(sp), allocatable, dimension(:,:) :: n_corr
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    character(len=10) :: processor_label   !< to have a nice output to screen

    ! shortcuts
    ntod = self%scans(scan_id)%ntod
    ndet = self%ndet

    ! Simulating 1/f noise
    !write(*,*) "Simulating correlated noise"
    nfft = 2 * ntod
    n = nfft / 2 + 1
    nomp = omp_get_max_threads()
    call sfftw_init_threads(omp_err)
    call sfftw_plan_with_nthreads(nomp)
    ! planning FFTW - in principle we should do both forward and backward FFTW,
    ! but in this case we can omit forward one and go directly with backward to
    ! save some time on a whole operation.
    allocate(dt(nfft), dv(0:n-1))
    call sfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    !$OMP PARALLEL PRIVATE(i, j, k, dt, dv, sigma0, nu)
    allocate(dt(nfft), dv(0:n-1), n_corr(ntod, ndet))
    !$OMP DO SCHEDULE(guided)
    do j = 1, ndet
      ! skipping iteration if scan was not accepted
      if (.not. self%scans(scan_id)%d(j)%accept) cycle
      ! getting gain for each detector (units, V / K)
      ! (gain is assumed to be CONSTANT for EACH SCAN)
      gain   = self%scans(scan_id)%d(j)%gain
      sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
      samprate = self%samprate
      ! used when adding fluctuation terms to Fourier coeffs (depends on Fourier convention)
      fft_norm = sqrt(1.d0 * nfft)
      !
      !dv(0) = dv(0) + fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0)
      dv(0) = 0. ! fft_norm * sigma0 * cmplx(rand_gauss(handle),rand_gauss(handle)) / sqrt(2.0) ! HKE: This expression is not correct for the monopole
      do k = 1, (n - 1)
        nu    = k * (samprate / 2) / (n - 1)
        N_c   = self%scans(scan_id)%d(j)%N_psd%eval_corr(nu)
        dv(k) = cmplx(rand_gauss(handle), rand_gauss(handle)) * sqrt(N_c) /sqrt(2.0)
      end do
      ! Executing Backward FFT
      call sfftw_execute_dft_c2r(plan_back, dv, dt)
      dt = dt / nfft
      n_corr(:, j) = dt(1:ntod)
      !write(*,*) "n_corr ", n_corr(:, j)
    end do
    !$OMP END DO
    deallocate(dt, dv)
    !$OMP END PARALLEL

    call sfftw_destroy_plan(plan_back)

    ! Allocating main simulations' array
    allocate(tod_per_detector(ntod, ndet))       ! Simulated tod

    ! Main simulation loop
    do i = 1, ntod
      do j = 1, ndet
        ! skipping iteration if scan was not accepted
        if (.not. self%scans(scan_id)%d(j)%accept) cycle
        ! getting gain for each detector (units, V / K)
        ! (gain is assumed to be CONSTANT for EACH SCAN)
        gain   = self%scans(scan_id)%d(j)%gain
        !write(*,*) "gain ", gain
        sigma0 = self%scans(scan_id)%d(j)%N_psd%sigma0
        !write(*,*) "sigma0 ", sigma0
        ! Simulating tods
        tod_per_detector(i,j) = gain * s_tot(i,j) + n_corr(i, j) + sigma0 * rand_gauss(handle)
        !tod_per_detector(i,j) = 0
      end do
    end do

    !----------------------------------------------------------------------------------
    ! Saving stuff to hdf file
    ! Getting the full path and name of the current hdf file to overwrite
    !----------------------------------------------------------------------------------
    mystring = trim(self%hdfname(scan_id))
    mysubstring = 'LFI_0'
    myindex = index(trim(mystring), trim(mysubstring))
    currentHDFFile = trim(self%sims_output_dir)//'/'//trim(mystring(myindex:))
    !write(*,*) "hdf5name "//trim(self%hdfname(scan_id))
    !write(*,*) "currentHDFFile "//trim(currentHDFFile)
    ! Converting PID number into string value
    call int2string(self%scanid(scan_id), pidLabel)
    call int2string(self%myid, processor_label)
    write(*,*) "Process: "//trim(processor_label)//" started writing PID: "//trim(pidLabel)//", into:"
    write(*,*) trim(currentHDFFile)
    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop

    dims(1) = ntod
    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_error)
    ! Open an existing file - returns hdf5_file_id
    call  h5fopen_f(currentHDFFile, H5F_ACC_RDWR_F, hdf5_file_id, hdf5_error)
    do j = 1, ndet
      detectorLabel = self%label(j)
      ! Open an existing dataset.
      call h5dopen_f(hdf5_file_id, trim(pidLabel)//'/'//trim(detectorLabel)//'/'//'tod', dset_id, hdf5_error)
      ! Write tod data to a dataset
      call h5dwrite_f(dset_id, H5T_IEEE_F32LE, tod_per_detector(:,j), dims, hdf5_error)
      ! Close the dataset.
      call h5dclose_f(dset_id, hdf5_error)
    end do
    ! Close the file.
    call h5fclose_f(hdf5_file_id, hdf5_error)
    ! Close FORTRAN interface.
    call h5close_f(hdf5_error)

    !write(*,*) "hdf5_error",  hdf5_error
    ! freeing memory up
    deallocate(n_corr, tod_per_detector)
    write(*,*) "Process:", self%myid, "finished writing PID: "//trim(pidLabel)//"."

    ! lastly, we need to copy an existing filelist.txt into simulation folder
    ! and change the pointers to new files
    !if (self%myid == 0) then
    !  call system("cp "//trim(filelist)//" "//trim(simsdir))
    !  !mystring = filelist
    !  !mysubstring = ".txt"
    !  !myindex = index(trim(mystring), trim(mysubstring))
    !end if

    ! For debugging
    !call MPI_Finalize(mpi_err)
    !stop
  end subroutine simulate_tod

end module comm_tod_driver_mod
