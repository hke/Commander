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
module comm_data_mod
  use comm_param_mod
  use comm_bp_mod
  use comm_noise_mod
  use comm_beam_mod
  use comm_map_mod
  use comm_tod_mod
  use comm_tod_LFI_mod
  use comm_tod_SPIDER_mod
  use comm_tod_WMAP_mod
  use comm_tod_LB_mod
  use locate_mod
  implicit none

  type comm_data_set
     character(len=512)           :: label, unit, comp_sens
     integer(i4b)                 :: period, id_abs
     logical(lgt)                 :: sample_gain
     real(dp)                     :: gain, gain_prior(2)
     character(len=128)           :: gain_comp
     integer(i4b)                 :: gain_lmin, gain_lmax
     integer(i4b)                 :: ndet
     character(len=128)           :: tod_type
     logical(lgt)                 :: pol_only

     class(comm_mapinfo), pointer :: info      => null()
     class(comm_map),     pointer :: map       => null()
     class(comm_map),     pointer :: map0      => null() !for TOD data if outputing to HDF
     class(comm_map),     pointer :: res       => null()
     class(comm_map),     pointer :: c_old     => null()
     class(comm_map),     pointer :: c_prop    => null()
     class(comm_map),     pointer :: mask      => null()
     class(comm_map),     pointer :: procmask  => null()
     class(comm_map),     pointer :: procmask2 => null()
     class(comm_map),     pointer :: gainmask  => null()
     class(comm_tod),     pointer :: tod       => null()
     class(comm_N),       pointer :: N         => null()
     class(B_ptr),         allocatable, dimension(:) :: B
     class(comm_bp_ptr),   allocatable, dimension(:) :: bp
     type(comm_B_bl_ptr),  allocatable, dimension(:) :: B_smooth
     type(comm_B_bl_ptr),  allocatable, dimension(:) :: B_postproc
     type(comm_N_rms_ptr), allocatable, dimension(:) :: N_smooth
   contains
     procedure :: RJ2data
     procedure :: chisq => get_chisq
     !procedure :: apply_proc_mask
  end type comm_data_set

  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data
  integer(i4b),        allocatable, dimension(:) :: ind_ds

  
contains

  subroutine initialize_data_mod(cpar, handle)
    implicit none
    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b)       :: i, j, n, nmaps, numband_tot, ierr
    character(len=512) :: dir, mapfile
    class(comm_N), pointer  :: tmp => null()
    class(comm_mapinfo), pointer :: info_smooth => null(), info_postproc => null()
    real(dp), allocatable, dimension(:)   :: nu
    real(dp), allocatable, dimension(:,:) :: regnoise, mask_misspix

    ! Read all data sets
    numband = count(cpar%ds_active)
    numband_tot = cpar%numband
    dir = trim(cpar%datadir) // '/'
    allocate(data(numband))
    n = 0
    do i = 1, numband_tot
       if (.not. cpar%ds_active(i)) cycle
       n                      = n+1
       data(n)%id_abs         = i
       data(n)%label          = cpar%ds_label(i)
       data(n)%period         = cpar%ds_period(i)
       data(n)%unit           = cpar%ds_unit(i)
       data(n)%sample_gain    = cpar%ds_sample_gain(i)
       data(n)%gain_comp      = cpar%ds_gain_calib_comp(i)
       data(n)%gain_prior     = cpar%ds_gain_prior(i,:)
       data(n)%gain_lmin      = cpar%ds_gain_lmin(i)
       data(n)%gain_lmax      = cpar%ds_gain_lmax(i)
       data(n)%comp_sens      = cpar%ds_component_sensitivity(i)
       data(n)%tod_type       = cpar%ds_tod_type(i)

       if (cpar%myid == 0 .and. cpar%verbosity > 0) &
            & write(*,fmt='(a,i5,a,a)') '  Reading data set ', i, ' : ', trim(data(n)%label)
       call update_status(status, "data_"//trim(data(n)%label))

       ! Initialize map structures
       nmaps = 1; if (cpar%ds_polarization(i)) nmaps = 3
       data(n)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & nmaps, cpar%ds_polarization(i))
       call get_mapfile(cpar, i, mapfile)
       data(n)%map  => comm_map(data(n)%info, trim(dir)//trim(mapfile), mask_misspix=mask_misspix)
       if (cpar%only_pol) data(n)%map%map(:,1) = 0.d0
       ! Read processing mask
       if (trim(cpar%ds_procmask) /= 'none') then
          data(n)%procmask => comm_map(data(n)%info, trim(cpar%datadir)//'/'//trim(cpar%ds_procmask), &
               & udgrade=.true.)
       end if
       data(n)%res  => comm_map(data(n)%map)
       call update_status(status, "data_map")

       if (data(n)%sample_gain) then
          ! Read calibration mask
          if (trim(cpar%ds_maskfile_calib(i)) /= 'fullsky') then
             data(n)%gainmask => comm_map(data(n)%info, trim(cpar%datadir)//'/'//trim(cpar%ds_maskfile_calib(i)), &
                  & udgrade=.true.)
          end if
       end if

       ! Initialize TOD structures
       data(n)%ndet = 0
       if (cpar%enable_TOD_analysis) then
          if (trim(data(n)%tod_type) == 'LFI') then
             data(n)%tod => comm_LFI_tod(cpar, i, data(n)%info, data(n)%tod_type)
             data(n)%ndet = data(n)%tod%ndet
          else if (trim(data(n)%tod_type) == 'WMAP') then
             data(n)%tod => comm_WMAP_tod(cpar, i, data(n)%info, data(n)%tod_type)
             data(n)%ndet = data(n)%tod%ndet
          else if (trim(data(n)%tod_type) == 'SPIDER') then
             data(n)%tod => comm_SPIDER_tod(cpar, i, data(n)%info)
             data(n)%ndet = data(n)%tod%ndet
          else if (trim(data(n)%tod_type) == 'LB') then
             data(n)%tod => comm_LB_tod(cpar, i, data(n)%info, data(n)%tod_type)
             data(n)%ndet = data(n)%tod%ndet
          else if (trim(cpar%ds_tod_type(i)) == 'none') then
          else
             write(*,*) 'Unrecognized TOD experiment type = ', trim(data(n)%tod_type)
             stop
          end if

          if (trim(cpar%ds_tod_type(i)) /= 'none') then
             data(n)%map0 => comm_map(data(n)%map) !copy the input map that has no added regnoise, for output to HDF
          end if
       end if

       ! Initialize beam structures
       allocate(data(n)%B(0:data(n)%ndet)) 
       select case (trim(cpar%ds_beamtype(i)))
       case ('b_l')
          data(n)%B(0)%p => comm_B_bl(cpar, data(n)%info, n, i)
          do j = 1, data(n)%ndet
             data(n)%B(j)%p => comm_B_bl(cpar, data(n)%info, n, i, fwhm=data(n)%tod%fwhm(j))
             ! MNG: I stripped mb_eff out of here to make it compile, if we need
             ! this ever we need to introduce it back in somehow
          end do
       case default
          call report_error("Unknown beam format: " // trim(cpar%ds_noise_format(i)))
       end select
       call update_status(status, "data_beam")
 
       ! Read default gain from instrument parameter file
       call read_instrument_file(trim(cpar%datadir)//'/'//trim(cpar%cs_inst_parfile), &
            & 'gain', cpar%ds_label(i), 1.d0, data(n)%gain)

       ! Initialize mask structures
       if (trim(cpar%ds_maskfile(i)) == 'fullsky') then
          data(n)%mask     => comm_map(data(n)%info)
          data(n)%mask%map =  1.d0
       else
          data(n)%mask  => comm_map(data(n)%info, trim(dir)//trim(cpar%ds_maskfile(i)))
          where (data(n)%mask%map > 0.5d0)
             data(n)%mask%map = 1.d0
          elsewhere
             data(n)%mask%map = 0.d0
          end where
       end if
       data(n)%mask%map = data(n)%mask%map * mask_misspix
       if (trim(cpar%ds_sourcemask) /= 'none') then
          call apply_source_mask(data(n)%mask, trim(cpar%datadir)//'/'//trim(cpar%ds_sourcemask), &
               & data(n)%B(0)%p%r_max)
       end if
       if (cpar%only_pol) data(n)%mask%map(:,1) = 0.d0
       call update_status(status, "data_mask")
       deallocate(mask_misspix)

       ! Initialize noise structures
       select case (trim(cpar%ds_noise_format(i)))
       case ('rms') 
          allocate(regnoise(0:data(n)%info%np-1,data(n)%info%nmaps))
          if (associated(data(n)%procmask)) then
             data(n)%N       => comm_N_rms(cpar, data(n)%info, n, i, 0, data(n)%mask, handle, regnoise, &
                  & data(n)%procmask)
          else
             data(n)%N       => comm_N_rms(cpar, data(n)%info, n, i, 0, data(n)%mask, handle, regnoise)
          end if
          data(n)%map%map = data(n)%map%map + regnoise  ! Add regularization noise
          deallocate(regnoise)
       case ('QUcov') 
          data(n)%N       => comm_N_QUcov(cpar, data(n)%info, n, i, 0, data(n)%mask, handle, regnoise, &
               & data(n)%procmask)
          data(n)%pol_only = .true.
       case default
          call report_error("Unknown noise format: " // trim(cpar%ds_noise_format(i)))
       end select
       data(n)%map%map  = data(n)%map%map * data(n)%mask%map
       data(n)%pol_only = data(n)%N%pol_only
       call update_status(status, "data_N")

       ! Initialize bandpass structures; 0 is full freq, i is detector
       allocate(data(n)%bp(0:data(n)%ndet))

       do j = 1, data(n)%ndet
          data(n)%bp(j)%p => comm_bp(cpar, n, i, detlabel=data(n)%tod%label(j))
       end do

       if (trim(cpar%ds_tod_type(i)) == 'none') then
          data(n)%bp(0)%p => comm_bp(cpar, n, i, detlabel=data(n)%label)
       else
          data(n)%bp(0)%p => comm_bp(cpar, n, i, subdets=cpar%ds_tod_dets(i))
       end if

!!$       if (trim(data(n)%label) == '070ds1' .or. trim(data(n)%label) == '070ds2' .or. trim(data(n)%label) == '070ds3') then
!!$          write(*,*) 'Check bp'
!!$          data(n)%bp(0)%p => comm_bp(cpar, n, i, '070')
!!$       else
!!$          
!!$       end if
       


       ! Initialize smoothed data structures
       allocate(data(n)%B_smooth(cpar%num_smooth_scales))
       allocate(data(n)%B_postproc(cpar%num_smooth_scales))
       allocate(data(n)%N_smooth(cpar%num_smooth_scales))
       do j = 1, cpar%num_smooth_scales
          if (cpar%fwhm_smooth(j) > 0.d0) then
             info_smooth => comm_mapinfo(data(n)%info%comm, data(n)%info%nside, cpar%lmax_smooth(j), &
                  & data(n)%info%nmaps, data(n)%info%pol)
             data(n)%B_smooth(j)%p => &
               & comm_B_bl(cpar, info_smooth, n, i, fwhm=cpar%fwhm_smooth(j), pixwin=cpar%pixwin_smooth(j), &
               & init_realspace=.false.)
          else
             nullify(data(n)%B_smooth(j)%p)
          end if
          if (cpar%fwhm_postproc_smooth(j) > 0.d0) then
             info_postproc => comm_mapinfo(data(n)%info%comm, cpar%nside_smooth(j), cpar%lmax_smooth(j), &
                  & data(n)%info%nmaps, data(n)%info%pol)
             data(n)%B_postproc(j)%p => &
                  & comm_B_bl(cpar, info_postproc, n, i, fwhm=cpar%fwhm_postproc_smooth(j),&
                  & init_realspace=.false.)
          else
             nullify(data(n)%B_postproc(j)%p)
          end if
          if (trim(cpar%ds_noise_rms_smooth(i,j)) == 'native') then
             tmp => data(n)%N
             select type (tmp)
             class is (comm_N_rms)
                data(n)%N_smooth(j)%p => tmp
             end select
          else if (trim(cpar%ds_noise_rms_smooth(i,j)) /= 'none') then
             data(n)%N_smooth(j)%p => comm_N_rms(cpar, data(n)%info, n, i, j, data(n)%mask, handle)
          else
             nullify(data(n)%N_smooth(j)%p)
          end if
       end do

    end do
    !numband = n

    ! Sort bands according to nominal frequency
    allocate(ind_ds(numband), nu(numband))
    do i = 1, numband
       ind_ds(i)  = i
       nu(i)        = data(i)%bp(0)%p%nu_c
    end do
    call QuickSort(ind_ds, nu)
    deallocate(nu)

    ! Dump unit conversion factors to file
    if (cpar%myid == 0) call dump_unit_conversion(cpar%outdir)

  end subroutine initialize_data_mod

  function get_chisq(self)
    implicit none
    class(comm_data_set), intent(in)           :: self
    real(dp)                                   :: get_chisq

    integer(i4b) :: ierr
    real(dp)     :: chisq
    class(comm_map), pointer :: invN_res => null()
    
    invN_res => comm_map(self%res)
    call self%N%invN(invN_res)
    chisq = sum(self%res%map*invN_res%map)
    call mpi_allreduce(chisq, get_chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%info%comm, ierr)
    call invN_res%dealloc(); deallocate(invN_res)

  end function get_chisq

  function RJ2data(self, det)
    implicit none

    class(comm_data_set), intent(in)           :: self
    integer(i4b),         intent(in), optional :: det
    real(dp)                                   :: RJ2data

    integer(i4b) :: d

    d = 0; if (present(det)) d = det

    select case (trim(self%unit))
    case ('uK_cmb') 
       RJ2data = self%bp(d)%p%a2t
    case ('mK_cmb') 
       RJ2data = self%bp(d)%p%a2t * 1d-3
    case ('K_cmb') 
       RJ2data = self%bp(d)%p%a2t * 1d-6
    case ('MJy/sr') 
       RJ2data = self%bp(d)%p%a2t / self%bp(d)%p%f2t
    case ('y_SZ') 
       RJ2data = self%bp(d)%p%a2sz
    case ('uK_RJ') 
       RJ2data = 1.d0
    case ('K km/s') 
       RJ2data = 1.d0
    case default
       RJ2data = 1.d0
    end select
    
  end function RJ2data

  subroutine dump_unit_conversion(dir)
    implicit none
    character(len=*), intent(in) :: dir
    integer(i4b) :: i, q, unit

    unit = getlun()
    open(unit, file=trim(dir)//'/unit_conversions.dat', recl=1024)
    write(unit,*) '# Band   BP type   Nu_c (GHz)  a2t [K_cmb/K_RJ]' // &
         & '  t2f [MJy/K_cmb] a2sz [y_sz/K_RJ]'
    do i = 1, numband
       q = ind_ds(i)
       write(unit,fmt='(a7,a10,f10.1,3e16.5)') trim(data(q)%label), trim(data(q)%bp(0)%p%type), &
            & data(q)%bp(0)%p%nu_c/1.d9, data(q)%bp(0)%p%a2t, 1.d0/data(q)%bp(0)%p%f2t*1e6, data(q)%bp(0)%p%a2sz * 1.d6
    end do
    close(unit)
  end subroutine dump_unit_conversion

  subroutine apply_source_mask(mask, sourcefile, r_max)
    implicit none
    class(comm_map),  intent(inout) :: mask
    character(len=*), intent(in)    :: sourcefile
    real(dp),         intent(in)    :: r_max

    integer(i4b)       :: i, j, l, unit, nlist, nmax=10000, itmp
    real(dp)           :: lon, lat, rad, vec(3)
    character(len=512) :: line
    integer(i4b), allocatable, dimension(:) :: listpix

    unit = getlun()
    open(unit,file=trim(sourcefile),recl=1024,status='old')
    do while (.true.)
       read(unit,fmt='(a)',end=23) line
       line = trim(adjustl(line))
       if (line(1:1) == '#' .or. line(1:1) == ' ') cycle
       read(line,*) lon, lat, rad
       call ang2vec(0.5d0*pi-lat*DEG2RAD, lon*DEG2RAD, vec)
       allocate(listpix(0:nmax-1))
       call query_disc(mask%info%nside, vec, r_max*rad, listpix, nlist)

       ! Sort pixel list according to increasing pixel number
       do j = 1, nlist-1
          itmp          = listpix(j)
          l             = j-1
          do while (l >= 0)
             if (listpix(l) <= itmp) exit
             listpix(l+1) = listpix(l)
             l            = l-1
          end do
          listpix(l+1)    = itmp
       end do

       ! Mask pixels belonging to current processor
       i    = 0
       j    = locate(mask%info%pix, listpix(i))
       if (j > 0) then
          do while (.true.)
             if (listpix(i) == mask%info%pix(j)) then
                mask%map(j,:) = 0.d0
                i             = i+1
                j             = j+1
             else if (listpix(i) < mask%info%pix(j)) then
                i               = i+1
             else
                j               = j+1
             end if
             if (i > nlist-1) exit
             if (j > mask%info%np) exit
          end do
       end if

       deallocate(listpix)
    end do
23  close(unit)

  end subroutine apply_source_mask

!!$  subroutine smooth_inside_procmask(data, fwhm)
!!$    implicit none
!!$    class(comm_data_set), intent(in) :: data
!!$    real(dp),             intent(in) :: fwhm
!!$
!!$    integer(i4b) :: i, j
!!$    real(dp)     :: w
!!$    class(comm_mapinfo), pointer :: info
!!$    class(comm_map),     pointer :: map
!!$
!!$    map => comm_map(data%map)
!!$    call map%smooth(fwhm)
!!$    do j = 1, data%map%info%nmaps
!!$       do i = 0, data%map%info%np-1
!!$          w = data%procmask%map(i,j)
!!$          data%map%map(i,j) = w * data%map%map(i,j) + (1.d0-w) * map%map(i,j)
!!$       end do
!!$    end do
!!$    deallocate(map)
!!$
!!$  end subroutine smooth_inside_procmask

!!$  subroutine apply_proc_mask(self, map)
!!$    implicit none
!!$    class(comm_data_set), intent(in)    :: self
!!$    class(comm_map),      intent(inout) :: map
!!$
!!$    if (.not. associated(self%procmask)) return
!!$    where (self%procmask%map < 1.d0)  ! Apply processing mask
!!$       map%map = -1.6375d30
!!$    end where
!!$
!!$  end subroutine apply_proc_mask

  subroutine smooth_map(info, alms_in, bl_in, map_in, bl_out, map_out)
    implicit none
    class(comm_mapinfo),                      intent(in),   target :: info
    logical(lgt),                             intent(in)           :: alms_in
    real(dp),            dimension(0:,1:),    intent(in)           :: bl_in, bl_out
    class(comm_map),                          intent(inout)        :: map_in
    class(comm_map),                          intent(out), pointer :: map_out

    integer(i4b) :: i, j, l, lmax

    map_out => comm_map(info)

    if (.not. alms_in) then
       !map_out%map = map_in%map
       call map_in%udgrade(map_out)
       call map_out%YtW
    else
       call map_in%alm_equal(map_out)
    end if


    ! Deconvolve old beam, and convolve with new beam
    lmax  = min(size(bl_in,1)-1, size(bl_out,1)-1)
    do i = 0, info%nalm-1
       l = info%lm(1,i)
       if (l > lmax) then
          map_out%alm(i,:) = 0.d0
          cycle
       end if
       do j = 1, map_out%info%nmaps
          if (bl_in(l,j) > 1.d-12) then
             map_out%alm(i,j) = map_out%alm(i,j) * bl_out(l,j) / bl_in(l,j)
          else
             map_out%alm(i,j) = 0.d0
          end if
       end do
    end do    

    ! Recompose map
    call map_out%Y

  end subroutine smooth_map

  subroutine get_mapfile(cpar, band, mapfile)
    implicit none
    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: band
    character(len=*),  intent(out)   :: mapfile

    integer(i4b)       :: i, n, unit
    character(len=512) :: filename

    filename = trim(adjustl(cpar%ds_mapfile(band)))
    n        = len(trim(adjustl(filename)))
    
    if (filename(n-3:n) == 'fits') then
       mapfile = filename
    else if (filename(n-2:n) == 'txt') then
       ! Assume filelist; pick the number given by mychain. 
       unit = getlun()
       open(unit, file=trim(cpar%datadir)//'/'//trim(filename), recl=1024)
       do i = 1, cpar%mychain
          read(unit,'(a)') mapfile
       end do
       close(unit)
    end if

  end subroutine get_mapfile

end module comm_data_mod
