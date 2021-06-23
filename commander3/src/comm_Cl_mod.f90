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
module comm_Cl_mod
  use comm_param_mod
  use comm_map_mod
  use math_tools
  use comm_hdf_mod
  use comm_bp_utils
  use InvSamp_mod
  implicit none

  private
  public comm_Cl, write_sigma_l

  integer(i4b), parameter :: TT = 1
  integer(i4b), parameter :: TE = 2
  integer(i4b), parameter :: TB = 3
  integer(i4b), parameter :: EE = 4
  integer(i4b), parameter :: EB = 5
  integer(i4b), parameter :: BB = 6
  logical(lgt), private   :: only_pol

  type :: Cl_bin
     integer(i4b)     :: lmin, lmax, spec, nsub, ntot, p1, p2
     real(dp)         :: sigma
     character(len=1) :: stat
     type(Cl_bin), pointer, dimension(:)   :: sub
     real(dp),     allocatable, dimension(:,:) :: M_prop
  end type Cl_bin
  
  type :: comm_Cl
     ! General parameters
     class(comm_mapinfo), pointer :: info
     character(len=512)           :: type  ! {none, binned, power_law, exp}
     character(len=512)           :: label ! {none, binned, power_law, exp}
     character(len=512)           :: unit
     character(len=512)           :: outdir
     integer(i4b)                 :: lmin, lmax, nmaps, nspec, l_apod, lmax_prior
     integer(i4b)                 :: lmin_lookup, lmax_lookup
     integer(i4b)                 :: poltype  ! {1 = {T+E+B}, 2 = {T,E+B}, 3 = {T,mE,B}}
     logical(lgt)                 :: only_pol, active_lookup(6)
     real(dp)                     :: nu_ref(3), RJ2unit(3)
     real(dp),         allocatable, dimension(:,:)   :: Dl
     real(dp),         allocatable, dimension(:)     :: par_lookup
     real(dp),         allocatable, dimension(:,:,:) :: Dl_lookup
     real(dp),         allocatable, dimension(:,:,:) :: sqrtS_mat, S_mat, sqrtInvS_mat

     ! Bin parameters
     integer(i4b) :: nbin
     character(len=1), allocatable, dimension(:,:) :: stat
     integer(i4b), allocatable, dimension(:,:)     :: bins

     ! Variable binsize structure
     integer(i4b) :: nbin2
     type(Cl_bin), allocatable, dimension(:)     :: bins2
     
     ! Power law/exponential parameters
     character(len=512) :: plfile
     integer(i4b)       :: iter = 0
     real(dp)           :: prior(2), lpiv
     real(dp), allocatable, dimension(:) :: amp, beta, theta
   contains
     ! Data procedures
     procedure :: S        => matmulS
     procedure :: sqrtS    => matmulSqrtS
     procedure :: sqrtInvS => matmulSqrtInvS
     procedure :: sampleCls
     procedure :: read_binfile
     procedure :: read_binfile2
     procedure :: binCls
     procedure :: binCls2
     procedure :: binCl2
     procedure :: writeFITS
     procedure :: initHDF
     procedure :: updatePowlaw
     procedure :: updatePowlawGauss
     procedure :: updateExponential
     procedure :: updateGaussian
     procedure :: updateS
     procedure :: getCl
     procedure :: set_Dl_bin
     procedure :: check_posdef
  end type comm_Cl

  interface comm_Cl
     procedure constructor
  end interface comm_Cl


contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs)
    implicit none
    type(comm_params),                  intent(in) :: cpar
    type(comm_mapinfo), target,         intent(in) :: info
    integer(i4b),                       intent(in) :: id, id_abs
    class(comm_Cl),                     pointer    :: constructor

    integer(i4b)       :: l, nmaps
    character(len=512) :: datadir, binfile
    logical(lgt)       :: pol

    allocate(constructor)
    
    ! General parameters
    constructor%type   = cpar%cs_cltype(id_abs)
    if (trim(constructor%type) == 'none') return
    
    constructor%info   => info
    constructor%label  = cpar%cs_label(id_abs)
    constructor%lmin   = cpar%cs_lmin_amp(id_abs)
    constructor%lmax   = cpar%cs_lmax_amp(id_abs)
    constructor%lmax_prior = cpar%cs_lmax_amp_prior(id_abs)
    constructor%unit   = cpar%cs_unit(id_abs)
    constructor%nu_ref = cpar%cs_nu_ref(id_abs,:)
    constructor%nmaps  = 1; if (cpar%cs_polarization(id_abs)) constructor%nmaps = 3
    constructor%nspec  = 1; if (cpar%cs_polarization(id_abs)) constructor%nspec = 6
    constructor%outdir = cpar%outdir
    datadir            = cpar%datadir
    nmaps              = constructor%nmaps
    constructor%only_pol = cpar%only_pol
    constructor%lmin_lookup = -1
    constructor%lmax_lookup = -1

    ! Set up conversion factor between RJ and native component unit
    ! D_l is defined in component units, while S, invS etc are defined in RJ
    do l = 1, 3
       select case (trim(constructor%unit))
       case ('uK_cmb')
          constructor%RJ2unit(l) = comp_a2t(constructor%nu_ref(l))
       case ('mK_cmb')
          constructor%RJ2unit(l) = comp_a2t(constructor%nu_ref(l)) * 1d-3
       case ('K_cmb')
          constructor%RJ2unit(l) = comp_a2t(constructor%nu_ref(l)) * 1d-6
       case ('MJy/sr') 
          constructor%RJ2unit(l) = comp_bnu_prime_RJ(constructor%nu_ref(l)) * 1e14
       case ('K km/s') 
          constructor%RJ2unit(l) = 1.d0
       case ('y_SZ') 
          constructor%RJ2unit(l) = 2.d0*constructor%nu_ref(l)**2*k_b/c**2 / &
               & (comp_bnu_prime(constructor%nu_ref(l)) * comp_sz_thermo(constructor%nu_ref(l)))
       case ('uK_RJ') 
          constructor%RJ2unit(l) = 1.d0
       case default
          call report_error('Unsupported unit: ' // trim(constructor%unit))
       end select
    end do

    allocate(constructor%Dl(0:constructor%lmax,constructor%nspec))
    allocate(constructor%sqrtS_mat(nmaps,nmaps,0:constructor%lmax))
    allocate(constructor%sqrtInvS_mat(nmaps,nmaps,0:constructor%lmax))
    allocate(constructor%S_mat(nmaps,nmaps,0:constructor%lmax))

    if (trim(constructor%type) == 'binned') then
       !call constructor%read_binfile(trim(datadir) // '/' // trim(cpar%cs_binfile(id_abs)))
       call constructor%read_binfile2(datadir, cpar%cs_binfile(id_abs))
       call read_Cl_file(trim(datadir) // '/' // trim(cpar%cs_clfile(id_abs)), &
            & constructor%Dl, 0, constructor%lmax, 'TT_TE_EE_BB')
       call constructor%binCls2
       if (cpar%only_pol) then
          constructor%stat(:,1:3) = '0'
          constructor%Dl(:,TT)     = 0.d0
          if (constructor%nspec == 6) constructor%Dl(:,[TE,TB]) = 0.d0
       end if
       write(*,*) "lmax:", constructor%lmax
       write(*,*) constructor%Dl(:10,BB)
    else if (trim(constructor%type) == 'power_law') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updatePowlaw(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps))
    else if (trim(constructor%type) == 'power_law_gauss') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps),constructor%theta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updatePowlawGauss(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps),cpar%cs_cl_theta_def(id_abs,1:nmaps))
    else if (trim(constructor%type) == 'exp') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updateExponential(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps))
    else if (trim(constructor%type) == 'gauss') then
       allocate(constructor%amp(nmaps), constructor%beta(nmaps))
       constructor%lpiv          = cpar%cs_lpivot(id_abs)
       constructor%prior         = cpar%cs_cl_prior(id_abs,:)
       constructor%poltype       = cpar%cs_cl_poltype(id_abs)
       call constructor%updateGaussian(cpar%cs_cl_amp_def(id_abs,1:nmaps), cpar%cs_cl_beta_def(id_abs,1:nmaps))
    else
       call report_error("Unknown Cl type: " // trim(constructor%type))
    end if
    if (constructor%only_pol) then
       constructor%Dl(:,TT)     = 0.d0
       if (constructor%nspec == 6) constructor%Dl(:,[TE,TB]) = 0.d0
    end if
    constructor%Dl(0:1,2:constructor%nspec) = 0.d0
    constructor%Dl(0:constructor%lmin-1,:)  = 0.d0
    !constructor%Dl = constructor%Dl / c%RJ2unit_    ! Define prior in output units
    call constructor%updateS
    
  end function constructor


  subroutine updatePowlaw(self, amp, beta)
    implicit none
    class(comm_Cl),                intent(inout)        :: self
    real(dp),       dimension(1:), intent(in), optional :: amp, beta

    integer(i4b) :: i, j, l, i_min
    
    if (present(amp))  self%amp   = amp
    if (present(beta)) self%beta  = beta
    self%Dl    = 0.d0
    i_min      = 1; if (self%only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = max(self%lmin,1), self%lmax
          if (i > 1 .and. l < 2) cycle
          self%Dl(l,j) = self%amp(i) * (real(l,dp)/real(self%lpiv,dp))**self%beta(i) 
       end do
       self%Dl(0,j) = self%Dl(1,j)
    end do

  end subroutine updatePowlaw

  subroutine updatePowlawGauss(self, amp, beta, theta)
    implicit none
    class(comm_Cl),                intent(inout)        :: self
    real(dp),       dimension(1:), intent(in), optional :: amp, beta, theta

    integer(i4b) :: i, j, l, i_min
    
    if (present(amp))  self%amp   = amp
    if (present(beta)) self%beta  = beta
    if (present(theta)) self%theta  = theta
    self%Dl    = 0.d0
    i_min      = 1; if (self%only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = max(self%lmin,1), self%lmax
          if (i > 1 .and. l < 2) cycle
          self%Dl(l,j) = self%amp(i) * (real(l,dp)/real(self%lpiv,dp))**self%beta(i) * max(exp(-l*(l+1)*(self%theta(i)*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2),1d-10)
       end do
       self%Dl(0,j) = self%Dl(1,j)
    end do

  end subroutine updatePowlawGauss

  subroutine updateExponential(self, amp, beta)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    real(dp),       dimension(1:), intent(in), optional    :: amp, beta

    integer(i4b) :: i, j, l, i_min
    
    if (present(amp))  self%amp   = amp
    if (present(beta)) self%beta  = beta
    self%Dl    = 0.d0
    i_min      = 1; if (self%only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = max(self%lmin,1), self%lmax
          if (i > 1 .and. l < 2) cycle
          self%Dl(l,j) = self%amp(i) * exp(-self%beta(i)*(real(l,dp)/real(self%lpiv,dp))) 
       end do
       if (self%lmin <= 0) self%Dl(0,j) = self%Dl(1,j)
    end do

  end subroutine updateExponential

  subroutine updateGaussian(self, amp, beta)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    real(dp),       dimension(1:), intent(in), optional :: amp, beta

    integer(i4b) :: i, j, l, i_min
    
    if (present(amp))  self%amp   = amp
    if (present(beta)) self%beta  = beta
    self%Dl    = 0.d0
    i_min      = 1; if (self%only_pol) i_min = 2
    do i = i_min, self%nmaps
       j = i*(1-i)/2 + (i-1)*self%nmaps + i
       do l = max(self%lmin,0), self%lmax
          if (i > 1 .and. l < 2) cycle
          self%Dl(l,j) = self%amp(i) * max(exp(-l*(l+1)*(self%beta(i)*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2),1d-10)
!          if (self%info%myid == 0) then
!             write(*,*) l, j, amp(i), beta(i), exp(-l*(l+1)*(beta(i)*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2), self%Dl(l,j)
!          end if
       end do
    end do

  end subroutine updateGaussian

  subroutine updateS(self)
    implicit none
    class(comm_Cl),                intent(inout) :: self

    integer(i4b) :: i, j, k, l
    logical(lgt) :: ok(self%nmaps)

    self%S_mat        = 0.d0
    self%sqrtS_mat    = 0.d0
    self%sqrtInvS_mat = 0.d0
    do l = 0, self%lmax
       k = 1
       do i = 1, self%nmaps
          do j = i, self%nmaps
             !k = i*(1-i)/2 + (i-1)*self%nmaps + j
             if (l < self%lmin) then
                self%sqrtS_mat(i,j,l) = 0
             else if (l == 0) then
                self%sqrtS_mat(i,j,l) = self%Dl(l,k)
             else
                self%sqrtS_mat(i,j,l) = self%Dl(l,k) / (l*(l+1)/(2.d0*pi))
             end if
             self%sqrtS_mat(i,j,l) = self%sqrtS_mat(i,j,l) / &
                  & (self%RJ2unit(i)*self%RJ2unit(j))
             self%sqrtS_mat(j,i,l) = self%sqrtS_mat(i,j,l) 
             if (i == j) ok(i) = self%Dl(l,k) > 0.d0
             k = k+1
          end do
       end do
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtS_mat(i,:,l) = 0.d0
             self%sqrtS_mat(:,i,l) = 0.d0
             self%sqrtS_mat(i,i,l) = 1.d0
          end if
       end do

!       if (self%info%myid == 0) write(*,*) l, self%sqrtS_mat(:,:,l)
       !if (self%info%myid == 0) write(*,*) l, self%sqrtS_mat(1,1,l)*self%sqrtS_mat(2,2,l) - self%sqrtS_mat(2,1,l)*self%sqrtS_mat(1,2,l)

       ! Change to RJ units
       !self%sqrtS_mat(:,:,l) = self%sqrtS_mat(:,:,l) / self%RJ2unit**2

       self%sqrtInvS_mat(:,:,l) = self%sqrtS_mat(:,:,l)
       call compute_hermitian_root(self%sqrtS_mat(:,:,l), 0.5d0)
!       if (self%info%myid == 0) write(*,*) l, self%sqrtS_mat(:,:,l), self%sqrtInvS_mat(:,:,l)
       
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtS_mat(i,:,l) = 0.d0
             self%sqrtS_mat(:,i,l) = 0.d0
          end if
       end do
       self%S_mat(:,:,l) = matmul(self%sqrtS_mat(:,:,l), self%sqrtS_mat(:,:,l))
       
       call compute_hermitian_root(self%sqrtInvS_mat(:,:,l), -0.5d0)
       do i = 1, self%nmaps
          if (.not. ok(i)) then
             self%sqrtInvS_mat(i,:,l) = 0.d0
             self%sqrtInvS_mat(:,i,l) = 0.d0
          end if
       end do

    end do

!!$    call mpi_finalize(i)
!!$    stop

  end subroutine updateS

  subroutine read_binfile(self, binfile)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)  :: binfile

    integer(i4b)       :: unit, l1, l2, n
    character(len=512) :: line
    character(len=1), dimension(6) :: stat

    unit = getlun()

    ! Find number of bins
    n = 0
    open(unit,file=trim(binfile))
    do while (.true.)
       read(unit,'(a)',end=1) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') cycle
       read(line,*) l1, l2
       if (l1 <= self%lmax .and. l2 <= self%lmax .and. l1 >= 0 .and. l2 >= 0) n = n+1
    end do
1   close(unit)

    if (n == 0) call report_error("Error: Cl bin file " // trim(binfile) &
         & // " does not contain any valid entries")

    self%nbin = n
    allocate(self%stat(n,6), self%bins(n,2))
    open(unit,file=trim(binfile))
    n = 0
    do while (.true.)
       read(unit,'(a)',end=2) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') cycle
       read(line,*) l1, l2
       if (l1 <= self%lmax .and. l2 <= self%lmax .and. l1 >= 0 .and. l2 >= 0) then
          n = n+1
          read(line,*) l1, l2, stat   !self%bins(n,1), self%bins(n,2), self%stat(n,:)
          self%bins(n,1) = l1
          self%bins(n,2) = l2
          self%stat(n,:) = stat
       end if
    end do
2   close(unit)

  end subroutine read_binfile

  recursive subroutine read_bin(unit, bin, ntot, pos, sigma)
    implicit none
    integer(i4b), intent(in)     :: unit
    type(Cl_bin), intent(inout)  :: bin
    integer(i4b), intent(out)    :: ntot
    integer(i4b), intent(inout)  :: pos
    real(dp), dimension(:), intent(inout) :: sigma

    integer(i4b)     :: i, nsub
    character(len=2) :: spec

    pos = pos+1
    read(unit,*) spec, bin%lmin, bin%lmax, bin%nsub, bin%stat, bin%sigma
    if (bin%sigma <= 0.d0 .and. (bin%stat == 'M' .or. bin%stat == 'S')) then
       write(*,*) 'Error -- sigma = 0 for active bin in binfile'
       stop
    end if
    sigma(pos) = bin%sigma
    select case (spec)
    case ('TT')
       bin%spec = 1; bin%p1 = 1; bin%p2 = 1
    case ('TE')
       bin%spec = 2; bin%p1 = 1; bin%p2 = 2
    case ('TB')
       bin%spec = 3; bin%p1 = 1; bin%p2 = 3
    case ('EE')
       bin%spec = 4; bin%p1 = 2; bin%p2 = 2
    case ('EB')
       bin%spec = 5; bin%p1 = 2; bin%p2 = 3
    case ('BB')
       bin%spec = 6; bin%p1 = 3; bin%p2 = 3
    end select

    ntot = 0; if (bin%stat == 'M') ntot = 1
    if (bin%nsub == 0) return

    allocate(bin%sub(bin%nsub))
    do i = 1, bin%nsub
       call read_bin(unit, bin%sub(i), nsub, pos, sigma)
       ntot = ntot + nsub
    end do

  end subroutine read_bin

  subroutine read_binfile2(self, datadir, binfile)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)    :: datadir, binfile

    integer(i4b)       :: unit, unit2, l1, l2, l, k, n, i, j, pos, nspline
    logical(lgt)       :: flag(6)
    real(dp), dimension(1000)           :: sigma
    real(dp), allocatable, dimension(:)     :: p_in
    real(dp), allocatable, dimension(:,:,:) :: Dl_in
    character(len=2)   :: spec
    character(len=512) :: line, filename
    character(len=1), dimension(6) :: stat
    type(spline_type) :: Dl_spline

    unit    = getlun()
    nspline = 100
    open(unit,file=trim(datadir)//'/'//trim(binfile))

    ! Check for parameter lookup table
    read(unit,'(a)') line
    if (trim(line) /= 'none') then
       read(line,*) filename, self%lmin_lookup, self%lmax_lookup, self%active_lookup, n
       allocate(p_in(n), Dl_in(self%lmin_lookup:self%lmax_lookup,6,n))
       allocate(self%par_lookup(nspline*(n-1)+1), self%Dl_lookup(self%lmin_lookup:self%lmax_lookup,6,nspline*(n-1)+1))
       unit2 = getlun()
       open(unit2,file=trim(datadir)//'/'//trim(filename))
       do i = 1, n
          read(unit2,*) p_in(i), filename
          call read_Cl_file(trim(datadir) // '/' // trim(filename), &
            & Dl_in(:,:,i), self%lmin_lookup, self%lmax_lookup, 'TT_EE_BB_TE')
          if (i == n/2) then
             do l1 = self%lmin_lookup, self%lmax_lookup
                where (self%active_lookup) 
                   self%Dl(l1,:) = Dl_in(l1,:,i) 
                end where
             end do
          end if
       end do
       close(unit2)
       !write(*,*) "Dl_in:", self%Dl
       ! Spline spectra
       do j = 1, 6
          if (.not. self%active_lookup(j)) cycle
          do l = self%lmin_lookup, self%lmax_lookup
             call spline(Dl_spline, p_in, Dl_in(l,j,:))
             do k = 1, nspline*(n-1)+1
                self%par_lookup(k) = p_in(1) + (p_in(n)-p_in(1)) * (k-1)/real(nspline*(n-1),dp)
                self%Dl_lookup(l,j,k) = splint(Dl_spline, self%par_lookup(k))
             end do
             call free_spline(Dl_spline)
          end do
       end do
    end if

    ! Find number of bins
    read(unit,*) self%nbin2
    allocate(self%bins2(self%nbin2))
    do i = 1, self%nbin2
       pos = 0
       call read_bin(unit, self%bins2(i), n, pos, sigma)
       self%bins2(i)%ntot = n
       allocate(self%bins2(i)%M_prop(n,n))
       self%bins2(i)%M_prop = 0.d0
       do j = 1, n
          self%bins2(i)%M_prop(j,j) = sigma(j)
       end do
    end do
    close(unit)

  end subroutine read_binfile2


  subroutine matmulS(self, map, alm, info)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    class(comm_mapinfo),               intent(in),    optional :: info
    integer(i4b) :: i, l
    real(dp)     :: f_apod
    if (trim(self%type) == 'none') then
       if (only_pol) then
          if (present(map)) then
             map%alm(:,1) = 0.d0
          else if (present(alm)) then
             alm(:,1) = 0.d0
          end if
       end if
       return
    end if

    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .true.)
          map%alm(i,:) = f_apod**2 * matmul(self%S_mat(:,:,l), map%alm(i,:))
       end do
    else if (present(alm)) then
       do i = 0, info%nalm-1
          l = info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .true.)
          if (l <= self%lmax) then
             alm(i,:) = f_apod**2 * matmul(self%S_mat(:,:,l), alm(i,:))
          else
             alm(i,:) = 0.d0
          end if
       end do
    end if
  end subroutine matmulS

  subroutine matmulSqrtS(self, map, alm, info, diag)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    class(comm_mapinfo),               intent(in),    optional :: info
    logical(lgt),                      intent(in),    optional :: diag
    integer(i4b) :: i, l, j
    real(dp)     :: f_apod
    if (trim(self%type) == 'none') then
       if (only_pol) then
          if (present(map)) then
             map%alm(:,1) = 0.d0
          else if (present(alm)) then
             alm(:,1) = 0.d0
          end if
       end if
       return
    end if
    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .true.)
          if (present(diag)) then
             do j = 1, self%info%nmaps
                map%alm(i,j) = f_apod * sqrt(self%S_mat(j,j,l)) * map%alm(i,j)
             end do
          else
             map%alm(i,:) = f_apod * matmul(self%sqrtS_mat(:,:,l), map%alm(i,:))
          end if
       end do
    else if (present(alm)) then
       do i = 0, info%nalm-1
          l = info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .true.)
          if (l <= self%lmax) then
             if (present(diag)) then
                do j = 1, self%info%nmaps
                   alm(i,j) = f_apod * sqrt(self%S_mat(j,j,l)) * alm(i,j)
                end do
             else
                alm(i,:) = f_apod * matmul(self%sqrtS_mat(:,:,l), alm(i,:))
             end if
          else
             alm(i,:) = 0.d0
          end if
          !if (self%info%myid == 0 .and. i == 4) write(*,*) l, real(f_apod,sp), real(self%sqrtS_mat(:,:,l),sp),  real(alm(i,:),sp), self%Dl(l,:), trim(self%label)
       end do
    end if
  end subroutine matmulSqrtS

  subroutine matmulSqrtInvS(self, map, alm, info)
    implicit none
    class(comm_Cl),                    intent(in)              :: self
    class(comm_map),                   intent(inout), optional :: map
    real(dp),        dimension(0:,1:), intent(inout), optional :: alm
    class(comm_mapinfo),               intent(in),    optional :: info
    integer(i4b) :: i, l
    real(dp)     :: f_apod
    if (trim(self%type) == 'none') then
       if (only_pol) then
          if (present(map)) then
             map%alm(:,1) = 0.d0
          else if (present(alm)) then
             alm(:,1) = 0.d0
          end if
       end if
       return
    end if
    if (present(map)) then
       do i = 0, self%info%nalm-1
          l = self%info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .false.)
          map%alm(i,:) = f_apod * matmul(self%sqrtInvS_mat(:,:,l), map%alm(i,:))
       end do
    else if (present(alm)) then
       do i = 0, info%nalm-1
          l = info%lm(1,i)
          f_apod = get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .false.)
          if (l <= self%lmax) then
             alm(i,:) = f_apod * matmul(self%sqrtInvS_mat(:,:,l), alm(i,:))
          else
             alm(i,:) = 0.d0
          end if
       end do
    end if
  end subroutine matmulSqrtInvS

  function get_Cl_apod(l, l_apod, lmax, lmax_prior, positive)
    implicit none
    integer(i4b), intent(in) :: l, l_apod, lmax, lmax_prior 
    logical(lgt), intent(in) :: positive
    real(dp)                 :: get_Cl_apod
    real(dp), parameter :: alpha = log(1d3)
    if (l_apod > 0) then
       if (l <= l_apod) then
          get_Cl_apod = 1.d0
       else if (l > lmax) then
          get_Cl_apod = 0.d0
       else
          get_Cl_apod = exp(-alpha * (l-l_apod)**2 / real(lmax-l_apod+1,dp)**2)
       end if
    else
       if (l >= abs(l_apod)) then
          get_Cl_apod = 1.d0
       else if (l == 0 .or. l > lmax) then
          get_Cl_apod = 0.d0
       else
          get_Cl_apod = exp(-alpha * (abs(l_apod)-l)**2 / real(abs(l_apod)-1,dp)**2)
       end if
    end if
    if (lmax_prior >= 0 .and. l < lmax_prior) then
       ! Apply cosine apodization between 0 and lmax_prior
       get_Cl_apod = get_Cl_apod * (0.5d0*(cos(pi*real(max(l,1)-lmax_prior,dp)/real(lmax_prior,dp))+1.d0))**2
    end if
    if (.not. positive .and. get_Cl_apod /= 0.d0) get_Cl_apod = 1.d0 / get_Cl_apod
  end function get_Cl_apod

  subroutine read_Cl_file(clfile, Dl, lmin, lmax, col_order)
    implicit none
    character(len=*),                      intent(in)  :: clfile, col_order
    real(dp),         dimension(lmin:,1:), intent(out) :: Dl
    integer(i4b),                          intent(in)  :: lmin, lmax

    integer(i4b)       :: unit, l, n, nspec
    character(len=512) :: line
    real(dp)           :: Dl_TT, Dl_TE, Dl_EE, Dl_BB
    real(dp),          allocatable, dimension(:,:) :: clin
    character(len=80),              dimension(120) :: header
    n     = len(trim(clfile))
    nspec = size(Dl,2)
    unit  = getlun()

    Dl = 0.d0
    if (clfile(n-3:n) == 'fits') then
       ! Assume input file is in standard Healpix FITS format
       allocate(clin(0:lmax,4))
       call fits2cl(clfile, clin, lmax, 4, header)
       do l = lmin, lmax
          Dl(l,TT) = clin(l,1) * (l*(l+1)/(2.d0*pi))
          if (nspec == 6) then
             Dl(l,EE) = clin(l,2)
             Dl(l,BB) = clin(l,3)
             Dl(l,TE) = clin(l,4)
          end if
       end do
       deallocate(clin)
    else
       ! Assume input file is in ASCII format with {l, TT, TE, EE, BB} ordering, and
       ! normalization is C_l * l(l+1)/2pi
       open(unit,file=trim(clfile))
       do while (.true.)
          read(unit,'(a)',end=3) line
          line = trim(line)
          if (line(1:1) == '#' .or. trim(line) == '') cycle
          if (nspec == 1) then
             read(line,*) l, Dl_TT
             if (l >= lmin .and. l <= lmax) then
                Dl(l,TT) = Dl_TT 
             end if
          else
             !write(*,*) trim(line)
             if (trim(col_order) == 'TT_TE_EE_BB') then
                read(line,*) l, Dl_TT, Dl_TE, Dl_EE, Dl_BB
             else if (trim(col_order) == 'TT_EE_BB_TE') then
                read(line,*) l, Dl_TT, Dl_EE, Dl_BB, Dl_TE
             else
                write(*,*) 'Unsupported column ordering = ', trim(col_order)
                stop
             end if
             if (l >= lmin .and. l <= lmax) then
                Dl(l,TT) = Dl_TT 
                Dl(l,TE) = Dl_TE 
                Dl(l,EE) = Dl_EE 
                Dl(l,BB) = Dl_BB 
             end if
          end if
       end do
3      close(unit)
    end if

    if (lmin == 0) then
       if (Dl(0,TT) == 0.d0) &
            & write(*,*) 'Warning: Input spectrum file ' // trim(clfile) &
            & // ' has vanishing monopole'
    end if
  end subroutine read_Cl_file

  subroutine binCls(self)
    implicit none
    class(comm_Cl), intent(inout) :: self

    real(dp)     :: val, n
    integer(i4b) :: j, k, l

    do k = 1, self%nbin
       do j = 1, self%nspec
          if (self%stat(k,j) == '0') then
             self%Dl(self%bins(k,1):self%bins(k,2),j) = 0.d0
             cycle
          end if
          val = 0.d0
          n   = 0.d0
          do l = self%bins(k,1), self%bins(k,2)
             val = val + (2*l+1) * self%Dl(l,j)
             n   = n   + 2*l+1
          end do
          self%Dl(self%bins(k,1):self%bins(k,2),j) = val/n
       end do
    end do
    
  end subroutine binCls

  recursive subroutine binCl2(self, bin)
    implicit none
    class(comm_Cl), intent(inout) :: self
    type(Cl_bin),   intent(in)    :: bin

    real(dp)     :: val, n
    integer(i4b) :: j, k, l

    val = 0.d0
    n   = 0.d0
    do l = bin%lmin, bin%lmax
       val = val + (2*l+1) * self%Dl(l,bin%spec)
       n   = n   + 2*l+1
    end do
    self%Dl(bin%lmin:bin%lmax,bin%spec) = val/n
    !write(*,*) bin%spec, bin%lmax, bin%lmax, self%Dl(bin%lmin:bin%lmax,bin%spec)
    
    do k = 1, bin%nsub
       call self%binCl2(bin%sub(k))
    end do

  end subroutine binCl2

  subroutine binCls2(self)
    implicit none
    class(comm_Cl), intent(inout) :: self

    integer(i4b) :: k

    do k = 1, self%nbin2
       call self%binCl2(self%bins2(k))
    end do
    
  end subroutine binCls2


  subroutine sampleCls(self, map, handle, ok)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map
    type(planck_rng), intent(inout) :: handle
    logical(lgt),     intent(inout) :: ok

!    return

    select case (trim(self%type))
    case ('none')
       return
    case ('binned')
       !call sample_Cls_inverse_wishart(self, map, handle, ok)
       call sample_Cls_inverse_wishart2(self, map, handle, ok)
    case ('power_law')
       call sample_Cls_powlaw(self, map, ok)
    case ('power_law_gauss')
       call sample_Cls_powlaw_gauss(self, map, ok)
    case ('exp')
       call sample_Cls_powlaw(self, map, ok)
    end select

    if (ok) call self%updateS
    
  end subroutine sampleCls

  subroutine sample_Cls_inverse_wishart(self, map, handle, ok)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    class(comm_map),  intent(in)    :: map
    type(planck_rng), intent(inout) :: handle
    logical(lgt),     intent(inout) :: ok

    integer(i4b) :: bin, b, i, j, k, l, m, n, p, ind, b1, b2, col, ierr, n_attempt
    logical(lgt) :: posdef
    logical(lgt), allocatable, dimension(:,:) :: pattern
    integer(i4b), allocatable, dimension(:) :: i2p
    real(dp), allocatable, dimension(:)     :: W
    real(dp), allocatable, dimension(:,:)   :: y, y_t
    real(dp), allocatable, dimension(:,:)   :: C_b
    real(dp), allocatable, dimension(:,:)   :: sigma, s, sigma_l

    allocate(sigma_l(0:self%lmax,self%nspec))

    call map%getSigmaL(sigma_l)

    if (self%info%myid == 0) then
       allocate(sigma(self%nmaps,self%nmaps), pattern(self%nmaps,self%nmaps))

       do bin = 1, self%nbin
          if (.not. any(self%stat(bin,:) == 'S') .and. .not. any(self%stat(bin,:) == 'M')) cycle
          
          sigma = 0.d0
          n     = 0.d0
          do l = self%bins(bin,1), self%bins(bin,2)
             n   = n + 2*l+1
             k = 1
             do i = 1, self%nmaps
                do j = i, self%nmaps
                   sigma(i,j) = sigma(i,j) + (2*l+1) * sigma_l(l,k) * l*(l+1)/(2.d0*pi)
                   k          = k+1
                end do
             end do
          end do

          do i = 1, self%nmaps
             do j = i, self%nmaps
                sigma(j,i) = sigma(i,j)
             end do
          end do

          ! Set up sampling pattern
          k = 1
          do i = 1, self%nmaps
             do j = i, self%nmaps
                pattern(i,j) = self%stat(bin,k) == 'S' .or. self%stat(bin,k) == 'M'
                pattern(j,i) = pattern(i,j)
                k            = k+1
             end do
          end do

          ! Keep sampling until all elements have been treated
          do while (any(pattern))
             ! Find which elements to include this time
             do col = 1, self%nmaps
                if (any(pattern(:,col))) exit
             end do

             ! Extract the appropriate segment
             p = count(pattern(:,col))
             allocate(s(p,p), y(p,1), y_t(1,p), i2p(p), C_b(p,p), W(p))
             j = 1
             do i = 1, self%nmaps
                if (pattern(i,col)) then
                   i2p(j) = i
                   j      = j+1
                end if
             end do
             s = sigma(i2p,i2p)
             if (all(s == 0.d0)) then
                pattern(i2p,i2p) = .false.
                deallocate(s, y, y_t, i2p, C_b, W)
                cycle
             end if

             call invert_matrix(s)
             call cholesky_decompose_single(s, ierr=ierr)
             if (ierr /= 0) then
                s = sigma(i2p,i2p)
                write(*,*) 'Error: Failed to cholesky decomp s, bin = ', bin, s, ierr
                deallocate(s, y, y_t, i2p, C_b, W)
                ok = .false.
                exit
             end if
             
             ! Draw sample
             C_b       = 0.d0
             posdef    = .false.
             n_attempt = 0
             do while (.not. posdef)
                do i = 1, n - p - 1
                   do j = 1, p
                      y(j,1) = rand_gauss(handle)
                   end do
                   y_t(1,:) = matmul(s, y(:,1))
                   y(:,1)   = y_t(1,:)
                   C_b      = C_b + matmul(y(:,1:1), y_t(1:1,:))
                end do
                call get_eigenvalues(C_b, W)
                posdef    = all(W > 1d-12)
                n_attempt = n_attempt+1
                if (n_attempt > 100) then
                   write(*,*) 'Error: Failed to sample positive definite C_b matrix in 100 attempts, bin = ', bin
                   ok = .false.
                   exit
                end if
             end do
             if (.not. ok) then
                deallocate(s, y, y_t, i2p, C_b, W)
                exit
             end if
             call invert_matrix(C_b)

             ! Copy information over to output Dl array
             do i = 1, p
                do j = i, p
                   ind = i2p(i)*(1-i2p(i))/2 + (i2p(i)-1)*self%nmaps + i2p(j)
                   self%Dl(self%bins(bin,1):self%bins(bin,2),ind) = C_b(i,j) 
                   !write(*,*) bin, i, j, C_b(i,j), sigma_l(self%bins(bin,1),ind) * self%bins(bin,1)*(self%bins(bin,1)+1)/2/pi
                end do
             end do
          
             ! Remove current elements from pattern, and prepare for next round
             pattern(i2p,i2p) = .false.
             deallocate(s, y, y_t, i2p, C_b, W)
             
          end do
          if (.not. ok) exit
       end do
       deallocate(sigma, pattern)
    end if

    call mpi_bcast(ok,      1,             MPI_LOGICAL,          0, self%info%comm, ierr)
    call mpi_bcast(self%Dl, size(self%Dl), MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)

    deallocate(sigma_l)

  end subroutine sample_Cls_inverse_wishart

  subroutine sample_Cls_inverse_wishart2(self, map, handle, ok)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    class(comm_map),  intent(in)    :: map
    type(planck_rng), intent(inout) :: handle
    logical(lgt),     intent(inout) :: ok

    integer(i4b) :: i, j, l, bin, ierr
    real(dp), allocatable, dimension(:)     :: ln_det_sigma_l
    real(dp), allocatable, dimension(:,:)   :: sigma, s
    real(dp), allocatable, dimension(:,:,:) :: sigma_l

    ! Local parameters for lnL
    integer(i4b) :: lmin, lmax, spec, p1, p2
    
    allocate(sigma_l(self%nmaps,self%nmaps,0:self%lmax))
    allocate(ln_det_sigma_l(0:self%lmax))

    call map%getSigmaL(sigma_l_mat=sigma_l)
    do l = 0, self%lmax
       do i = 1, self%nmaps
          if (sigma_l(i,i,l) == 0.d0) sigma_l(i,i,l) = 1.d0
       end do
       !write(*,*) l, sum(abs(sigma_l(:,:,l)))
       ln_det_sigma_l(l) = log_det(sigma_l(:,:,l))
    end do

    if (self%info%myid == 0) then
       if (self%lmin_lookup >= 0) then
          call sample_Dl_lookup(handle, ok)
       end if

       do i = 1, self%nbin2
          call sample_Dl_bin(self%bins2(i), handle, ok)
!          write(*,*) i, self%bins2(i)%lmin, self%Dl(self%bins2(i)%lmin,1)
!          write(*,*) "ok", ok
       end do
    end if

    call mpi_bcast(self%Dl, size(self%Dl), MPI_DOUBLE_PRECISION, 0, self%info%comm, ierr)
    !call mpi_allreduce(MPI_IN_PLACE, self%Dl, size(self%Dl), MPI_DOUBLE_PRECISION, &
    !     & MPI_SUM, self%info%comm, ierr)
    !call self%updateS()

!!$    if (self%info%myid == 0) then
!!$       do i = 1, self%nbin2
!!$          write(*,*) self%bins2(i)%lmin, self%bins2(i)%lmax, self%Dl(self%bins2(i)%lmin,1)
!!$       end do
!!$    end if
!!$    call mpi_finalize(ierr)
!!$    stop
    
    deallocate(sigma_l, ln_det_sigma_l)
    
  contains

    subroutine sample_Dl_lookup(handle, ok)
      implicit none
      type(planck_rng),                   intent(inout) :: handle
      logical(lgt),                       intent(inout) :: ok
    
      integer(i4b) :: i, j, k, l, m, n, status
      real(dp)     :: S(3,3), eta, w, ln_det_S
      real(dp), allocatable, dimension(:) :: lnL
      
      if (.not. ok) return

      ! Compute log-likelihood for each model
      n        = size(self%par_lookup)
      allocate(lnL(n))
      lnL = 0.d0
      do i = 1, n
         do l = self%lmin_lookup, self%lmax_lookup
            S = self%S_mat(:,:,l)
            m = 0
            do j = 1, self%nmaps
               do k = j, self%nmaps
                  m = m+1
                  if (self%active_lookup(m)) then
                     S(j,k) = self%Dl_lookup(l,m,i) / &
                          & (l*(l+1)/2.d0/pi*self%RJ2unit(j)*self%RJ2unit(k))
                     S(k,j) = S(j,k)
                  end if
               end do
            end do
            do j = 1, self%nmaps
               if (S(j,j) == 0.d0) S(j,j) = 1.d0
            end do
            call invert_matrix(S, cholesky=.true., status=status, ln_det=ln_det_S)
            if (status /= 0) then
               lnL(i) = -1.d30
            else
!!$               write(*,*) 
!!$               write(*,*) i, l, real(S,sp)
!!$               write(*,*) i, l, real(sigma_l(:,:,l),sp)
!!$               write(*,*) i, l, ln_det_S, trace_sigma_inv_Cl(sigma_l(:,:,l), S)
               lnL(i) = lnL(i) - 0.5d0 * (real(2*l+1,dp) * ln_det_S + &
                    & real(2*l+1,dp)*trace_sigma_inv_Cl(sigma_l(:,:,l), S))
            end if
         end do
      end do

      if (all(lnL == -1d30)) then
         ok = .false.
         return
      end if

!      write(*,*) 'lnL = ', lnL 

      ! Compute probabilities for each model spectrum
      where (lnL > -1d30)
         lnL = exp(lnL - maxval(lnL))
      elsewhere
         lnL = 0.d0
      end where
      lnL = lnL / sum(lnL)

      !write(*,*) 'P = ', lnL 

      ! Draw random model given respective probability
      eta = rand_uni(handle)
      w   = 0.d0
      i   = 0
      do while (w < eta)
         w = w + lnL(i+1)  ! lnL is now probability
         i = i + 1
      end do

      ! Update
      do l = self%lmin_lookup, self%lmax_lookup
         where (self%active_lookup)
            self%Dl(l,:) = self%Dl_lookup(l,:,i)
         end where
      end do
      write(*,*) '   Current lookup sample parameter = ', self%par_lookup(i), lnL(i)

      deallocate(lnL)

    end subroutine sample_Dl_lookup

    recursive subroutine sample_Dl_bin(bin, handle, ok)
      implicit none
      type(Cl_bin),                       intent(in)    :: bin
      type(planck_rng),                   intent(inout) :: handle
      logical(lgt),                       intent(inout) :: ok
    
      integer(i4b) :: i, j, k, l, status
      real(dp)     :: Dl_prop, Dl_in(3)
      real(dp)     :: prior(2)
      
      !write(*,*) bin%lmin, bin%lmax, bin%spec

      if (.not. ok) return
      if (bin%stat == 'S') then
         ! Set up prior = positive definite Dl; thus is not exact for TEB
         if (self%nspec == 1) then
            prior    = [0.d0, 1.d5]
         else
            prior    = [-1.d5, 1.d5]
            do l = bin%lmin, bin%lmax 
               write(*,*) "l:", l
               write(*,*) "lmin:", bin%lmin, "lmax:", bin%lmax 
               !write(*,*) "bin%spec", bin%spec
               select case (bin%spec)
               case (1)
                  prior(1) = max(prior(1), self%Dl(l,2)**2/self%Dl(l,4))
               case (2)
                  prior(1) = max(prior(1), -sqrt(self%Dl(l,1)*self%Dl(l,4)))
                  prior(2) = min(prior(2),  sqrt(self%Dl(l,1)*self%Dl(l,4)))
               case (3)
                  prior(1) = max(prior(1), -sqrt(self%Dl(l,1)*self%Dl(l,6)))
                  prior(2) = min(prior(2),  sqrt(self%Dl(l,1)*self%Dl(l,6)))
               case (4)
                  if (self%Dl(l,1) /= 0) then
                     prior(1) = max(prior(1), self%Dl(l,2)**2/self%Dl(l,1))
!                  else 
!                     prior(1) = 0.
                  end if
               case (5)
                  prior(1) = max(prior(1), -sqrt(self%Dl(l,4)*self%Dl(l,6)))
                  prior(2) = min(prior(2),  sqrt(self%Dl(l,4)*self%Dl(l,6)))
               case (6)
                  !write(*,*) "case(6) prior(1)", prior(1)
                  prior(1) = max(prior(1), 0.d0)
               end select
            end do
         end if
         Dl_in(2) = self%Dl(bin%lmin,bin%spec)
         Dl_in(1) = max(Dl_in(2) - 3*bin%sigma, 0.5d0*(Dl_in(2)+prior(1)))
         Dl_in(3) = min(Dl_in(2) + 3*bin%sigma, 0.5d0*(Dl_in(2)+prior(2)))

         ! Draw sample
         lmin=bin%lmin; lmax=bin%lmax; spec=bin%spec; p1=bin%p1; p2=bin%p2
         !write(*,*) "prior(1)=", real(prior(1), sp), "prior(2)=", real(prior(2),sp)
         !write(*,*) "Dl:", real(DL_in, sp)
         Dl_prop = sample_InvSamp(handle, Dl_in, lnL_invWishart, prior, status)
      
         ! Update
         !write(*,*) lmin, lmax, Dl_prop, status
         !stop
         if (status == 0) then
            self%Dl(bin%lmin:bin%lmax,bin%spec) = Dl_prop
         else
            ok = .false.
            write(*,*) lmin, lmax
         end if
          write(*,*) lmin, lmax, ok
      end if

      ! Sample all sub-bins
      do i = 1, bin%nsub
         call sample_Dl_bin(bin%sub(i), handle, ok)
      end do

    end subroutine sample_Dl_bin

    function lnL_invWishart(x)
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_invWishart
      
      integer(i4b) :: i, l, n, status
      real(dp)     :: ln_det_S
      real(dp), allocatable, dimension(:,:)   :: S
      
      allocate(S(self%nmaps,self%nmaps))
      
      lnL_invWishart = 0.d0
      do l = lmin, lmax
         S        = self%S_mat(:,:,l)
         S(p1,p2) = x / (l*(l+1)/2.d0/pi*self%RJ2unit(p1)*self%RJ2unit(p2))
         S(p2,p1) = S(p1,p2)
         do i = 1, self%nmaps
            if (S(i,i) == 0.d0) S(i,i) = 1.d0
         end do
         n = count(S(:,p1)/= 0.d0)
         call invert_matrix(S, cholesky=.true., status=status, ln_det=ln_det_S)
         if (status /= 0) then
            deallocate(S)
            lnL_invWishart = -1.d30
            return
         end if
         lnL_invWishart = lnL_invWishart - 0.5d0 * (real(2*l+1,dp) * ln_det_S + &
              & real(2*l+1,dp)*trace_sigma_inv_Cl(sigma_l(:,:,l), S))
      end do

      !write(*,*) x, lnL_invWishart
      
      deallocate(S)

    end function lnL_invWishart
    
  end subroutine sample_Cls_inverse_wishart2

  subroutine sample_Cls_powlaw(self, map, ok)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map
    logical(lgt),     intent(inout) :: ok    
    
  end subroutine sample_Cls_powlaw

  subroutine sample_Cls_powlaw_gauss(self, map, ok)
    implicit none
    class(comm_Cl),  intent(inout) :: self
    class(comm_map), intent(in)    :: map
    logical(lgt),     intent(inout) :: ok    
    
  end subroutine sample_Cls_powlaw_gauss

  subroutine writeFITS(self, chain, iter, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    integer(i4b),     intent(in)    :: chain, iter
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath

    character(len=4) :: ctext
    character(len=6) :: itext

    if (trim(self%type) == 'none') return
    if (self%info%myid /= 0) return

    call int2string(chain, ctext)
    call int2string(iter,  itext)
    
    select case (trim(self%type))
    case ('none')
       return
    case ('binned')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext, hdffile=hdffile, hdfpath=hdfpath)
    case ('power_law')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext, hdffile=hdffile, hdfpath=hdfpath)
    case ('power_law_gauss')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext, hdffile=hdffile, hdfpath=hdfpath)
    case ('exp')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext, hdffile=hdffile, hdfpath=hdfpath)
    case ('gauss')
       call write_Dl_to_FITS(self, 'c'//ctext//'_k'//itext)
       call write_powlaw_to_FITS(self, 'c'//ctext, hdffile=hdffile, hdfpath=hdfpath)
    end select

  end subroutine writeFITS

  subroutine write_Dl_to_FITS(self, postfix, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(in) :: self
    character(len=*), intent(in) :: postfix
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath
    
    integer(i4b) :: i, l, unit
    character(len=512) :: filename
    real(dp), allocatable, dimension(:,:) :: sigmal
    real(dp), allocatable, dimension(:) :: out

    if (trim(self%outdir) == 'none') return

    allocate(out(self%nspec))
    unit = getlun()
    filename = trim(self%outdir) // '/cls_' // trim(self%label) // '_' // trim(postfix) // '.dat'
    open(unit,file=trim(filename), recl=1024)
    if (self%nspec == 1) then
       write(unit,*) '# Columns are {l, Dl_TT}'
    else
       write(unit,*) '# Columns are {l, Dl_TT, Dl_TE, Dl_TB, Dl_EE, Dl_EB, Dl_BB}'
    end if
    do l = 0, self%lmax
       out = self%Dl(l,:)
       if (self%nspec == 1) then
          write(unit,fmt='(i6,e16.8)') l, out
       else
          write(unit,fmt='(i6,6e16.8)') l, out
       end if
    end do
    close(unit)
    deallocate(out)

    if (present(hdffile)) call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl', self%Dl)       
    
  end subroutine write_Dl_to_FITS

  subroutine write_powlaw_to_FITS(self, postfix, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    character(len=*), intent(in)    :: postfix
    type(hdf_file),   intent(in), optional :: hdffile
    character(len=*), intent(in), optional :: hdfpath

    integer(i4b) :: unit
    character(len=512) :: filename

    if (trim(self%outdir) == 'none') return

    unit    = getlun()
    filename = trim(self%outdir) // '/cls_' // trim(self%label) // '_powlaw_' // &
         & trim(postfix) // '.dat' 

    if (self%iter == 0) then
       open(unit,file=trim(filename), recl=1024)
       if (self%nspec == 1) then
          if (trim(self%type)=='power_law_gauss') then
             write(unit,*) '# Columns are {iteration, A_TT, beta_TT, theta_TT}'
          else
             write(unit,*) '# Columns are {iteration, A_TT, beta_TT}'
          end if
       else
          if (trim(self%type)=='power_law_gauss') then
             write(unit,*) '# Columns are {iteration, A_TT, A_EE, A_BB, beta_TT, beta_EE, beta_BB, theta_TT, theta_EE, theta_BB}'
          else
             write(unit,*) '# Columns are {iteration, A_TT, A_EE, A_BB, beta_TT, beta_EE, beta_BB}'
          end if
       end if
       close(unit)
    end if

    self%iter = self%iter + 1
    open(unit,file=trim(filename), recl=1024, position='append')
    if (self%nspec == 1) then
       if (trim(self%type)=='power_law_gauss') then
          write(unit,fmt='(i8,3e16.8)') self%iter, self%amp, self%beta, self%theta
       else
          write(unit,fmt='(i8,2e16.8)') self%iter, self%amp, self%beta
       end if
    else
       if (trim(self%type)=='power_law_gauss') then
          write(unit,fmt='(i8,9e16.8)') self%iter, self%amp, self%beta, self%theta
       else
          write(unit,fmt='(i8,6e16.8)') self%iter, self%amp, self%beta
       end if
    end if
    close(unit)

    if (present(hdffile)) then
       call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_amp',  self%amp)
       call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_beta', self%beta)
       if (trim(self%type)=='power_law_gauss') call write_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_theta', self%theta)
    end if
    
  end subroutine write_powlaw_to_FITS


  subroutine initHDF(self, hdffile, hdfpath)
    implicit none
    class(comm_Cl),   intent(inout) :: self
    type(hdf_file),   intent(in)    :: hdffile
    character(len=*), intent(in)    :: hdfpath

    if (trim(self%type) == 'none') return
    select case (trim(self%type))
    case ('none')
       return
    case ('binned')
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl', self%Dl)
    case ('power_law')
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_amp',  self%amp)
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_beta', self%beta)
       call self%updatePowLaw()
    case ('power_law_gauss')
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_amp',  self%amp)
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_beta', self%beta)
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_theta', self%theta)
       call self%updatePowLawGauss()
    case ('exp')
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_amp',  self%amp)
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_beta', self%beta)
       call self%updateExponential()
    case ('gauss')
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_amp',  self%amp)
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'/Dl_beta', self%beta)
       call self%updateGaussian()
    end select
    call self%updateS()

  end subroutine initHDF


  subroutine write_sigma_l(filename, sigma_l)
    implicit none
    character(len=*),                   intent(in) :: filename
    real(dp),         dimension(0:,1:), intent(in) :: sigma_l

    integer(i4b) :: unit, l, lmax, nspec

    lmax  = size(sigma_l,1)-1
    nspec = size(sigma_l,2)

    unit = getlun()
    open(unit,file=trim(filename), recl=1024)
    if (nspec == 1) then
       write(unit,*) '# Columns are {l, Dl_TT}'
    else
       write(unit,*) '# Columns are {l, Dl_TT, Dl_TE, Dl_TB, Dl_EE, Dl_EB, Dl_BB}'
    end if
    do l = 0, lmax
       if (nspec == 1) then
          write(unit,fmt='(i6,e16.8)') l, sigma_l(l,:) * l*(l+1)/(2.d0*pi)
       else
          write(unit,fmt='(i6,6e16.8)') l, sigma_l(l,:) * l*(l+1)/(2.d0*pi)
       end if
    end do
    close(unit)

  end subroutine write_sigma_l

  function getCl(self, l, p)
    implicit none
    class(comm_Cl),  intent(in) :: self    
    integer(i4b),    intent(in) :: l, p
    real(dp)                    :: getCl

    integer(i4b) :: j

    j = p*(1-p)/2 + (p-1)*self%nmaps + p
    if (l == 0) then
       getCl = self%Dl(l,j)
    else
       getCl = self%Dl(l,j) / (l*(l+1)/2.d0/pi)
    end if
    getCl = getCl * get_Cl_apod(l, self%l_apod, self%lmax, self%lmax_prior, .true.)**2

  end function getCl

  recursive subroutine set_Dl_bin(self, bin, Dl, pos, put)
    implicit none
    class(comm_Cl),                intent(inout) :: self
    type(Cl_bin),                  intent(in)    :: bin
    real(dp),       dimension(1:), intent(inout) :: Dl
    integer(i4b),                  intent(inout) :: pos
    logical(lgt),                  intent(in)    :: put
    
    integer(i4b) :: i

    if (bin%stat /= 'M') return
    pos     = pos+1
!    write(*,*) pos, bin%lmin, bin%spec, Dl(pos)
    if (put) then
       self%Dl(bin%lmin:bin%lmax,bin%spec) = Dl(pos)
    else
       Dl(pos) = self%Dl(bin%lmin,bin%spec)
    end if

    do i = 1, bin%nsub
       call self%set_Dl_bin(bin%sub(i), Dl, pos, put)
    end do
    
  end subroutine set_Dl_bin

  function check_posdef(self, bin, Dl)
    implicit none
    class(comm_Cl),                intent(in) :: self
    integer(i4b),                  intent(in) :: bin
    real(dp),       dimension(1:), intent(in) :: Dl
    logical(lgt)                              :: check_posdef

    integer(i4b) :: i, l, pos
    real(dp), allocatable, dimension(:)     :: w
    real(dp), allocatable, dimension(:,:,:) :: S

    allocate(S(self%nmaps, self%nmaps, self%bins2(bin)%lmin:self%bins2(bin)%lmax))
    allocate(w(self%nmaps))

    do l = self%bins2(bin)%lmin, self%bins2(bin)%lmax
       S(:,:,l) = self%S_mat(:,:,l)
    end do

    pos = 0
    call update_S(self%bins2(bin), pos, Dl, self%bins2(bin)%lmin, S)
    
    check_posdef = .true.
    do l = self%bins2(bin)%lmin, self%bins2(bin)%lmax
       do i = 1, self%nmaps
          if(S(i,i,l) == 0.d0) S(i,i,l) = 1.d0
       end do
!!$       write(*,*) l
!!$       do i = 1, self%nmaps
!!$          write(*,*) real(S(i,:,l),sp)
!!$       end do
       call get_eigenvalues(S(:,:,l), w)
!!$       write(*,*) any(w<0.d0)
!!$       write(*,*) 
       if (any(w < 0.d0)) then
          check_posdef = .false.
          exit
       end if
    end do

    deallocate(S, w)

  contains

    recursive subroutine update_S(bin, pos, Dl, lmin, S)
      implicit none
      type(Cl_bin),                         intent(in)    :: bin
      integer(i4b),                         intent(inout) :: pos
      real(dp),     dimension(1:),          intent(in)    :: Dl
      integer(i4b),                         intent(in)    :: lmin
      real(dp),     dimension(1:,1:,lmin:), intent(inout) :: S

      integer(i4b) :: l, i, j

      select case (bin%spec)
      case (1)
         i = 1; j = 1
      case (2)
         i = 1; j = 2
      case (3)
         i = 1; j = 3
      case (4)
         i = 2; j = 2
      case (5)
         i = 2; j = 3
      case (6)
         i = 3; j = 3
      case default
         write(*,*) 'Unsupported spectrum id'
         stop
      end select

      if (bin%stat == 'M') then
         pos    = pos+1
         do l = bin%lmin, bin%lmax
            S(i,j,l) = Dl(pos) / (l*(l+1)/(2.d0*pi)) / (self%RJ2unit(i)*self%RJ2unit(j))
            S(j,i,l) = S(i,j,l)
         end do
      end if

      do i = 1, bin%nsub
         call update_S(bin%sub(i), pos, Dl, lmin, S)
      end do

    end subroutine update_S

  end function check_posdef
    


end module comm_Cl_mod
