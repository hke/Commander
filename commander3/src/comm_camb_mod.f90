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
module comm_camb_mod
  use CAMB
  use params
  use comm_utils
  implicit none

  private
  public comm_camb, comm_camb_sample

  type comm_camb_sample
     real(dp), dimension(6)                :: theta    ! Cosmological parameters
     real(dp), dimension(:,:), allocatable :: c_l      ! Angular power spectrum
     real(dp), dimension(:,:), allocatable :: s_lm     ! Mean-field alms
     real(dp), dimension(:,:), allocatable :: f_lm     ! Fluctuation alms
   contains
     procedure :: dealloc => dealloc_camb_sample
  end type comm_camb_sample

  type comm_camb
     type(comm_camb_sample) :: theta_curr   ! Current sample 

     integer(i4b)     :: l_min              ! Minimum multipole
     integer(i4b)     :: l_max              ! Maximum multipole
     integer(i4b)     :: nalm               ! Number of alms = (l_lmax+1)**2
     integer(i4b)     :: nmaps              ! Number of polarization maps (typically {T,E,B})
     integer(i4b)     :: nspec              ! Number of power spectra (typically six, {TT,TE,TB,EE,EB,BB})
     integer(i4b)     :: nr_of_samples      
     character(len=4) :: dat_length
     real(dp)         :: proposal_multiplier
     real(dp),         dimension(2)   :: noise_l
     real(dp),         dimension(6)   :: correct_cosmo_param ! [Ombh2, omch2, H0, tau, ln(10^10As), n_s]
     real(dp),         dimension(6)   :: sigma_cosmo_param   ! [Ombh2, omch2, H0, tau, ln(10^10As), n_s]
     real(dp),         dimension(6,6) :: covariance_matrix   ! Parameter covariance matrix
     real(dp),         dimension(6,6) :: L_mat               ! Cholesky factor
     character(len=2), dimension(3)   :: spectra_list
   contains
     procedure sample => sample_camb_params
  end type comm_camb

  interface comm_camb
     procedure constructor
  end interface comm_camb

  interface comm_camb_sample
     procedure constructor_camb_sample
  end interface comm_camb_sample

contains

  function constructor(lmin, lmax)
    ! 
    ! Constructor for CAMB object
    ! 
    ! Arguments
    ! --------- 
    ! lmin:    int (scalar)
    !          Minimum multipole to be used for parameter fitting
    ! lmax:    int (scalar)
    !          Maximum multipole to be used for parameter fitting
    !
    ! Returns
    ! -------
    ! constructor: 
    !          pointer to CAMB object
    ! 
    implicit none
    integer(i4b),              intent(in) :: lmin, lmax
    class(comm_camb), pointer             :: constructor

    allocate(constructor)
    constructor%l_min = lmin
    constructor%l_max = lmax
    constructor%nmaps = 3
    constructor%nspec = constructor%nmaps*(constructor%nmaps+1)/2
    constructor%nalm  = (lmax+1)**2
    allocate(constructor%m_old%c_l(0:constructor%lmax,constructor%nspec))
    allocate(constructor%m_old%s_lm(0:constructor%nalm,constructor%nmaps))
    allocate(constructor%m_old%f_lm(0:constructor%nalm,constructor%nmaps))
    allocate(constructor%m_new%c_l(0:constructor%lmax,constructor%nspec))
    allocate(constructor%m_new%s_lm(0:constructor%nalm,constructor%nmaps))
    allocate(constructor%m_new%f_lm(0:constructor%nalm,constructor%nmaps))

    ! Static variables
    constructor%nr_of_samples       = 10
    constructor%dat_length          = '1501'
    constructor%proposal_multiplier = 0.3d0
    constructor%noise_l             = [3.04d-16, 6.08d-16]
    constructor%correct_cosmo_param = [0.02202d0, 0.1153d0, 68.2d0, 0.066d0, 3.035d0, 0.960d0]
    constructor%%sigma_cosmo_param  = [0.001d0,   0.005d0,  1.d0,   0.005d0, 0.005d0, 0.005d0]
    constructor%spectra_list        = ['TT', 'EE', 'TE']

  end function constructor

  function constructor_camb_sample(lmin, lmax, nmaps, s_init)
    ! 
    ! Constructor for CAMB sample object
    ! 
    ! Arguments
    ! --------- 
    ! lmin:    int (scalar)
    !          Minimum multipole to be used for parameter fitting
    ! lmax:    int (scalar)
    !          Maximum multipole to be used for parameter fitting
    ! nmaps:   int (scalar)
    !          Number of polarization maps (typically 3; {T,E,B})
    ! s_init:  derived type (comm_camb_sample; optional)
    !          Existing object; if passed, the returned object will be set equal to this
    !
    ! Returns
    ! -------
    ! constructor: 
    !          pointer to CAMB sample object
    ! 
    implicit none
    integer(i4b),                      intent(in)           :: lmin, lmax, nmaps
    class(comm_camb_sample),           intent(in), optional :: s_init
    class(comm_camb_sample),   pointer                      :: constructor_camb_sample
  
    integer(i4b) :: nalm, nspec

    nspec = nmaps*(nmaps+1)/2
    nalm  = (lmax+1)**2

    allocate(constructor_camb_sample%c_l(0:lmax,nspec))
    allocate(constructor_camb_sample%s_lm(0:nalm,nmaps))    
    allocate(constructor_camb_sample%f_lm(0:nalm,nmaps))

    if (present(s_init)) then
       constructor_camb_sample%c_l  = s_init%c_l
       constructor_camb_sample%s_lm = s_init%s_lm
       constructor_camb_sample%f_lm = s_init%f_lm
    end if

  end function constructor_camb_sample

  subroutine dealloc_camb_sample(self)
    ! 
    ! Cleanup routine for comm_camb_sample object
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived type (comm_camb_sample)
    !          Object to be deallocated
    ! 
    implicit none
    class(comm_camb_sample), intent(inout) :: self
  
    if (allocated(constructor_camb_sample%c_l))  deallocate constructor_camb_sample%c_l
    if (allocated(constructor_camb_sample%s_lm)) deallocate constructor_camb_sample%s_lm
    if (allocated(constructor_camb_sample%f_lm)) deallocate constructor_camb_sample%f_lm

  end subroutine dealloc_camb_sample


  subroutine sample_camb_params(self, rng_handle)
    ! 
    ! Routine for drawing new samples with Metropolis MCMC and a joint (C_l, s) accept-reject ratio;
    ! see Racine et al. (2016) for details.
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_camb)
    !          CAMB object
    ! rng_handle: derived type (planck_rng)
    !          Random number handle
    ! 
    implicit none
    class(comm_camb), intent(inout) :: self
    type(planck_rng), intent(inout) :: rng_handle

    real(dp), dimension(:, :),    allocatable :: list_of_cosmo_param
    real(dp), dimension(:, :, :), allocatable :: list_of_sigma_l

    real(dp), dimension(:, :),    allocatable :: cur_sigma_l, hat_c_l, d_lm, correct_c_l

    integer(i4b) :: sample_nr, i, j, accepted_samples
    logical(lgt) :: accept
    real(dl), dimension(6) :: average, var
    class(comm_camb_sample), pointer :: old_sample, new_sample, correct_sample
    
    allocate(d_lm(2, (l_max+1)**2))
    allocate(list_of_cosmo_param(6, nr_of_samples))
    allocate(list_of_sigma_l(3, 0:l_max, nr_of_samples))
    allocate(cur_sigma_l(3, 0:l_max))
    allocate(hat_c_l(3, 0:l_max))
    allocate(old_sample%c_l(3, 0:l_max))
    allocate(new_sample%c_l(3, 0:l_max))
    list_of_sigma_l = 0._dl
    cur_sigma_l = 0._dl
    hat_c_l = 0._dl
    accepted_samples = 0

    ! Initialize
    correct_sample%theta = correct_cosmo_param
    call rand_init(rng_handle, 12345, 1245689)
    call get_c_l_from_camb(correct_sample)
    call init_CMB_and_noise(d_lm, correct_sample, rng_handle)
    call init_covariance_matrix(L_mat, covariance_matrix)
    
    ! First sample
    old_sample%theta = correct_cosmo_param
    call get_c_l_from_camb(old_sample)
    call get_s_lm_f_lm(old_sample, d_lm, rng_handle)
    
    DO sample_nr = 1, nr_of_samples
       ! Get new theta, get c_l, s_lm and f_lm from new theta. Then rescale f_lm
       call get_new_sample(accept, new_sample, old_sample, d_lm, L_mat, rng_handle)
       
       print *, 'Sample:', sample_nr, 'Out of', nr_of_samples
       print *, 'New Sample OmbH2', new_sample%theta(1)
       print *, 'Old Sample OmbH2', old_sample%theta(1)
       
       ! Save information about new sample that will be printed to dat files
       list_of_cosmo_param(:, sample_nr) = new_sample%theta
       call get_c_l_from_a_lm(cur_sigma_l, new_sample%s_lm + new_sample%f_lm)
       list_of_sigma_l(:, :, sample_nr) = cur_sigma_l
       
       if (accept) then
          ! Sample was accepted
          accepted_samples = accepted_samples + 1 
          old_sample = new_sample
       end if
    END DO

    ! Done sampling, save data
    print*, 'Accepted Samples: ', accepted_samples, 'Ratio:', real(accepted_samples)/real(nr_of_samples)
    average = sum(list_of_cosmo_param, dim = 2)/nr_of_samples
    print *, 'Average:', average
    var = 0._dl
    DO i = 1, nr_of_samples
       var = var + (list_of_cosmo_param(:, i) - average)**2
    END DO
    var = var / (nr_of_samples - 1)
    print *, 'Deviation Param:', sqrt(var)
    
    call get_c_l_from_a_lm(hat_c_l, d_lm)

    DO i = 1, 3
       open(unit=1, file='sigma_'//spectra_list(i)//'_l_out.dat', status='replace', action='write')
       open(unit=2, file='hat_'//spectra_list(i)//'_c_l_out.dat', status='replace', action='write')
       DO j = 1, nr_of_samples  
          write(1, '( '//dat_length//'(2X, ES14.6) )') list_of_sigma_l(i, :, j)
       END DO
       write(2, '( '//dat_length//'(2X, ES14.6) )') hat_c_l(i, :)
       close(1)
       close(2)
    END DO

    open(unit=1, file='cosmo_param_out.dat', status='replace', action='write')
    do i = 1, nr_of_samples
       write(1, '( 6(2X, ES14.6) )') list_of_cosmo_param(:, i)
    end do
    close(1)

  end subroutine sample_camb_params


  subroutine get_new_sample(accept, new_sample, old_sample, d_lm, L_mat, rng_handle)
    ! 
    ! Gets new theta sample. Finds c_l from theta, and then s_lm and scaled f_lm
    ! 
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    ! 
    implicit none
    logical(lgt),                                 intent(out) :: accept
    type(sample),                                 intent(out) :: new_sample
    type(sample),                                 intent(in)  :: old_sample
    real(dp),         dimension(2, (l_max+1)**2), intent(in)  :: d_lm
    real(dp),         dimension(6, 6),            intent(in)  :: L_mat
    type(planck_rng),                             intent(in)  :: rng_handle
     
    real(dp), dimension(2, (l_max+1)**2) :: scaled_f_lm
    
    call cosmo_param_proposal(new_sample, old_sample, L_mat, rng_handle)
    call get_c_l_from_camb(new_sample)
    call get_s_lm_f_lm(new_sample, d_lm, rng_handle)
    call get_scaled_f_lm(scaled_f_lm, new_sample, old_sample) 
    accept = acceptance(scaled_f_lm, new_sample, old_sample, d_lm, rng_handle)
  end subroutine get_new_sample

  subroutine init_covariance_matrix(L, covariance_matrix)
    ! 
    ! Caclulates a multivariate Gaussian for the proposal function w
    ! If there exists a cosmo_param_out.dat file in the parent folder it uses
    ! that to calculate a covariance matrix and then uses Cholesky decomposition
    ! to find the L (L*L^T=Covariance matrix), and so theta_new = theta_old + L*z
    ! where z is a random Gaussian vector
    ! 
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    ! 
    implicit none
    real(dp), dimension(6, 6), intent(out) :: L, covariance_matrix

    integer(i4b) :: i, j, k, nlines
    real(dp)     :: sum
    logical(lgt) :: previous_sample
    real(dp),    dimension(:, :), allocatable :: old_samples 
    real(dp),    dimension(6)                 :: averages
    
    INQUIRE(FILE="../cosmo_param_out.dat", EXIST=previous_sample)

    if (previous_sample) then
       ! Caclulate covariance matrix
       print *, 'Found previous sample chain. Calculating covariance matrix and use Cholesky decomposition in proposals for this run'
       open(1, file="../cosmo_param_out.dat", status='old', action='read')
       
       nlines = 0
       DO
          READ (1, *, END=10)
          nlines = nlines + 1 
       END DO
10     close(1)
       
       allocate(old_samples(nlines, 6))
       
       open(1, file="../cosmo_param_out.dat", status='old', action='read')
       DO k = 1, nlines
          READ (1, *) old_samples(k, :)
       END DO
       close(1)
       
       averages = sum(old_samples, dim = 1) / nlines
       print *, averages
       DO j = 1,6
          DO k = 1,6
             covariance_matrix(j, k) = 0._dl
             DO i = 1, nlines
                covariance_matrix(j, k) = covariance_matrix(j, k) + (old_samples(i, j) - averages(j))*(old_samples(i, k) - averages(k))
             END DO
          END DO
       END DO
       covariance_matrix = covariance_matrix / (nlines - 1)
    else
       print *, 'Did not find previous sample chain. Using hard-coded diagonal covariance matrix.'
       covariance_matrix = 0._dl
       DO i = 1, 6
          covariance_matrix(i, i) = sigma_cosmo_param(i)**2
       END DO
    end if
    
    L = 0._dl
    DO i = 1, 6
       DO j = 1, i
          sum = 0._dl
          
          if (i == j) then
             DO k = 1, j
                sum = sum + L(j, k)**2
             END DO
             L(j, j) = sqrt(covariance_matrix(j, j) - sum)
             
          else
             DO k = 1, j
                sum = sum + L(i, k) * L(j, k)
             END DO
             L(i, j) = (covariance_matrix(i, j) - sum) / L(j, j)
          end if
       END DO
    END DO
  end subroutine init_covariance_matrix

  subroutine get_c_l_from_a_lm(c_l, a_lm)  
    ! 
    ! Calculates power spectra from a_lm. TT, EE, and TE
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    implicit none
    real(dp), dimension(3, 0:l_max) :: c_l
    real(dp), dimension(2, (l_max+1)**2), intent(in) :: a_lm
    
    integer(i4b) :: k, l, m, index
    real(dp)     :: cur_c_l

    c_l = 0.d0
    ! Do TT and EE
    DO k = 1, 2
       DO l = l_min, l_max
          cur_c_l = a_lm(k, l**2 + l + 1)**2
          DO m = 1, l
             index = l**2 + l + m + 1
             cur_c_l = cur_c_l + 2*a_lm(k, index)**2
          END DO
          c_l(k, l) = cur_c_l / (2*l+1)
       END DO
    END DO
    
    ! Do TE
    DO l = l_min, l_max
       cur_c_l = a_lm(1, l**2 + l + 1)*a_lm(2, l**2 + l + 1)
       DO m = 1, l
          index = l**2 + l + m + 1
          cur_c_l = cur_c_l + 2*a_lm(1, index)*a_lm(2, index)
       END DO
       c_l(3, l) = cur_c_l / (2*l+1)
    END DO

  end subroutine get_c_l_from_a_lm


  subroutine cosmo_param_proposal(new_sample, old_sample, L, rng_handle)
    ! 
    ! Proposal function w. Finds new sample based on covariance
    ! matrix L*L^T. Proposal_multiplier = 0.3 to make sure the proposal theta
    ! becomes too large.
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    implicit none
    type(comm_camb_sample),                  intent(out)   :: new_sample
    type(comm_camb_sample),                  intent(in)    :: old_sample
    real(dp),               dimension(6, 6), intent(in)    :: L
    type(planck_rng),                        intent(inout) :: rng_handle
    
    real(dp), dimension(6) :: z
    integer(i4b) :: i
    
    DO i = 1, 6
       z(i) = proposal_multiplier * rand_gauss(rng_handle)
    END DO
    
    new_sample%theta = old_sample%theta + mat_vec_mul(L, z)
    
  end subroutine cosmo_param_proposal


  function acceptance(scaled_f_lm, new_sample, old_sample, d_lm, rng_handle)
    ! 
    ! This function determines if the new sample should be accepted or not.
    ! Assumes no priors and that the proposal is symmetric. Hence 
    ! A = min(1, pi(theta^{i+1})/pi(theta^i))
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    implicit none
    real(dp),               dimension(2, (l_max+1)**2), intent(in)    :: scaled_f_lm, d_lm
    type(comm_camb_sample),                             intent(in)    :: old_sample, new_sample
    type(planck_rng),                                   intent(inout) :: rng_handle

    real(dp), dimension(2, (l_max+1)**2) :: old_s_lm, old_f_lm, new_s_lm
    real(dp), dimension(3, 0: l_max) :: old_c_l, new_c_l
    real(dp) :: ln_pi_ip1, ln_pi_i, probability, uni
    real(dp), dimension(2, 2) :: new_S, old_S
    integer(i4b) :: k, i, l, m
    logical(i4b) :: acceptance

    old_s_lm = old_sample%s_lm
    old_f_lm = old_sample%f_lm
    old_c_l = old_sample%c_l
    new_s_lm = new_sample%s_lm
    new_c_l = new_sample%c_l

    ln_pi_ip1 = 0.d0
    ln_pi_i = 0.d0
  
    DO l = l_min, l_max
       DO m = 0, l
          i = l**2 + l + m + 1
          
          ! This part is a bit ugly. Everything is diagonal except c_l (because of
          ! C^TE) and so that is done after the k loop. k=1 is a^T_lm and k=2 is
          ! a^E_lm
          DO k = 1, 2   
             ln_pi_ip1 = ln_pi_ip1 + (d_lm(k, i) - new_s_lm(k, i))**2 / noise_l(k) + scaled_f_lm(k, i)**2 / noise_l(k)
             ln_pi_i   = ln_pi_i   + (d_lm(k, i) - old_s_lm(k, i))**2 / noise_l(k) + old_f_lm(k, i)**2 / noise_l(k)
          END DO
          new_S     = reshape((/ new_c_l(1, l), new_c_l(3, l), new_c_l(3, l), new_c_l(2, l) /), shape(new_S))
          old_S     = reshape((/ old_c_l(1, l), old_c_l(3, l), old_c_l(3, l), old_c_l(2, l) /), shape(old_S))
          ln_pi_ip1 = ln_pi_ip1 + dot_product(new_s_lm(:, i), mat_vec_mul(inv_mat(new_S), new_s_lm(:, i)))
          ln_pi_i   = ln_pi_i + dot_product(old_s_lm(:, i), mat_vec_mul(inv_mat(old_S), old_s_lm(:, i)))
       END DO
    END DO
  
    probability = exp(-(ln_pi_ip1 - ln_pi_i) / 2.0_dl)
    print *, 'prob:', probability
  
    uni = rand_uni(rng_handle) 
    if (uni < probability) then
       acceptance = .true.
       print *, '---------- ACCEPTED -----------'
    else
       acceptance = .false.
       print *, '-------- NOT ACCEPTED ---------'
    end if
  end function acceptance

  subroutine init_CMB_and_noise(d_lm, cur_sample, rng_handle)
    ! 
    ! Initializes simulated d_lm which is a_lm from CMB plus noise.
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    ! 
    implicit none
    real(dp),               dimension(2, (l_max+1)**2), intent(out)   :: d_lm
    type(comm_camb_sample),                             intent(in)    :: cur_sample
    type(planck_rng),                                   intent(inout) :: rng_handle  
    
    integer(i4b) :: index, l, m, k
    real(dp), dimension(4) :: z
    real(dp), dimension(3, 0: l_max) :: c_l  
    
    ! Do a^T_lm and a^E_lm
    c_l = cur_sample%c_l
    d_lm = 0._dl
    DO l = l_min, l_max
       DO m = 0, l
          index = l**2 + l + m + 1
          z = (/ rand_gauss(rng_handle), rand_gauss(rng_handle), rand_gauss(rng_handle), rand_gauss(rng_handle) /)  
          
          d_lm(1, index) = sqrt(c_l(1, l)) * z(1) + sqrt(noise_l(1)) * z(2)
          d_lm(2, index) = c_l(3, l) / sqrt(c_l(1, l)) * z(1) + sqrt(c_l(2, l) - c_l(3, l)**2 / c_l(1, l)) * z(3) + sqrt(noise_l(2)) * z(4)
       END DO
    END DO
    DO i = 1, 3
       open(unit=1, file='lcdm_c_'//spectra_list(i)//'_l_out.dat', status='replace', action='write')
       write(1, '( '//dat_length//'(2X, ES14.6) )') c_l(i, :)
       close(1)
    END DO

  end subroutine init_CMB_and_noise

  subroutine get_s_lm_f_lm(cur_sample, d_lm, rng_handle)
    ! 
    ! Calculates mean field s_lm and fluctuation f_lm from c_l from cur_sample
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    ! 
    implicit none
    type(comm_camb_sample),                             intent(inout) :: cur_sample
    real(dp),               dimension(2, (l_max+1)**2), intent(in)    :: d_lm
    type(planck_rng),                                   intent(inout) :: rng_handle
    
    integer(i4b) :: l, m, k, index
    real(dp), dimension(2, (l_max+1)**2) :: s_lm, f_lm
    real(dp), dimension(3, 0: l_max) :: c_l  
    real(dp), dimension(2, 2) :: common_matrix, S_mat, N, S_inv, N_inv, S_sqrt_inv, N_sqrt_inv
    real(dp), dimension(2) :: d_vector, omega_1, omega_2 

    c_l = cur_sample%c_l
    s_lm = 0._dl
    f_lm = 0._dl 
    
    DO l = l_min, l_max
       S_mat = reshape((/ c_l(1, l), c_l(3, l), c_l(3, l), c_l(2, l) /), shape(S_mat))
       N = reshape((/ noise_l(1), 0._dl, 0._dl, noise_l(2) /), shape(N))
       N_inv = inv_mat(N)
       S_inv = inv_mat(S_mat)
       S_sqrt_inv = sqrt_mat(S_inv)
       N_sqrt_inv = sqrt_mat(N_inv)
       common_matrix = inv_mat(S_inv + N_inv)
       DO m = 0, l
          index = l**2 + l + m + 1
          d_vector =  (/ d_lm(1, index), d_lm(2, index) /)
          omega_1 = (/ rand_gauss(rng_handle), rand_gauss(rng_handle) /)
          omega_2 = (/ rand_gauss(rng_handle), rand_gauss(rng_handle) /)
          s_lm(:, index) = mat_vec_mul(matmul(common_matrix, N_inv), d_vector)
          f_lm(:, index) = mat_vec_mul(matmul(common_matrix, S_sqrt_inv), omega_1) + mat_vec_mul(matmul(common_matrix, N_sqrt_inv), omega_2) 
       END DO
    END DO
    cur_sample%s_lm = s_lm
    cur_sample%f_lm = f_lm 
  end subroutine get_s_lm_f_lm
  
  subroutine get_scaled_f_lm(scaled_f_lm, new_sample, old_sample)
    ! 
    ! Scaled f_lm from new_sample. f_scaled = sqrt(c_l^{i+1}/c_l^i)f^{i+1}
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    implicit none
    
    type(comm_camb_sample), intent(in) :: new_sample, old_sample
    integer(i4b) :: index, l, m, k
    real(dp) :: prefactor
    real(dp), dimension(3, 0: l_max)     :: old_c_l, new_c_l
    real(dp), dimension(2, (l_max+1)**2) :: old_f_lm
    real(dp), dimension(2, (l_max+1)**2) :: scaled_f_lm
    real(dp), dimension(2, 2)            :: new_S, old_S, inv_old_S
    
    old_f_lm = old_sample%f_lm
    old_c_l = old_sample%c_l
    new_c_l = new_sample%c_l
    DO l = l_min, l_max
       !Scaling is non-trivial when C^TE != 0, then we need to do matrix operations
       new_S = reshape((/ new_c_l(1, l), new_c_l(3, l), new_c_l(3, l), new_c_l(2, l) /), shape(new_S))
       old_S = reshape((/ old_c_l(1, l), old_c_l(3, l), old_c_l(3, l), old_c_l(2, l) /), shape(old_S))
       inv_old_S = inv_mat(old_S)
       DO m = 0, l
          index = l**2 + l + m + 1
          scaled_f_lm(:, index) = mat_vec_mul( matmul(sqrt_mat(new_S)  ,  sqrt_mat(inv_old_S)), old_f_lm(:, index)) 
       END DO
    END DO
  end subroutine get_scaled_f_lm
  
  subroutine get_c_l_from_camb(cur_sample)
    ! 
    ! Gets TT, EE, and TE power spectra from camb using the cosmological
    ! parameters in theta.
    !
    ! Arguments
    ! ---------
    ! (COMPLETE THIS)
    ! 
    implicit none
  
    type(comm_camb_sample) :: cur_sample 
    real(dp), dimension(6) :: cosmo_param
    type(CAMBparams) P
    type(CAMBdata) camb_data
    integer(i4b) :: l, k
    real(dp), dimension(3, 0: l_max) :: c_l
    cosmo_param = cur_sample%theta
    call CAMB_SetDefParams(P)
    
    P%ombh2 = cosmo_param(1)
    P%omch2 = cosmo_param(2)
    P%omk = 0._dl
    P%H0= cosmo_param(3)
    !call P%set_H0_for_theta(0.0104_dl)!cosmo_param(4)
    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
       InitPower%As = exp(cosmo_param(5))*1e-10
       InitPower%ns = cosmo_param(6)
       InitPower%r = 0.0
    end select
    
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
       Reion%use_optical_depth = .true.
       Reion%optical_depth = cosmo_param(4)
    end select
    
    P%WantScalars = .true.
    P%WantTensors = .true.
    P%Accuracy%AccurateBB  = .true.
    
    P%Max_l=2500
    P%Max_eta_k=6000
    P%Max_l_tensor=2500
    P%Max_eta_k_tensor=6000
    
    call CAMB_GetResults(camb_data, P)
    
    ! Set TT, EE, and TE
    c_l = 0._dl
    DO k = 1, 3
       c_l(k, :) = camb_data%CLData%Cl_scalar(0:l_max, k)
       
       c_l(k, 0) = 0._dl
       c_l(k, 1) = 0._dl
       DO l = 2, l_max
          c_l(k, l) = 2._dl * 3.14159265_dl / (l * (l + 1)) * c_l(k, l)
       END DO
    END DO
    cur_sample%c_l = c_l
  end subroutine get_c_l_from_camb

  ! HKE: Replace these function calls with
  ! invert_matrix (from math_tools)
  ! compute_hermitian_root (from math_tools)
  ! matmul (native Fortran routine)

!!$  ! TOOL FUNCTIONS
!!$  function inv_mat(A) result(inv_A)
!!$    real(dl), dimension(2, 2), intent(in) :: A
!!$    real(dl), dimension(2, 2) :: inv_A
!!$    real(dl) :: det
!!$    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
!!$    
!!$    inv_A(1, 1) = A(2, 2)
!!$    inv_A(2, 2) = A(1, 1)
!!$    inv_A(1, 2) = -A(2, 1)
!!$    inv_A(2, 1) = -A(1, 2)
!!$    inv_A = inv_A / det
!!$  end function inv_mat
!!$  
!!$  function sqrt_mat(A) result(sqrt_A)
!!$    ! Assumes the matrix is 2x2 and symmetric
!!$    real(dl), dimension(2, 2), intent(in) :: A
!!$    real(dl), dimension(2, 2) :: sqrt_A
!!$    real(dl) :: s, t
!!$    s = sqrt(A(1,1)*A(2,2) - A(1,2)**2)
!!$    t = sqrt(A(1,1) + A(2,2) + 2*s)
!!$    
!!$    sqrt_A(1, 1) = A(1, 1) + s
!!$    sqrt_A(2, 2) = A(2, 2) + s
!!$    sqrt_A(1, 2) = A(1, 2)
!!$    sqrt_A(2, 1) = A(2, 1)
!!$    sqrt_A = sqrt_A / t
!!$  end function sqrt_mat
!!$  
!!$  function mat_vec_mul(A, v) result(v_out)  
!!$    real(dl), dimension(:, :), intent(in) :: A
!!$    real(dl), dimension(:), intent(in) :: v
!!$    
!!$    real(dl), dimension(:), allocatable :: v_out
!!$    integer :: k, i, v_dim
!!$    integer, dimension(1) :: v_shape
!!$    v_shape = shape(v)
!!$    v_dim = v_shape(1)
!!$    allocate(v_out(v_dim))
!!$    v_out = 0._dl
!!$    DO i = 1, v_dim
!!$       DO k = 1, v_dim
!!$          v_out(i) = v_out(i) + A(i, k) * v(k)
!!$       END DO
!!$    END DO
!!$  end function mat_vec_mul

end program comm_camb_mod

