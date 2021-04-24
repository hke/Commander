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
module comm_tod_gain_mod
  use comm_tod_mod
  use comm_tod_noise_mod
  use comm_utils
  use comm_fft_mod
  implicit none



contains

  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  !  subroutine sample_gain_per_scan(self, handle, det, scan_id, n_corr, mask, s_ref)
  !  subroutine calculate_gain_mean_std_per_scan(self, det, scan_id, s_tot, invn, mask)
   subroutine calculate_gain_mean_std_per_scan(tod, scan_id, s_invN, mask, s_ref, s_tot, handle, mask_lowres, tod_arr)
    implicit none
    class(comm_tod),                      intent(inout) :: tod
    real(sp),             dimension(:,:), intent(in)    :: s_invN, mask, s_ref, s_tot
    integer(i4b),                         intent(in)    :: scan_id
    type(planck_rng),                     intent(inout)  :: handle
    real(sp),             dimension(:,:), intent(in), optional :: mask_lowres
    real(sp),             dimension(:,:), intent(in), optional :: tod_arr


    real(sp), allocatable, dimension(:,:) :: residual
    real(sp), allocatable, dimension(:)   :: r_fill
    integer(i4b) :: ext(2), j, i, ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_tot,1)

!    allocate(residual(size(stot_invN)))
      !g_old = tod%scans(scan_id)%d(det)%gain 

!   residual = (tod%scans(scan_id)%d(det)%tod - (tod%gain0(0) + &
!         & tod%gain0(det)) * s_tot) * mask

    allocate(r_fill(size(s_tot,1)))
    call tod%downsample_tod(s_tot(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet))
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          residual(:,j) = 0.d0
          cycle
       end if
       if (present(tod_arr)) then
         r_fill = tod_arr(:, j) - tod%scans(scan_id)%d(j)%baseline & 
           & - (tod%gain0(0) + tod%gain0(j)) * s_tot(:,j)
       else
         r_fill = tod%scans(scan_id)%d(j)%tod - tod%scans(scan_id)%d(j)%baseline & 
           & - (tod%gain0(0) + tod%gain0(j)) * s_tot(:,j)
       end if
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', real(tod%scans(scan_id)%d(j)%N_psd%sigma0, sp), handle, tod%scans(scan_id)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do
    call multiply_inv_N(tod, scan_id, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    ! Get a proper invn multiplication from Haavard's routine here
    do j = 1, ndet
       if (.not. tod%scans(scan_id)%d(j)%accept) then
          tod%scans(scan_id)%d(j)%gain  = 0.d0
          tod%scans(scan_id)%d(j)%dgain = 0.d0
       else
          if (present(mask_lowres)) then
             tod%scans(scan_id)%d(j)%dgain         = sum(s_invN(:,j) * residual(:,j) * mask_lowres(:,j))
             tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invN(:,j) * s_ref(:,j)    * mask_lowres(:,j))
          else
             tod%scans(scan_id)%d(j)%dgain         = sum(s_invN(:,j) * residual(:,j))
             tod%scans(scan_id)%d(j)%gain_invsigma = sum(s_invN(:,j) * s_ref(:,j))
          end if
          if (tod%scans(scan_id)%d(j)%gain_invsigma < 0.d0) then
             write(*,*) 'Warning: Not positive definite invN = ', tod%scanid(scan_id), j, tod%scans(scan_id)%d(j)%gain_invsigma
          end if
       end if
    end do
!    tod%scans(scan_id)%d(det)%dgain      = sum(stot_invN * residual)
!    tod%scans(scan_id)%d(det)%gain_invsigma = sum(stot_invN * s_tot)
!    if (tod%scans(scan_id)%d(det)%gain_invsigma < 0) then
!       write(*,*) 's', sum(mask * stot_invN * s_tot), sum(stot_invN * s_tot)
!    end if

!!$    if
!(tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_invsigma > 5.) then
!!$       open(58, file='tod.dat')
!!$       do i = p, q
!!$          if (mask(i) > 0.5) write(58,*) i, residual(i), s_tot(i)*0.001 
!!$       end do
!!$       close(58)
!!$       stop
!!$    end if

    !write(*,*) det, scan_id, real(tod%scans(scan_id)%d(det)%dgain/tod%scans(scan_id)%d(det)%gain_invsigma,sp), real(tod%gain0(0),sp), real(tod%gain0(det),sp), real(g_old,sp)

   ! write(*,*) tod%scanid(scan_id), real(tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(3)%gain_invsigma,sp), real(tod%gain0(0) + tod%gain0(3) + tod%scans(scan_id)%d(3)%dgain/tod%scans(scan_id)%d(3)%gain_invsigma,sp), '# deltagain'

    !if (.false. .and. trim(tod%freq) == '030' .and. mod(tod%scanid(scan_id),100) == 0) then
    if (.false.) then
       call int2string(tod%scanid(scan_id), itext)
       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%scans(scan_id)%d(1)%dgain/tod%scans(scan_id)%d(1)%gain_invsigma
       open(58,file='gain_delta_'//itext//'.dat')
       do i = ext(1), ext(2)
          write(58,*) i, residual(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_ref(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_invN(i,1)
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          if (present(tod_arr)) then
            write(58,*) i, tod_arr(i, 1) - tod%scans(scan_id)%d(1)%baseline &
            & - (tod%gain0(0) +  tod%gain0(1)) * s_tot(i,1)
          else
            write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - tod%scans(scan_id)%d(1)%baseline &
            & - (tod%gain0(0) +  tod%gain0(1)) * s_tot(i,1)
          end if
       end do
       write(58,*)
       do i = 1, size(s_tot,1)
          if (present(tod_arr)) then
            write(58,*) i, tod_arr(i, 1) - tod%scans(scan_id)%d(1)%baseline
          else
            write(58,*) i, tod%scans(scan_id)%d(1)%tod(i) - tod%scans(scan_id)%d(1)%baseline
          end if
       end do
       close(58)
    end if


    deallocate(residual, r_fill)

  end subroutine calculate_gain_mean_std_per_scan


  ! Compute gain as g = (d-n_corr-n_temp)/(map + dipole_orb), where map contains an 
  ! estimate of the stationary sky
  ! Eirik: Update this routine to sample time-dependent gains properly; results should be stored in self%scans(i)%d(j)%gain, with gain0(0) and gain0(i) included
  subroutine sample_smooth_gain(tod, handle, dipole_mods)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),   dimension(:, :),       intent(in)     :: dipole_mods


!    real(sp), dimension(:, :) :: inv_gain_covar ! To be replaced by proper matrix, and to be input as an argument
    integer(i4b) :: i, j, k, ndet, nscan_tot, ierr, ind(1)
    integer(i4b) :: currstart, currend, window, i1, i2, pid_id, range_end
    real(dp)     :: mu, denom, sum_inv_sigma_squared, sum_weighted_gain, g_tot, g_curr, sigma_curr, fknee, sigma_0, alpha
    real(dp), allocatable, dimension(:)     :: temp_gain, temp_invsigsquared
    real(dp), allocatable, dimension(:)     :: summed_invsigsquared, smoothed_gain
    real(dp), allocatable, dimension(:,:,:) :: g
    real(dp), allocatable, dimension(:,:)   :: pidrange_gainarr
    integer(i4b),   allocatable, dimension(:, :) :: window_sizes
    integer(i4b), save :: count = 0
    character(len=128)  :: kernel_type

    ndet       = tod%ndet
    nscan_tot  = tod%nscan_tot

    ! Collect all gain estimates on the root processor
    allocate(g(nscan_tot,ndet,2))
    allocate(temp_gain(nscan_tot))
    g = 0.d0
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          if (.not. tod%scans(i)%d(j)%accept) cycle
          g(k,j,1) = tod%scans(i)%d(j)%dgain
          g(k,j,2) = tod%scans(i)%d(j)%gain_invsigma
       end do
    end do
!    if (tod%myid == 0) then
       call mpi_allreduce(mpi_in_place, g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
            & tod%comm, ierr)
!    else
!       call mpi_reduce(g,            g, size(g), MPI_DOUBLE_PRECISION, MPI_SUM, &
!            & 0, tod%comm, ierr)
!    end if

!    if (tod%myid == 0) then
!!$       open(58,file='tmp.unf', form='unformatted')
!!$       read(58) g
!!$       close(58)

 !       count = count+1
        !write(*, *) "FREQ IS ", trim(tod%freq), count
       !nbin = nscan_tot / binsize + 1

!!$!        open(58,file='gain_' // trim(tod%freq) // '.dat', recl=1024)
!!$       do j = 1, ndet
!!$          do k = 1, nscan_tot
!!$             !if (g(k,j,2) /= 0) then
!!$             !if (g(k,j,2) /= g(k,j,2)) write(*,*) j,k, real(g(k,j,1),sp), real(g(k,j,2),sp), real(g(k,j,1)/g(k,j,2),sp)
!!$             if (g(k,j,2) > 0) then
!!$                if (abs(g(k, j, 1)) > 1e10) then
!!$                   write(*, *) 'G1'
!!$                   write(*, *) g(k, j, 1)
!!$                end if
!!$                if (abs(g(k, j, 2)) > 1e10) then
!!$                   write(*, *) 'G2'
!!$                   write(*, *) g(k, j, 2)
!!$                end if
!!$                !if (abs(dipole_mods(k, j) > 1e10)) then
!!$                !   write(*, *) 'DIPOLE_MODS'
!!$                !   write(*, *) dipole_mods(k, j)
!!$                !else
!!$                   write(58,*) j, k, real(g(k,j,1)/g(k,j,2),sp), real(g(k,j,1),sp), real(g(k,j,2),sp), real(dipole_mods(k, j), sp)
!!$                !end if
!!$             else
!!$                write(58,*) j, k, 0., 0.0, 0., 0.
!!$             end if
!!$          end do
!!$          write(58,*)
!!$       end do
!!$       close(58)

       do j = 1+tod%myid, ndet, tod%numprocs
         if (all(g(:, j, 1) == 0)) continue
          fknee = 0.002d0 / (60.d0 * 60.d0) ! In seconds
          alpha = -1.d0
          temp_gain = 0.d0
          ! This is not completely correct - should probably truncate, or
          ! something.
          do k = 1, nscan_tot
            if (g(k, j, 2) > 0.d0) then
               temp_gain(k) = g(k, j, 1) / g(k, j, 2)
            end if
          end do
          sigma_0 = calc_sigma_0(temp_gain)
!          sigma_0 = 0.002d0
          call wiener_filtered_gain(g(:, j, 1), g(:, j, 2), sigma_0, alpha, &
             & fknee, trim(tod%operation)=='sample', handle)
       end do
!    end if

    ! Distribute and update results
       do j = 1, ndet
          call mpi_bcast(g(:,j,:), size(g(:,j,:)),  MPI_DOUBLE_PRECISION, mod(j-1,tod%numprocs), tod%comm, ierr)    
       end do
    do j = 1, ndet
       do i = 1, tod%nscan
          k        = tod%scanid(i)
          !if (g(k, j, 2) <= 0.d0) cycle
          tod%scans(i)%d(j)%dgain = g(k,j,1)
          tod%scans(i)%d(j)%gain  = tod%gain0(0) + tod%gain0(j) + g(k,j,1)
          !write(*,*) j, k,  tod%gain0(0), tod%gain0(j), g(k,j,1), tod%scans(i)%d(j)%gain 
       end do
    end do

!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(g)

  end subroutine sample_smooth_gain


  ! Haavard: Remove current monopole fit, and replace inv_sigmasq with a proper invN(alpha,fknee) multiplication
!  subroutine accumulate_abscal(self, scan, det, mask, s_sub, &
!       & s_orb, A_abs, b_abs)
   ! This is implementing equation 16, adding up all the terms over all the sums
   ! the sum i is over the detector.
   subroutine accumulate_abscal(tod, scan, mask, s_sub, s_ref, s_invN, A_abs, b_abs, handle, out, s_highres, mask_lowres, tod_arr)
    implicit none
    class(comm_tod),                   intent(in)     :: tod
    integer(i4b),                      intent(in)     :: scan
    real(sp),          dimension(:,:), intent(in)     :: mask, s_sub, s_ref
    real(sp),          dimension(:,:), intent(in)     :: s_invN
    real(dp),          dimension(:),   intent(inout)  :: A_abs, b_abs
    type(planck_rng),                  intent(inout)  :: handle
    logical(lgt), intent(in) :: out
    real(sp),          dimension(:,:), intent(in), optional :: s_highres
    real(sp),          dimension(:,:), intent(in), optional :: mask_lowres
    real(sp),          dimension(:,:), intent(in), optional :: tod_arr
 
    real(sp), allocatable, dimension(:,:)     :: residual
    real(sp), allocatable, dimension(:)       :: r_fill
    real(dp)     :: A, b, scale
    integer(i4b) :: i, j, ext(2), ndet, ntod
    character(len=5) :: itext

    ndet = tod%ndet
    ntod = size(s_sub,1)
!    p = 5000 ! Buffer width
!    q = size(mask)-p
    call tod%downsample_tod(s_sub(:,1), ext)    
    allocate(residual(ext(1):ext(2),ndet), r_fill(size(s_sub,1)))
    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) then
          residual(:,j) = 0.
          cycle
       end if
       if (present(tod_arr)) then
         r_fill = tod_arr(:,j)-s_sub(:,j) - tod%scans(scan)%d(j)%baseline
         !if (tod%scanid(scan) == 30 .and. out) write(*,*) tod%scanid(scan), sum(abs(tod_arr(:,j))), sum(abs(s_sub(:,j))), tod%scans(scan)%d(j)%baseline
       else
         r_fill = tod%scans(scan)%d(j)%tod - s_sub(:,j) - tod%scans(scan)%d(j)%baseline
         !if (tod%scanid(scan) == 30 .and. out) write(*,*) tod%scanid(scan), sum(abs(tod%scans(scan)%d(j)%tod)), sum(abs(s_sub(:,j))), tod%scans(scan)%d(j)%baseline
       end if
       call fill_all_masked(r_fill, mask(:,j), ntod, trim(tod%operation) == 'sample', abs(real(tod%scans(scan)%d(j)%N_psd%sigma0, sp)), handle, tod%scans(scan)%chunk_num)
       call tod%downsample_tod(r_fill, ext, residual(:,j))
    end do

    call multiply_inv_N(tod, scan, residual, sampfreq=tod%samprate_lowres, pow=0.5d0)

    do j = 1, ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle
       if (present(mask_lowres)) then
          A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j)    * mask_lowres(:,j))
          b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j) * mask_lowres(:,j))
       else
          A_abs(j) = A_abs(j) + sum(s_invN(:,j) * s_ref(:,j))
          b_abs(j) = b_abs(j) + sum(s_invN(:,j) * residual(:,j))
       end if
       if (tod%scanid(scan) == 30 .and. out) then
         !write(*,*) tod%scanid(scan), real(sum(s_invN(:,j) * residual(:,j))/sum(s_invN(:,j) * s_ref(:,j)),sp), real(1/sqrt(sum(s_invN(:,j) * s_ref(:,j))),sp), '  # absK', j
         !write(*,*) tod%scanid(scan), sum(abs(s_invN(:,j))), sum(abs(residual(:,j))), sum(abs(s_ref(:,j))), '  # absK', j
       end if
    end do

    if (.false. .and. out) then
    !if (sum(s_invN(:,4) * residual(:,4) * mask_lowres(:,4))/sum(s_invN(:,4) * s_ref(:,4) * mask_lowres(:,4)) < -0.5d0) then
       call int2string(tod%scanid(scan), itext)
       !write(*,*) 'gain'//itext//'   = ', tod%gain0(0) + tod%gain0(1), tod%gain0(0), tod%gain0(1)
       open(58,file='gainfit3_'//itext//'.dat')
       do i = ext(1), ext(2)
          write(58,*) i-ext(1)+1, residual(i,1)
       end do
       write(58,*)
       do i = 1, size(s_ref,1)
          write(58,*) i, s_ref(i,1)
       end do
       write(58,*)
       scale = sum(s_invN(:,1) * residual(:,1)) / sum(s_invN(:,1) * s_ref(:,1))
       residual(:,1) = residual(:,1) - scale * s_ref(:,1)
       do i = 1, size(s_ref,1)
          write(58,*) i+ext(1)-1, residual(i+ext(1)-1,1) 
       end do
       if (present(mask_lowres)) then
          write(58,*)
          do i = 1, size(s_ref,1)
             write(58,*) i+ext(1)-1, mask_lowres(i,1) 
          end do
       end if
       close(58)

       open(58,file='gainfit4_'//itext//'.dat')       
       do i = 1, size(s_sub,1)
          if (present(tod_arr)) then
            write(58,*) i, tod_arr(i, 4) - tod%scans(scan)%d(4)%baseline
          else
            write(58,*) i, tod%scans(scan)%d(4)%tod(i) - tod%scans(scan)%d(4)%baseline
          end if
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, r_fill(i)
       end do
       write(58,*)
       if (present(s_highres)) then
          do i = 1, size(s_sub,1)
             write(58,*) i, s_highres(i,4)
          end do
          write(58,*)
       end if
       do i = 1, size(s_sub,1)
          write(58,*) i, s_sub(i,4)
       end do
       write(58,*)
       do i = 1, size(s_sub,1)
          write(58,*) i, mask(i,4)*10
       end do
       close(58)
    end if



    !if (det == 1) write(*,*) tod%scanid(scan), b/A, '  # absgain', real(A,sp), real(b,sp)

    deallocate(residual, r_fill)

  end subroutine accumulate_abscal

  ! Sample absolute gain from orbital dipole alone 
  subroutine sample_abscal_from_orbital(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (tod%myid == 0) then
       tod%gain0(0) = sum(b)/sum(A)
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%gain0(0) = tod%gain0(0) + 1.d0/sqrt(sum(A)) * rand_gauss(handle)
       end if
       if (tod%verbosity > 1) then
         write(*,*) 'abscal = ', tod%gain0(0)
         !write(*,*) 'sum(b), sum(A) = ', sum(b), sum(A)
       end if
    end if
    call mpi_bcast(tod%gain0(0), 1,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(A, b)

  end subroutine sample_abscal_from_orbital

  ! Sample absolute gain from orbital dipole alone 
  ! Eirik: Change this routine to sample the constant offset per radiometer; put the result in tod%gain0(i)
  subroutine sample_relcal(tod, handle, A_abs, b_abs)
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, k, ierr
    integer(i4b), allocatable, dimension(:) :: ind
    real(dp),     allocatable, dimension(:) :: A, b, rhs, x, tmp
    real(dp),     allocatable, dimension(:, :) :: coeff_matrix

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet), rhs(tod%ndet+1), x(tod%ndet+1), ind(tod%ndet+1), tmp(tod%ndet+1))
    allocate(coeff_matrix(tod%ndet+1, tod%ndet+1))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    coeff_matrix = 0.d0
    rhs          = 0.d0
    x            = 0.d0
    ind          = 0     ! List of indices for active (non-zero) elements in linear system
    k            = 0
    if (tod%myid == 0) then
       do j = 1, tod%ndet
         coeff_matrix(j, j) = A(j)
         rhs(j) = b(j) 
         if (trim(tod%operation) == 'sample') rhs(j) = rhs(j) + sqrt(A(j)) * rand_gauss(handle)
         coeff_matrix(j, tod%ndet+1) = 0.5d0
         coeff_matrix(tod%ndet+1, j) = 1
         if (coeff_matrix(j,j) > 0.d0) then
            k      = k+1
            ind(k) = j
         end if
       end do
       k = k+1
       ind(k) = tod%ndet+1
       coeff_matrix(tod%ndet+1, tod%ndet+1) = 0.d0
       rhs(tod%ndet+1) = 0.d0
       call solve_system_real(coeff_matrix(ind(1:k),ind(1:k)), tmp(1:k), rhs(ind(1:k)))
       x(ind(1:k)) = tmp(1:k)
       if (tod%verbosity > 1) then
         !write(*,*) 'A =', A
         !write(*,*) 'b =', b
         write(*,*) 'relcal = ', real(x,sp)
       end if
    end if
    call mpi_bcast(x, tod%ndet+1, MPI_DOUBLE_PRECISION, 0, &
       & tod%info%comm, ierr)

    tod%gain0(1:tod%ndet) = x(1:tod%ndet)
    deallocate(coeff_matrix, rhs, A, b, x, ind, tmp)

    do j = 1, tod%nscan
       do i = 1, tod%ndet
          tod%scans(j)%d(i)%gain = tod%gain0(0) + tod%gain0(i) + tod%scans(j)%d(i)%dgain 
       end do
    end do

  end subroutine sample_relcal

  subroutine sample_imbal_cal(tod, handle, A_abs, b_abs)
    !  Subroutine to sample the transmission imbalance parameters, defined in
    !  the WMAP data model as the terms x_im; given the definition
    !  d_{A/B} = T_{A/B} \pm Q_{A/B} cos(2 gamma_{A/B}) \pm U_{A/B} sin(2 gamma_{A/B})
    !  we have
    !  d = g[(1+x_im)*d_A - (1-x_im)*d_B]
    !  Returns x_{im,1} for detectors 13/14, and x_{im,2} for detectors 23/24.
    !
    !
    !  Arguments (fixed):
    !  ------------------
    !  A_abs: real(dp)
    !     Accumulated A_abs = s_ref^T N^-1 s_ref for all scans
    !  b_abs: real(dp)
    !     Accumulated b_abs = s_ref^T N^-1 s_sub for all scans
    !
    !  
    !  Arguments (modified):
    !  ---------------------
    !  tod: comm_WMAP_tod
    !     The entire tod object. tod%x_im estimated and optionally sampled
    !  handle: planck_rng derived type 
    !     Healpix definition for random number generation
    !     so that the same sequence can be resumed later on from that same point
    !
    implicit none
    class(comm_tod),                   intent(inout)  :: tod
    type(planck_rng),                  intent(inout)  :: handle
    real(dp),            dimension(:), intent(in)     :: A_abs, b_abs

    integer(i4b) :: i, j, ierr
    real(dp), allocatable, dimension(:) :: A, b

    ! Collect contributions from all cores
    allocate(A(tod%ndet), b(tod%ndet))
    call mpi_reduce(A_abs, A, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)
    call mpi_reduce(b_abs, b, tod%ndet, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & tod%info%comm, ierr)

    ! Compute gain update and distribute to all cores
    if (tod%myid == 0) then
       tod%x_im(1) = sum(b(1:2))/sum(A(1:2))
       tod%x_im(3) = sum(b(3:4))/sum(A(3:4))
       if (trim(tod%operation) == 'sample') then
          ! Add fluctuation term if requested
          tod%x_im(1) = tod%x_im(1) + 1.d0/sqrt(sum(A(1:2))) * rand_gauss(handle)
          tod%x_im(3) = tod%x_im(3) + 1.d0/sqrt(sum(A(3:4))) * rand_gauss(handle)
       end if
       tod%x_im(2) = tod%x_im(1)
       tod%x_im(4) = tod%x_im(3)
       if (tod%verbosity > 1) then
         write(*,*) 'b', sum(b(1:2)), sum(b(3:4))
         write(*,*) 'A', sum(A(1:2)), sum(A(3:4))
         write(*,*) 'imbal =', tod%x_im(1), tod%x_im(3)
       end if
    end if
    call mpi_bcast(tod%x_im, 4,  MPI_DOUBLE_PRECISION, 0, &
         & tod%info%comm, ierr)


    deallocate(A, b)

  end subroutine sample_imbal_cal

  subroutine get_smoothing_windows(tod, windows, dipole_mods)
     implicit none

     integer(i4b), dimension(:, :), intent(out)         :: windows
     class(comm_tod),               intent(in)          :: tod
     real(dp), dimension(:, :),     intent(in)          :: dipole_mods

     real(sp), allocatable, dimension(:)                :: nonzero_mask
     integer(i4b)  :: window_size_dipole_minimum, window_size_dipole_maximum
     real(dp)   :: low_dipole_thresh, curr_dipole
     real(dp)   :: high_dipole_thresh

     real(dp)   :: reduce_fac 

     integer(i4b)   :: i, j, window_size


     high_dipole_thresh = 5d-6
     low_dipole_thresh = 1d-6
     reduce_fac = 1.d0

     select case (trim(tod%freq))
         case ('030')
            window_size_dipole_minimum = 1200
            window_size_dipole_maximum = 400
         case ('044')
            window_size_dipole_minimum = 1500
            window_size_dipole_maximum = 300
         case ('070')
            window_size_dipole_minimum = 1880
            window_size_dipole_maximum = 400
         case ('023-WMAP_K')
             window_size_dipole_minimum = 1 !120
             window_size_dipole_maximum = 1 ! 40
         case default
            window_size_dipole_minimum = 1 !1800
            window_size_dipole_maximum = 1 !400
      end select

      do i = 1, tod%nscan_tot
         do j = 1, tod%ndet
            curr_dipole = dipole_mods(i, j)
            if (curr_dipole == 0) then
               windows(j, i) = 0
            else if (curr_dipole < low_dipole_thresh) then
               windows(j, i) = window_size_dipole_minimum
            else if (curr_dipole > high_dipole_thresh) then
               windows(j, i) = window_size_dipole_maximum
            else
               ! Interpolate
               windows(j, i) = window_size_dipole_minimum + &
                  & int(real(window_size_dipole_maximum - & 
                  & window_size_dipole_minimum, dp) / (high_dipole_thresh - &
                  & low_dipole_thresh) * (curr_dipole - low_dipole_thresh))
            end if
            windows(j, i) = windows(j, i) * reduce_fac
         end do
      end do

  end subroutine get_smoothing_windows


  subroutine wiener_filtered_gain(b, inv_N_wn, sigma_0, alpha, fknee, sample, &
     & handle)
     implicit none

     real(dp), dimension(:), intent(inout)   :: b
     real(dp), dimension(:), intent(in)      :: inv_N_wn
     real(dp), intent(in)                    :: sigma_0, alpha, fknee
     logical(lgt), intent(in)                :: sample
     type(planck_rng)                        :: handle

     real(dp), allocatable, dimension(:)     :: freqs, dt, inv_N_corr
     complex(dpc), allocatable, dimension(:) :: dv
     real(dp), allocatable, dimension(:)     :: fluctuations, temp
     real(dp), allocatable, dimension(:)     :: precond
     complex(dpc), allocatable, dimension(:) :: fourier_fluctuations
     integer*8          :: plan_fwd, plan_back
     integer(i4b)       :: nscan, nfft, n, nomp, err
     real(dp)           :: samprate

     integer(i4b)   :: i


     nscan = size(b)
!     nfft = 2 * nscan
     nfft = nscan
     n = nfft / 2 + 1
     samprate = 1.d0 / (60.d0 * 60.d0) ! Just assuming a pid per hour for now
     allocate(freqs(n-1))
     do i = 1, n - 1
        freqs(i) = i * (samprate * 0.5) / (n - 1)
     end do

     nomp = 1
     call dfftw_init_threads(err)
     call dfftw_plan_with_nthreads(nomp)
     allocate(dt(nfft), dv(0:n-1))
     call dfftw_plan_dft_r2c_1d(plan_fwd, nfft, dt, dv, fftw_estimate + fftw_unaligned)
     call dfftw_plan_dft_c2r_1d(plan_back, nfft, dv, dt, fftw_estimate + fftw_unaligned)
     ! We use temporary arrays inside the transformation function, so no need
     ! for these anymore
     deallocate(dt, dv)

     allocate(inv_N_corr(n))
     allocate(fluctuations(nscan))
!     allocate(temp(nscan))
     allocate(precond(nscan))
     allocate(fourier_fluctuations(n))

     write(*, *) 'Sigma_0: ', sigma_0
     write(*, *) 'alpha: ', alpha
     write(*, *) 'fknee: ', fknee

     inv_N_corr = calculate_invcov(sigma_0, alpha, fknee, freqs)
     if (sample) then
        do i = 1, nscan
            fluctuations(i) = rand_gauss(handle)
        end do
        call timev_to_fourier(fluctuations, fourier_fluctuations, plan_fwd)
        do i = 1, n
            fourier_fluctuations(i) = fourier_fluctuations(i) * sqrt(inv_N_corr(i))
        end do
        call fourierv_to_time(fourier_fluctuations, fluctuations, plan_back)

        do i = 1, nscan
            fluctuations(i) = fluctuations(i) + sqrt(inv_N_wn(i)) * rand_gauss(handle)
            b(i) = b(i) + fluctuations(i)
            precond(i) = inv_N_wn(i) + 1 / sigma_0 ** 2
         end do
      end if

!      temp = solve_cg_gain(inv_N_wn, inv_N_corr, b, precond, plan_fwd, plan_back, .true.)
!      precond = 1.d0
!      b = solve_cg_gain(inv_N_wn, inv_N_corr, b, precond, plan_fwd, plan_back, .false.)
      b = solve_cg_gain(inv_N_wn, inv_N_corr, b, precond, plan_fwd, plan_back)
      deallocate(inv_N_corr, freqs, fluctuations, fourier_fluctuations, precond)
          
      call dfftw_destroy_plan(plan_fwd)                                           
      call dfftw_destroy_plan(plan_back) 

  end subroutine wiener_filtered_gain

  subroutine timev_to_fourier(vector, fourier_vector, plan_fwd)
     implicit none
     real(dp), dimension(:), intent(in)         :: vector
     complex(dpc), dimension(:), intent(out)    :: fourier_vector
     integer*8             , intent(in)     :: plan_fwd

     real(dp), allocatable, dimension(:)                 :: dt
     complex(dpc), allocatable, dimension(:)             :: dv
     integer(i4b)   :: nfft, n

!     nfft = size(vector) * 2
     nfft = size(vector)
     n = nfft / 2 + 1

     allocate(dt(nfft), dv(0:n-1))
     dt(1:nfft) = vector
!     dt(nfft:nfft/2+1:-1) = dt(1:nfft/2)
     call dfftw_execute_dft_r2c(plan_fwd, dt, dv)
     !Normalization
     dv = dv / sqrt(real(nfft, dp))
     fourier_vector(1:n) = dv(0:n-1)
     deallocate(dt, dv)

  end subroutine timev_to_fourier

  subroutine fourierv_to_time(fourier_vector, vector, plan_back)
     implicit none
     complex(dpc), dimension(:), intent(in)     :: fourier_vector
     real(dp), dimension(:), intent(out)    :: vector
     integer*8             , intent(in)     :: plan_back

     real(dp), allocatable, dimension(:)                 :: dt
     complex(dpc), allocatable, dimension(:)             :: dv
     integer(i4b)   :: nfft, n

     nfft = size(vector)
     n = nfft / 2 + 1

     allocate(dt(nfft), dv(0:n-1))
     dv(0:n-1) = fourier_vector(1:n)
     call dfftw_execute_dft_c2r(plan_back, dv, dt)
     ! Normalization
     dt = dt / sqrt(real(nfft, dp))
     vector = dt(1:nfft)
     deallocate(dt, dv)

  end subroutine fourierv_to_time

  function solve_cg_gain(inv_N_wn, inv_N_corr, b, precond, plan_fwd, plan_back)!, &
!     & with_precond)
     implicit none
     real(dp), dimension(:), intent(in) :: inv_N_wn, inv_N_corr, b, precond
     integer*8             , intent(in) :: plan_fwd, plan_back
!     logical(lgt)                       :: with_precond

     real(dp), dimension(size(b))       :: solve_cg_gain
     real(dp), allocatable, dimension(:)    :: initial_guess, prop_sol, residual
     real(dp), allocatable, dimension(:)    :: p, Abyp, new_residual, z, new_z
     real(dp)       :: alpha, beta, orig_residual
     integer(i4b)   :: iterations, nscan, i
     logical(lgt)   :: converged
     character(len=8)   :: itext

!      open(58, file='gain_cg_invcov.dat')
!      do i = 1, size(inv_N_corr)
!         write(58, *) inv_N_corr(i)
!      end do
!      close(58)
!      open(58, file='gain_cg_invnwn.dat')
!      do i = 1, size(inv_N_wn)
!         write(58, *) inv_N_wn(i)
!      end do
!      close(58)


     nscan = size(b)
     allocate(initial_guess(nscan), prop_sol(nscan), residual(nscan), p(nscan))
     allocate(Abyp(nscan), new_residual(nscan), z(nscan), new_z(nscan))

     initial_guess = 0.d0
     prop_sol = initial_guess
     residual = b - tot_mat_mul_by_vector(inv_N_wn, inv_N_corr, prop_sol, &
        & plan_fwd, plan_back)
     orig_residual = sum(abs(residual))
     z = residual / precond
!      open(58, file='gain_cg_residual.dat')
!      do i = 1, nscan
!         write(58, *) residual(i)
!      end do
!      close(58)

     p = z
     iterations = 0
     converged = .false.
     do while ((.not. converged) .and. (iterations < 100000))
         Abyp = tot_mat_mul_by_vector(inv_N_wn, inv_N_corr, p, plan_fwd, &
            & plan_back, filewrite=(iterations == 0))
!         if (iterations == 0) then
!            open(58, file='gain_cg_abyp.dat')
!            do i = 1, nscan
!               write(58, *) Abyp(i)
!            end do
!            close(58)
!         end if

         alpha = sum(residual * z) / sum(p * Abyp) 
         prop_sol = prop_sol + alpha * p
!         if (iterations == 0) then
!            write(*,*) 'Alpha:', alpha
!            open(58, file='gain_cg_prop_sol.dat')
!            do i = 1, nscan
!               write(58, *) prop_sol(i)
!            end do
!            close(58)
!         end if

         new_residual = residual - alpha * Abyp
!         if (iterations == 0) then
!            open(58, file='gain_cg_new_residual.dat')
!            do i = 1, nscan
!               write(58, *) new_residual(i)
!            end do
!            close(58)
!         end if

         if (sum(abs(new_residual))/orig_residual < 1e-12) then 
            converged = .true.
            exit
         end if
         new_z = new_residual / precond
         beta = sum(new_residual * new_z) / sum(residual * z)
         p = new_z + beta * p
!         if (iterations == 0) then
!            write(*,*) 'Beta:', beta
!            open(58, file='gain_cg_p.dat')
!            do i = 1, nscan
!               write(58, *) p(i)
!            end do
!            close(58)
!         end if

         residual = new_residual
         z = new_z
         iterations = iterations + 1
         if (mod(iterations, 100) == 0) then
            write(*, *) "Gain CG search res: ", sum(abs(new_residual))/orig_residual, sum(abs(new_residual))
!            call int2string(iterations, itext)
!            open(58, file='gain_cg_' // itext // '.dat')
!            do i = 1, nscan
!               write(58, *) prop_sol(i)
!            end do
!            close(58)
         end if
!         if (iterations == 1) then
!            open(58, file='gain_cg_' // itext // '.dat')
!            call int2string(iterations, itext)
!            do i = 1, nscan
!               write(58, *) prop_sol(i)
!            end do
!            close(58)
!         end if
     end do
     if (.not. converged) then
        write(*, *) "Gain CG search did not converge."
     end if
!     if (with_precond) then
!        write(*, *) "With preconditioner"
!     else
!        write(*, *) "Without preconditioner"
!     end if
     write(*, *) "Gain CG iterations: ", iterations
     solve_cg_gain = prop_sol

     deallocate(initial_guess, prop_sol, residual, p, Abyp, new_residual, z, new_z)

  end function solve_cg_gain

  function tot_mat_mul_by_vector(time_mat, fourier_mat, vector, plan_fwd, &
     & plan_back, filewrite)
     implicit none

     real(dp), dimension(:)                            :: vector
     real(dp), dimension(size(vector))                 :: tot_mat_mul_by_vector
     real(dp), dimension(size(vector)), intent(in)     :: time_mat
     real(dp), dimension(:) , intent(in)     :: fourier_mat
     integer*8              , intent(in)     :: plan_fwd, plan_back
     logical(lgt)           , optional       :: filewrite

     logical(lgt)       :: write_file
     integer(i4b)       :: i
     real(dp), dimension(size(vector))       :: temp_vector
     complex(dpc), dimension(size(fourier_mat))        :: fourier_vector

     write_file = .false.
     if (present(filewrite)) then
         write_file = filewrite
      end if

     if (all(vector .eq. 0)) then
        tot_mat_mul_by_vector = 0.d0
        return
     end if
!     if (write_file) then
!         open(58, file='gain_cg_vector_in.dat')
!         do i = 1, size(vector)
!            write(58, *) vector(i)
!         end do
!         close(58)
!      end if

     call timev_to_fourier(vector, fourier_vector, plan_fwd)
!     if (write_file) then
!         open(58, file='gain_cg_fourier_vector.dat')
!         do i = 1, size(fourier_vector)
!            write(58, *) abs(fourier_vector(i))
!         end do
!         close(58)
!      end if

     fourier_vector = fourier_vector * fourier_mat
!     if (write_file) then
!         open(58, file='gain_cg_matmul_fourier_vector.dat')
!         do i = 1, size(fourier_vector)
!            write(58, *) abs(fourier_vector(i))
!         end do
!         close(58)
!      end if


     call fourierv_to_time(fourier_vector, temp_vector, plan_back)
!     if (write_file) then
!         open(58, file='gain_cg_vector_after_fourier.dat')
!         do i = 1, size(temp_vector)
!            write(58, *) temp_vector(i)
!         end do
!         close(58)
!      end if

!     if (write_file) then
!         open(58, file='gain_cg_vector_time_matmul.dat')
!         do i = 1, size(vector)
!            write(58, *) vector(i) * time_mat(i)
!         end do
!         close(58)
!      end if

     tot_mat_mul_by_vector = vector * time_mat + temp_vector

  end function tot_mat_mul_by_vector

  function calculate_invcov(sigma_0, alpha, fknee, freqs)
     implicit none
     real(dp), dimension(:), intent(in)     :: freqs
     real(dp), intent(in)                   :: sigma_0, fknee, alpha
     real(dp), dimension(size(freqs)+1)     :: calculate_invcov

     calculate_invcov(1) = 0.d0
     calculate_invcov(2:size(freqs)+1) = 1.d0 / (sigma_0 ** 2 * (freqs/fknee) ** alpha)

  end function calculate_invcov

  function calc_sigma_0(gain)
     implicit none
     real(dp)       :: calc_sigma_0
     real(dp), dimension(:) :: gain

     real(dp), dimension(size(gain))    :: res
     real(dp)   :: std, mean
     integer(i4b)   :: i

     mean = gain(1) 
     res = 0.d0
     do i = 2, size(gain)
         res(i) = gain(i) - gain(i-1)
         mean = mean + res(i)
      end do
      mean = mean / size(gain)
      std = 0.d0
      do i = 2, size(gain)
         std = std + (res(i) - mean) ** 2
      end do
      if (std == 0) then
         open(58, file='gain_cg_std_gain.dat')
         do i = 1, size(gain)
            write(58, *) gain(i)
         end do
         close(58)
      end if
         
      calc_sigma_0 = sqrt(std / (2 * (size(gain) - 1)))
  end function calc_sigma_0

end module comm_tod_gain_mod
