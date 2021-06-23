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
module InvSamp_mod
  use healpix_types
  use rngmod
  use spline_1D_mod
  implicit none

  integer(i4b), parameter          :: INVSAMP_MAX_NUM_EVALS = 1000
  integer(i4b), parameter, private :: N_SPLINE              = 10000
  real(dp),     parameter, private :: DELTA_LNL             = 12.5d0 ! Five sigma
  integer(i4b), parameter, private :: MIN_NUM_ACTIVE_POINT  = 50
  real(dp),     parameter, private :: TOLERANCE             = 1d-2

contains

  function sample_InvSamp(handle, x_in, lnL, prior, status, n_eval, lnL_in, optimize, use_precomputed_grid, &
       & tolerance_)
    implicit none

    type(planck_rng)                               :: handle
    real(dp), dimension(1:), intent(in)            :: x_in
    real(dp)                                       :: sample_InvSamp
    real(dp), dimension(2),               optional :: prior
    integer(i4b),            intent(out), optional :: status, n_eval
    logical(lgt),            intent(in),  optional :: optimize, use_precomputed_grid
    real(dp), dimension(1:), intent(in),  optional :: lnL_in
    real(dp),                intent(in),  optional :: tolerance_
    interface
       function lnL(x)
         use healpix_types
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: lnL
       end function lnL
    end interface

    integer(i4b) :: i, j, n, m, iter, stat, x_peak(1), a, b
    logical(lgt) :: optimize_, use_precomputed_grid_
    real(dp)     :: prior_(2), x_new, y_new, y_new_spline, lnL_peak, epsilon
    real(dp)     :: dx, x_min, x_max, eta, tol, lnL0, dpr
    real(dp), dimension(INVSAMP_MAX_NUM_EVALS) :: x_n, x_spline
    real(dp), dimension(INVSAMP_MAX_NUM_EVALS) :: S_n, S_n2, S_spline
    real(dp), dimension(N_SPLINE)      :: x, P, F

    sample_InvSamp = 1.d30

    stat = 0
    if (present(prior)) then
       if (prior(2) < prior(1)) then
          write(*,*) 'InvSamp_mod -- error: Upper prior is lower than lower prior', prior
          stop
       end if
       prior_ = prior
    else
       prior_ = [-1d100, 1d100]
    end if
    dpr = prior_(2)-prior_(1)
    use_precomputed_grid_ = .false.
    if (present(use_precomputed_grid)) use_precomputed_grid_ = use_precomputed_grid
    tol = TOLERANCE; if (present(tolerance_)) tol = tolerance_
    if (use_precomputed_grid_) then
       n    = size(x_in)
       !write(*,*) "precompute: n=",n
       x_n(1:n) = x_in
       S_n(1:n) = lnL_in
    else
       n        = size(x_in)
       !write(*,*) "not precomputed: n=",n
       x_n(1:n) = x_in
       if (present(lnL_in)) then
          ! Use pre-computed values to start grid
          if (size(lnL_in) /= size(x_in)) then
             write(*,*) 'InvSamp_mod error: lnL_in and x_in have different lengths'
             stop
          end if
          S_n(1:n) = lnL_in
       else
          ! Evaluate initial points
          do i = 1, n
             S_n(i) = lnL(x_n(i))
          end do
       end if

       ! Check if priors are peak
       lnL0  = lnL(prior_(1))
       x_new = prior_(1) + 1d-6*dpr
       y_new = lnL(x_new)
       if (lnL0 > y_new) then
          call update_InvSamp_sample_set(prior_(1), lnL0, x_n, S_n, n, stat)
!          write(*,*) "n after update_InvSamp_sample_set prior 1: ", n
       end if

       lnL0  = lnL(prior_(2))
       x_new = prior_(1) - 1d-6*dpr
       y_new = lnL(x_new)
       if (lnL0 > y_new) then
          call update_InvSamp_sample_set(prior_(2), lnL0, x_n, S_n, n, stat)
!          write(*,*) "update_InvSamp_sample_set prior_(2)", stat
!          write(*,*) "n after update_InvSamp_sample_set: ", n
       end if


       ! Check that peak is bounded; if not do a golden ratio search
       if (x_n(1) /= prior_(1)) then
          do while (S_n(1) > S_n(2) .and. x_n(1) > prior_(1))
             !x_new = max(x_n(1) - 1.61803d0*(x_n(2)-x_n(1)), 0.5d0*(x_n(1)+prior_(1)))
             x_new = 0.5d0*(x_n(1)+prior_(1))
             y_new = lnL(x_new)
             call update_InvSamp_sample_set(x_new, y_new, x_n, S_n, n, stat)
!             write(*,*) "n after update_InvSamp_sample_set peak bound: ", n
          end do
       end if

       if (x_n(n) /= prior_(2)) then
          do while (S_n(n) > S_n(n-1) .and. x_n(n) < prior_(2))
             x_new = min(x_n(n) + 1.61803d0*(x_n(n)-x_n(n-1)), prior_(2))
             y_new = lnL(x_new)
             call update_InvSamp_sample_set(x_new, y_new, x_n, S_n, n, stat)
             !write(*,*) "n after update_InvSamp_sample_set peak bound 2: ", n
          end do
       end if
       if (stat /= 0) then
          sample_InvSamp = 1.d30
          return
       end if

       ! Check that we bound the 5-sigma range, ie., delta_chisq = 25 => delta lnL = 12.5
       lnL_peak = maxval(S_n(1:n))
       do while (lnL_peak-S_n(1) < DELTA_LNL .and. x_n(1) > prior_(1) + 1d-6*dpr)
          x_new = 0.5d0*(x_n(1)+prior_(1))
!          x_new = max(x_n(1) - 1.61803d0*(x_n(2)-x_n(1)), prior_(1))
          y_new = lnL(x_new)
!          if (maxval(S_n)-y_new < 100.d0) then
!             prior_(1) = x_new 
!          else
             call update_InvSamp_sample_set(x_new, y_new, x_n, S_n, n, stat)
!             write(*,*) "n after update_InvSamp_sample_set wrong1: ", n
!          end if
       end do

       do while (lnL_peak-S_n(n) < DELTA_LNL .and. x_n(n) < prior_(2) - 1.d-6*dpr)
          x_new = min(x_n(n) + 1.61803d0*(x_n(n)-x_n(n-1)), prior_(2))
          y_new = lnL(x_new)
!          if (maxval(S_n)-y_new < 100.d0) then
!             prior_(2) = x_new 
!          else
             call update_InvSamp_sample_set(x_new, y_new, x_n, S_n, n, stat)
!             write(*,*) "n after update_InvSamp_sample_set wrong 2: ", n
!          end if
       end do
       if (stat /= 0) then
          sample_InvSamp = 1.d30
          return
       end if

       !write(*,*) x_n
       !write(*,*) S_n

       ! Refine grid until we have a sufficiently accurate solution
       epsilon = 1.d30
       iter = 0
       do while (epsilon > tol)
          iter          = iter+1
          m             = n
          lnL_peak      = maxval(S_n(1:n))
          x_spline(1:n) = x_n(1:n)
          S_spline(1:n) = S_n(1:n)
          epsilon       = 0.d0
!!$          if (.true. .or. any(S_spline(1:n) /= S_spline(1:n))) then
!!$             open(78,file='S.dat')
!!$             do i = 1, n
!!$                write(78,*) x_spline(i), S_spline(i)
!!$             end do
!!$             close(78)
!!$          end if
          
          call spline(x_spline(1:n), S_spline(1:n), 1.d30, 1.d30, S_n2(1:n))
          do i = m, 2, -1
             if (lnL_peak-S_n(i-1) < DELTA_LNL .or. lnL_peak-S_n(i) < DELTA_LNL) then
                x_new        = 0.5d0 * (x_spline(i-1)+x_spline(i))
                y_new        = lnL(x_new)
                y_new_spline = splint(x_spline(1:m), S_spline(1:m), S_n2(1:m), x_new)
                epsilon      = max(abs(y_new-y_new_spline), epsilon)
                !write(*,*) "y_new-y_new_spline:", y_new-y_new_spline, tol
                if (abs(y_new-y_new_spline) > tol) then
                   call update_InvSamp_sample_set(x_new, y_new, x_n, S_n, n, stat)
                   !write(*,*) "x_new, y_new, x_n, S_n, n, stat: ", x_new, y_new, x_n, S_n, n, stat
                   !write(*,*) "n after update_InvSamp_sample_set DELTA_LNL: ", n, i, m
                end if
             end if
             if (stat /= 0) exit
          end do
!          write(*,*) iter, epsilon, n
          if (iter > 100) then
             stat = stat+1
             open(69,file='p_inv.dat')
             do i = 1, n
                write(69,*) x_n(i), S_n(i)
             end do
             close(69)
             stop
             exit
          end if
          if (stat /= 0) exit
       end do
    end if

!    write(*,*) x_n(1:n)
!    write(*,*) S_n(1:n)

    if (.false.) then
       write(*,*) 'n = ', n
       write(*,*) 'x', real(x_n(1:n),sp)
       write(*,*) 'S', real(S_n(1:n),sp)
       write(*,*)
       
       open(60,file='P_InvSamp.dat')
       do i = 1, n
          write(60,*) x_n(i), S_n(i)
       end do
       close(60)
       stop
    end if

    ! Spline probability function
    if (stat == 0) then
       call spline(x_n(1:n), S_n(1:n), 1.d30, 1.d30, S_n2(1:n))

       optimize_ = .false.; if (present(optimize)) optimize_ = optimize
       if (optimize_) then
          x_peak         = maxloc(S_n(1:n))
          sample_InvSamp = x_n(x_peak(1))
       else
          
          lnL_peak = maxval(S_n(1:n))
          a        = 1
          b        = n
          do while (lnL_peak-S_n(a+1) > DELTA_LNL .and. S_n(a+1) > S_n(a))
             a = a+1
          end do
          do while (lnL_peak-S_n(b-1) > DELTA_LNL .and. S_n(b-1) > S_n(b))
             b = b-1
          end do
          x_min = x_n(a)
          x_max = x_n(b)
          dx    = (x_max-x_min) / (N_SPLINE-1.d0)
          do i = 1, N_SPLINE
             x(i) = x_min + dx * (i-1.d0)
             P(i) = splint(x_n(1:n), S_n(1:n), S_n2(1:n), x(i))
          end do
          P = exp(P-maxval(P))
          
          ! Compute cumulative distribution
          F(1) = 0.d0
          F(2) = dx * 0.5d0*(P(1) + P(2))
          do j = 3, N_SPLINE
             !F(j) = F(j-2) + dx * (P(j-2)/3.d0 + 4.d0*P(j-1)/3.d0 + P(j)/3.d0)
             F(j) = F(j-1) + dx * 0.5d0*(P(j-1) + P(j))
          end do
          P = P / F(N_SPLINE)
          F = F / F(N_SPLINE)
          
          if (.false.) then
             open(61,file='P_InvSamp.dat')
             do j = 1, N_SPLINE
                write(61,*) x(j), P(j)
             end do
             close(61) 
             
             open(61,file='F_InvSamp.dat')
             do j = 1, N_SPLINE
                write(61,*) x(j), F(j)
             end do
             close(61) 
             !stop
          end if
          
          ! Draw a uniform variate between 0 and 1
          eta = rand_uni(handle)
          do while (eta < F(1) .or. eta > F(N_SPLINE))
             eta = rand_uni(handle)
          end do
          
          ! Solve F(x) ~ F(x_1) + (x-x_1) * (F(x_2)-F(x_1))/(x_2 - x_1) = eta
          i = 2
          do while (eta > F(i) .and. i < N_SPLINE)
             i = i+1
          end do
          if (i == N_SPLINE) then
             sample_InvSamp = x(N_SPLINE)
          else
             sample_InvSamp = x(i-1) + (eta-F(i-1)) * (x(i)-x(i-1))/(F(i)-F(i-1))
          end if
          
       end if

       if (sample_InvSamp /= sample_InvSamp) then
          stat = stat+1
          sample_InvSamp = 1.d30
       end if
    end if
    
    if (stat == 0) then
       sample_InvSamp = max(min(sample_InvSamp, prior_(2)), prior_(1))
    else
       sample_InvSamp = 1.d30
    end if

    if (present(n_eval)) n_eval = n
    if (present(status)) status = stat
    

  end function sample_InvSamp


  subroutine update_InvSamp_sample_set(x_new, y_new, x, y, n, stat)
    implicit none

    real(dp),                   intent(in)    :: x_new, y_new
    real(dp),     dimension(:), intent(inout) :: x, y
    integer(i4b),               intent(inout) :: n, stat

    integer(i4b) :: i

    if (n == INVSAMP_MAX_NUM_EVALS) then
       !write(*,*) 'InvSamp_mod -- sample set full. Increase INVSAMP_MAX_NUM_EVALS.'
       stat = stat+1
!       write(*,*) "INVSAMP_MAX_NUM_EVALS", INVSAMP_MAX_NUM_EVALS
!       write(*,*) "x_new, y_new, x, y, n, stat:", x_new, y_new, x, y, n, stat
       !stop
       return
    end if

    do i = 1, n
       if (x_new == x(i)) return
       if (x_new < x(i)) exit
    end do
    x(i+1:n+1) = x(i:n); x(i) = x_new
    y(i+1:n+1) = y(i:n); y(i) = y_new
    n      = n+1

  end subroutine update_InvSamp_sample_set

end module InvSamp_mod
