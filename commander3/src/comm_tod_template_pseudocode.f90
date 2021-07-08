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
module comm_tod_TEMPLATE_mod
   subroutine process_TEMPLATE_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)

      main_it: do main_iter = 1, n_main_iter
         if (self%myid == 0) write(*,*) '  Performing main iteration = ', main_iter
         ! Select operations for current iteration
         do_oper(samp_acal)    = (main_iter == n_main_iter-3) 
         do_oper(samp_rcal)    = (main_iter == n_main_iter-2) 
         do_oper(samp_G)       = (main_iter == n_main_iter-1) 
         do_oper(samp_N)       = (main_iter >= n_main_iter-0) 
         do_oper(bin_map)      = (main_iter == n_main_iter  )

         dipole_mod = 0

         if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
            A_abscal = 0.d0; b_abscal = 0.d0
         end if

         do i = 1, self%nscan

            ! Decompress pointing, psi and flags for current scan
            do j = 1, ndet
               call self%decompress_pointing_and_flags(i, j, pix(:, j, :), &
                    & psi(:, j, :), flag(:, j))
            end do

            ! Construct sky signal template
            call project_sky(self, map_sky(:, :, :, 1), pix, psi, flag, &
                   & i, s_sky, mask)

            ! Construct orbital dipole template
            call self%orb_dp%p%compute_orbital_dipole_4pi(i, pix(:,:,1), psi(:,:,1), s_orb)

            ! Add orbital dipole to total signal
            s_buf = 0.d0
            do j = 1, ndet
               s_tot(:, j) = s_sky(:, j) + s_orb(:,j)
               s_buf(:, j) = s_tot(:, j)
            end do

            ! Determine which data chunks are "bad"
            if (main_iter == 1 .and. self%first_call) then
               do j = 1, ndet
                  self%scans(i)%d(j)%accept = .true.
                  if (all(mask(:,j) == 0)) self%scans(i)%d(j)%accept = .false.
                  if (self%scans(i)%d(j)%sigma0 <= 0.d0) self%scans(i)%d(j)%accept = .false.
               end do
            end if

            ! Precompute filtered signal for calibration
            if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. do_oper(samp_acal)) then
               call self%downsample_tod(s_orb(:,1), ext)
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  if (do_oper(samp_G) .or. do_oper(samp_rcal) .or. .not. self%orb_abscal) then
                     s_buf(:,j) = s_tot(:,j)
                     call fill_all_masked(s_buf(:,j), mask(:,j), ntod, trim(self%operation)=='sample', real(self%scans(i)%d(j)%sigma0, sp), handle, self%scans(i)%chunk_num)
                     call self%downsample_tod(s_buf(:,j), ext, &
                          & s_lowres(:,j))
                  else
                     call self%downsample_tod(s_orb(:,j), ext, &
                          & s_lowres(:,j))
                  end if
               end do
               s_invN = s_lowres
               call multiply_inv_N(self, i, s_invN,   sampfreq=self%samprate_lowres, pow=0.5d0)
               call multiply_inv_N(self, i, s_lowres, sampfreq=self%samprate_lowres, pow=0.5d0)
            end if

            ! Prepare for absolute calibration
            if (do_oper(samp_acal) .or. do_oper(samp_rcal)) then
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
                     s_buf(:,j) = real(self%gain0(0) + self%scans(i)%d(j)%dgain,sp) * s_tot(:, j)
                  end if
               end do
               call accumulate_abscal(self, i, mask, s_buf, s_invN, A_abscal, b_abscal, handle)

            end if

            ! Fit gain
            if (do_oper(samp_G)) then
               call calculate_gain_mean_std_per_scan(self, i, s_invN, mask, s_tot, handle)
            end if

            ! Fit correlated noise
            if (do_oper(samp_N)) then
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  if (do_oper(samp_mono)) then
                     s_buf(:,j) = s_tot(:,j)-s_mono(:,j)
                  else
                     s_buf(:,j) = s_tot(:,j)
                  end if
               end do
               call sample_n_corr(self, handle, i, mask, s_buf, n_corr, pix(:,:,1), .false.)
            else
               n_corr = 0.
            end if

            ! Compute noise spectrum
            if (do_oper(samp_N_par)) then
               call sample_noise_psd(self, handle, i, mask, s_tot, n_corr)
            end if

            ! Get calibrated map
            if (do_oper(bin_map)) then
               d_calib = 0
               do j = 1, ndet
                  if (.not. self%scans(i)%d(j)%accept) cycle
                  inv_gain = 1.0/real(self%scans(i)%d(j)%gain, sp)
                  d_calib(1, :, j) = (self%scans(i)%d(j)%tod - n_corr(:, j))* &
                     & inv_gain - s_tot(:, j) + s_sky(:, j)
               end do

               ! Bin the calibrated map
               call bin_TOD(self, d_calib, pix,  &
                      & psi, flag, self%x_im, sprocmask%a, b_map, A_map, i)
            end if

            do j = 1, ndet
               if (.not. self%scans(i)%d(j)%accept) cycle
               dipole_mod(self%scanid(i), j) = masked_variance(s_sky(:, j), mask(:, j))
            end do

         end do

         if (do_oper(samp_acal)) then
            call sample_abscal_from_orbital(self, handle, A_abscal, b_abscal)
         end if

         if (do_oper(samp_rcal)) then
            call sample_relcal(self, handle, A_abscal, b_abscal)
         end if

         if (do_oper(samp_G)) then
            call sample_smooth_gain(self, handle, dipole_mod)
         end if

      end do main_it


      call finalize_binned_map(self, handle, A_map, b_map, rms_out, outmaps=outmaps)

   end subroutine process_TEMPLATE_tod

end module comm_tod_TEMPLATE_mod
