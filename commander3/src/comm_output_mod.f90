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
module comm_output_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_chisq_mod
  use comm_hdf_mod
  implicit none

contains

  subroutine init_chain_file(cpar, iter)
    implicit none
    
    type(comm_params), intent(in)  :: cpar
    integer(i4b),      intent(out) :: iter

    integer(i4b)                 :: i, j, hdferr, ierr
    logical(lgt)                 :: exist
    character(len=4)             :: ctext
    character(len=6)             :: itext
    character(len=512)           :: postfix, chainfile, hdfpath, fg_file
    type(hdf_file)   :: file
    TYPE(h5o_info_t) :: object_info

    call int2string(cpar%mychain, ctext)
    chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
         & '_c' // trim(adjustl(ctext)) // '.h5'

    ! Delete existing chain file if necessary; create new file if necessary; open file
    iter = 1
    if (cpar%myid_chain == 0) then
       inquire(file=trim(chainfile), exist=exist)
       if (trim(cpar%chain_status) == 'new' .or. .not. exist) then
          if (exist) call rm(trim(chainfile))
          call open_hdf_file(chainfile, file, 'w')

          ! testing hdf parameter output
          call write_params_to_hdf(cpar, file)

          call close_hdf_file(file)
          iter = -1
       else if (trim(cpar%chain_status) == 'append') then
          call open_hdf_file(chainfile, file, 'r')
          exist = .true.
          do while (exist)
             call int2string(iter, itext)
             call h5eset_auto_f(0, hdferr)
             call h5oget_info_by_name_f(file%filehandle, itext, object_info, hdferr)
             exist = (hdferr == 0)
             if (exist) iter = iter+1
          end do
          iter = max(1,iter-1)
          write(*,*) '  Continuing chain '//ctext// ' on iteration ', iter
          call close_hdf_file(file)          
       else
          write(*,*) 'Unsupported chain mode =', trim(cpar%chain_status)
          call mpi_finalize(ierr)
          stop
       end if

       !delete fg_ind_mean_cXXXX.dat if it exists
       fg_file=trim(cpar%outdir)//'/fg_ind_mean_c' // trim(adjustl(ctext))//'.dat'
       inquire(file=fg_file, exist=exist)
       if (exist) call rm(trim(fg_file))


    end if
    call mpi_bcast(iter, 1, MPI_INTEGER, 0, cpar%comm_chain, ierr)

  end subroutine init_chain_file

  subroutine output_FITS_sample(cpar, iter, output_hdf)
    implicit none
    
    type(comm_params), intent(in) :: cpar
    integer(i4b),      intent(in) :: iter
    logical(lgt),      intent(in) :: output_hdf

    integer(i4b)                 :: i, j, hdferr, ierr, unit, p_min, p_max
    real(dp)                     :: chisq, t1, t2, t3, t4, theta_sum
    logical(lgt)                 :: exist, init, new_header
    character(len=4)             :: ctext
    character(len=6)             :: itext
    character(len=512)           :: postfix, chainfile, hdfpath, fg_file, temptxt, fg_txt
    character(len=2048)          :: outline, fg_header
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map),     pointer :: map => null(), chisq_map => null(), chisq_sub => null()
    class(comm_comp),    pointer :: c => null()
    class(comm_N),      pointer :: N => null()
    type(hdf_file) :: file
    TYPE(h5o_info_t) :: object_info

    call update_status(status, "output_start")
    
    call int2string(cpar%mychain, ctext)
    call int2string(iter,         itext)
    postfix = 'c'//ctext//'_k'//itext
    chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
         & '_c' // trim(adjustl(ctext)) // '.h5'

!!$    ! Delete existing chain file if necessary; create new file if necessary; open file
!!$    if (first_call .and. cpar%myid_chain == 0 .and. output_hdf) then
!!$       if (trim(cpar%chain_status) == 'new') then
!!$          inquire(file=trim(chainfile), exist=exist)
!!$          if (exist) call rm(trim(chainfile))
!!$       else if (trim(cpar%chain_status) == 'append') then
!!$          
!!$       else
!!$          write(*,*) 'Unsupported chain mode =', trim(cpar%chain_status)
!!$          call mpi_finalize(ierr)
!!$          stop
!!$       end if
!!$       inquire(file=trim(chainfile), exist=exist)
!!$       if (.not. exist) then
!!$          call open_hdf_file(chainfile, file, 'w')
!!$          call close_hdf_file(file)
!!$       end if
!!$       first_call = .false.
!!$    end if
    if (cpar%myid_chain == 0 .and. output_hdf) then
       inquire(file=trim(chainfile), exist=exist)
       call open_hdf_file(chainfile, file, 'b')
       ! Delete group if it already exists
       call h5eset_auto_f(0, hdferr)
       call h5oget_info_by_name_f(file%filehandle, trim(adjustl(itext)), object_info, hdferr)
       if (hdferr == 0) call h5gunlink_f(file%filehandle, trim(adjustl(itext)), hdferr)
       !write(*,*) 'group ', trim(adjustl(itext))
       call create_hdf_group(file, trim(adjustl(itext)))
       !if (.not. cpar%resamp_CMB) call create_hdf_group(file, trim(adjustl(itext))//'/md')
       call create_hdf_group(file, trim(adjustl(itext))//'/md')
    end if
    call update_status(status, "output_chain")

    !Prepare mean foregrounds values print to file
    if (cpar%myid_chain == 0) then
       fg_file=trim(cpar%outdir)//'/fg_ind_mean_c' // trim(adjustl(ctext))//'.dat'
       inquire(file=fg_file, exist=exist)
       unit = getlun()       
       if (.not. exist) then
          new_header=.true.
          open(unit, file=trim(fg_file), status='new',action='write', recl=1024)
       else
          new_header=.false.
          open(unit, file=trim(fg_file), status='old', position='append',action='write', recl=1024)
       end if

       if (new_header) write(fg_header,fmt='(a10)') '# Sample'
       write(outline,fmt='(a10)') itext
    end if


    ! Output component results
    c => compList
    call wall_time(t1)
    do while (associated(c))
!       if (.not. c%output .or. (cpar%resamp_CMB .and. trim(c%type) /= 'cmb')) then
       if (.not. c%output) then
          c => c%next()
          cycle
       end if
       call c%dumpFITS(iter, file, output_hdf, postfix, cpar%outdir)
       select type (c)
       class is (comm_diffuse_comp)
          if (output_hdf) then
             hdfpath = trim(adjustl(itext))//'/'//trim(adjustl(c%label))
             call c%Cl%writeFITS(cpar%mychain, iter, hdffile=file, hdfpath=hdfpath)
          else
             call c%Cl%writeFITS(cpar%mychain, iter)
          end if

          !get mean values (and component labels) for fg mean print
          do i = 1,c%npar
             if (cpar%myid_chain == 0) then 
                if (new_header) then
                   fg_txt=''
                   if (c%poltype(i) == 1) then
                      if (cpar%only_pol) then
                         if (c%nmaps > 1) write(fg_txt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(QU)'
                      else
                         write(fg_txt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(IQU)'
                      end if
                   else if (c%poltype(i) == 2) then
                      if (cpar%only_pol) then
                         if (c%nmaps > 1) then
                            write(fg_txt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(QU)'
                         end if
                      else
                         if (c%nmaps > 1) then
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(I)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(QU)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                         else
                            write(fg_txt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(I)'
                         end if
                      end if
                   else if (c%poltype(i) == 3) then
                      if (cpar%only_pol) then
                         if (c%nmaps > 2) then
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(Q)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(U)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                         else if (c%nmaps > 1) then
                            write(fg_txt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(Q)'
                         end if
                      else
                         if (c%nmaps > 2) then
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(I)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(Q)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(U)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                         else if (c%nmaps > 1) then
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(I)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                            write(temptxt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(Q)'
                            fg_txt=trim(fg_txt)//trim(temptxt)
                         else
                            write(fg_txt,fmt='(a20)') trim(c%label)//'_'//trim(c%indlabel(i))//'(I)'
                         end if
                      end if
                   end if
                end if
                if (.not. trim(fg_txt)=='') fg_header = trim(fg_header)//trim(fg_txt)
             end if

             do j = 1,c%poltype(i)
                if (c%poltype(i) == 1) then
                   p_min = 1; p_max = c%nmaps
                   if (cpar%only_pol) p_min = 2
                else if (c%poltype(i) == 2) then
                   if (j == 1) then
                      if (cpar%only_pol) cycle
                      p_min = 1; p_max = 1
                   else
                      p_min = 2; p_max = c%nmaps
                   end if
                else if (c%poltype(i) == 3) then
                   if (j == 1 .and. cpar%only_pol) cycle
                   p_min = j
                   p_max = j
                end if

                !all maps in the same poltype has the same theta, only evaluate p_min
                if (p_min > c%nmaps) cycle
                call mpi_reduce(sum(c%theta(i)%p%map(:,p_min)), theta_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, cpar%comm_chain, ierr)
                if (cpar%myid_chain == 0) then 
                   write(temptxt,fmt='(e20.7)') theta_sum/c%theta(i)%p%info%npix
                   outline = trim(outline)//trim(temptxt)
                end if
             end do
          end do

       end select
       call update_status(status, "output_"//trim(c%label))
       c => c%next()
    end do

    

    if (cpar%resamp_CMB) then
       if (cpar%myid_chain == 0 .and. output_hdf) call close_hdf_file(file)    
       return
    end if

    ! Output instrumental parameters
    call output_inst_params(cpar, file, itext, &
         & trim(cpar%outdir)//'/instpar_'//trim(postfix)//'.dat', output_hdf)
    call update_status(status, "output_inst")

    ! Output channel-specific residual maps
    if (cpar%output_residuals .or. cpar%output_chisq) then
       if (cpar%output_chisq) then
          info      => comm_mapinfo(cpar%comm_chain, cpar%nside_chisq, 0, cpar%nmaps_chisq, cpar%pol_chisq)
          chisq_map => comm_map(info)
       end if
       do i = 1, numband
          call wall_time(t3)
          map => compute_residual(i)
          call update_status(status, "output_res1_"//trim(data(i)%label))
          !call data(i)%apply_proc_mask(map)
          !map%map = map%map * data(i)%mask%map ! Apply frequency mask to current residual
          if (cpar%output_residuals) then
             call map%writeFITS(trim(cpar%outdir)//'/res_'//trim(data(i)%label)//'_'// &
                  & trim(postfix)//'.fits')
             call wall_time(t4)
          end if
          call update_status(status, "output_res2_"//trim(data(i)%label))
          if (cpar%output_chisq) then
             call data(i)%N%sqrtInvN(map)
             map%map = map%map**2
             
             info  => comm_mapinfo(data(i)%info%comm, chisq_map%info%nside, 0, data(i)%info%nmaps, data(i)%info%nmaps==3)
             chisq_sub => comm_map(info)
             call map%udgrade(chisq_sub)
             do j = 1, data(i)%info%nmaps
                chisq_map%map(:,j) = chisq_map%map(:,j) + chisq_sub%map(:,j) * (map%info%npix/chisq_sub%info%npix)
             end do
             call chisq_sub%dealloc(); deallocate(chisq_sub)
          end if
          call map%dealloc(); deallocate(map)
          call update_status(status, "output_res3_"//trim(data(i)%label))
       end do
       
       if (cpar%output_chisq) then
          call mpi_reduce(sum(chisq_map%map), chisq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, cpar%comm_chain, ierr)
          call chisq_map%writeFITS(trim(cpar%outdir)//'/chisq_'// trim(postfix) //'.fits')
          if (cpar%myid_chain == 0) write(*,fmt='(a,i4,a,e16.8)') &
               & '    Chain = ', cpar%mychain, ' -- chisq = ', chisq
          call chisq_map%dealloc(); deallocate(chisq_map)
       end if
       call update_status(status, "output_chisq")
    end if

    ! get chisq for fg_mean file 
    if (cpar%myid_chain == 0) then
       if (new_header) fg_header=trim(fg_header)//'          full_chisq           avg_chisq       chisq_highlat      avg_reduced_chisq'
       write(temptxt,fmt='(e20.8,e20.8,a25,a25)') chisq, chisq/(12*cpar%nside_chisq**2), '(too be implemented)', '(too be implemented)'
       outline = trim(outline)//trim(temptxt)
       !need to find a nice way of only gathering high latitude chisq

       !write fg_mean info to file and close file
       if (new_header) write(unit,fmt='(a)') trim(fg_header)
       write(unit,fmt='(a)') trim(outline)
       close(unit)
    end if

    ! Output signal components per band
    if (cpar%output_sig_per_band) call output_signals_per_band(cpar%outdir, postfix)
    call update_status(status, "output_sig")

    ! Output TOD parameters
    if (output_hdf .and. cpar%enable_TOD_analysis) then
       if (cpar%myid_chain == 0) then
          call create_hdf_group(file, trim(adjustl(itext))//'/tod')
       end if
       do i = 1, numband  
          if (trim(data(i)%tod_type) == 'none') cycle
          !write(*,*) 'associated', i, associated(data(i)%tod)
          if (associated(data(i)%tod)) then
             N => data(i)%N
             select type (N)
             class is (comm_N_rms)
                call data(i)%tod%dumpToHDF(file, iter, &
                     & data(i)%map0, N%rms0)
             end select
          end if
       end do
    end if

    if (cpar%myid_chain == 0 .and. output_hdf) call close_hdf_file(file)    
  end subroutine output_FITS_sample

  subroutine output_inst_params(cpar, chainfile, iter, filename, output_hdf)
    implicit none
    type(comm_params), intent(in) :: cpar
    type(hdf_file),    intent(in) :: chainfile
    character(len=*),  intent(in) :: iter
    character(len=*),  intent(in) :: filename
    logical(lgt),      intent(in) :: output_hdf
    
    integer(i4b) :: unit, i, j
    real(dp), allocatable, dimension(:,:) :: bp_delta

    if (cpar%myid_chain /= 0) return

    if (output_hdf) then
       call create_hdf_group(chainfile, trim(adjustl(iter))//'/gain')
       call create_hdf_group(chainfile, trim(adjustl(iter))//'/bandpass')
    end if
    
    unit = getlun()
    open(unit, file=trim(filename), recl=1024)
    write(unit,*) '#       Band          Gain         delta_bp(0)%p'
    do i = 1, numband
       write(unit,'(a15,f12.5,f10.2)') adjustl(trim(data(i)%label)), data(i)%gain, data(i)%bp(0)%p%delta
       if (output_hdf) then
          call write_hdf(chainfile, trim(adjustl(iter))//'/gain/'//trim(adjustl(data(i)%label)), &
               & data(i)%gain)

          allocate(bp_delta(0:data(i)%ndet,data(i)%bp(0)%p%npar))
          do j = 0, data(i)%ndet
             bp_delta(j,:) = data(i)%bp(j)%p%delta
          end do
          call write_hdf(chainfile, trim(adjustl(iter))//'/bandpass/'//&
               & trim(adjustl(data(i)%label)//'_det'), bp_delta)
          deallocate(bp_delta)
       end if
    end do
    close(unit)

  end subroutine output_inst_params

  ! =============== testing hdf parameter output ===============
  subroutine write_params_to_hdf(cpar, chainfile)
      implicit none 
      type(comm_params), intent(in) :: cpar
      type(hdf_file),    intent(in) :: chainfile
      integer(i4b) :: n, i
      character(len=512) :: hdf_path, polarization

      call create_hdf_group(chainfile, 'parameters')
      n = size(cpar%cs_label)

      do i = 1, n
         hdf_path = 'parameters/'//trim(adjustl(cpar%cs_label(i)))
         call create_hdf_group(chainfile, trim(hdf_path))

         if (cpar%cs_polarization(i)) then
            polarization = "True"
         else
            polarization = "False"
         endif

         ! Model parameters required by cosmoglobe not included in the gibbs iter part tof h5 file
         call write_hdf_0d_char(chainfile, trim(adjustl(hdf_path))//'/type', trim(cpar%cs_type(i)))
         call write_hdf(chainfile, trim(adjustl(hdf_path))//'/class', trim(cpar%cs_class(i)))
         call write_hdf(chainfile, trim(adjustl(hdf_path))//'/nside', cpar%cs_nside(i))
         call write_hdf_0d_char(chainfile, trim(adjustl(hdf_path))//'/polarization', trim(polarization))
         call write_hdf(chainfile, trim(adjustl(hdf_path))//'/unit', trim(cpar%cs_unit(i)))
         call write_hdf(chainfile, trim(adjustl(hdf_path))//'/nu_ref', cpar%cs_nu_ref(i,:))
         call write_hdf(chainfile, trim(adjustl(hdf_path))//'/fwhm', cpar%cs_fwhm(i))
      end do 

  end subroutine write_params_to_hdf
  ! ============================================================
end module comm_output_mod
