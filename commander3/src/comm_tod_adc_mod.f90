!===============================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University
! of Oslo.
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
!===============================================================================

! This module handles corrections to time ordered data due to adc issues

module comm_tod_adc_mod

  use spline_1D_mod


  implicit none

  private
  public comm_adc, adc_pointer

  type :: comm_adc
    real(dp), dimension(:), allocatable :: adc_in, adc_out
  contains
    procedure :: adc_correct
  end type comm_adc

  interface comm_adc
    procedure constructor
  end interface comm_adc

  type adc_pointer
    class(comm_adc), pointer :: p => null()
  end type adc_pointer

contains

   function constructor(adc_i, adc_o)
     ! ====================================================================
     ! Sets up an adc correction object that maps input and output voltages
     !
     ! Inputs:
     ! 
     ! adc_i : float array
     !   The array of input voltages
     !
     ! adc_o : float array
     !   The array of output voltages

     implicit none
     real(dp), dimension(:), intent(in)          :: adc_i, adc_o
     class(comm_adc), pointer                    :: constructor    

     allocate(constructor)

     allocate(constructor%adc_in(size(adc_i)), constructor%adc_out(size(adc_o)))

     constructor%adc_in  = adc_i
     constructor%adc_out = adc_o
 
   end function constructor

   subroutine adc_correct(self, tod_in, correct_tod)
     !=========================================================================
     ! Adc corrects a timestream 
     ! 
     ! Inputs:
     !
     ! self : comm_adc object
     !    Defines the adc correction that should be applied
     ! tod_in : float array
     !    The tod that is to be corrected
     ! 
     ! Outputs : 
     !
     ! correct_tod : float array
     !    The adc corrected version of tod_in

     implicit none
     class(comm_adc),                 intent(inout) :: self
     real(sp), dimension(:),          intent(in)    :: tod_in
     real(sp), dimension(:),          intent(out)   :: correct_tod
     type(spline_type)                              :: sadc
     integer(i4b)                                   :: i, len

     len = size(self%adc_in)

     call spline(sadc, self%adc_in, self%adc_out)

     write(*,*) tod_in

     stop

     !TODO: figure out the correct algorithm and implement it
     !-------------------------------------------------------

     ! To start we'll the just spline the DPC adc correction table
     ! Must check units!!

     ! correct_tod = splint(sadc,tod_in)
     correct_tod = tod_in

   end subroutine adc_correct

end module comm_tod_adc_mod
