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
module comm_huffman_mod
  use comm_utils
  implicit none

  private
  public huffcode, huffman_decode, huffman_decode2, hufmak, get_bitstring, hufmak_precomp, huffman_decode3


  !TODO: this needs to be generalized somehow to support a sp symbols case 
  type huffcode
     integer(i4b) :: nch, nodemax
     integer(i4b), allocatable, dimension(:) :: icode,left,iright,ncode,symbols,nfreq
  end type huffcode

  interface huffman_decode2
     procedure huffman_decode2_int, huffman_decode2_sp
  end interface huffman_decode2

contains

  ! Public routines
  subroutine huffman_decode2_int(hcode, x_in, x_out, imod)
    implicit none
    type(huffcode),               intent(in)  :: hcode
    byte,           dimension(:), intent(in)  :: x_in
    integer(i4b),   dimension(:), intent(out) :: x_out
    integer(i4b),                 intent(in), optional :: imod

    integer(i4b) :: i, j, k, n, nb, ich, l,nc,node
    integer(i4b), allocatable, dimension(:) :: buf

    n  = size(x_out)
    
    i = 2 ! Byte counter
    j = 7 ! Bit counter
    do k = 1, n
       node=hcode%nodemax
       do while (node > hcode%nch)
          if (btest(x_in(i),j)) then 
             node=hcode%iright(node) 
          else
             node=hcode%left(node)
          end if
          j = j-1
          if (j == -1) then
             i = i+1
             j = 7
          end if
       end do
       x_out(k) = hcode%symbols(node) 
       if (k > 1)         x_out(k) = x_out(k-1) + x_out(k)
       if (present(imod)) x_out(k) = iand(x_out(k),imod)
    end do

  end subroutine huffman_decode2_int

  !subroutine huffman_decode_vlen(hcode, x_in, x_out)
  !  implicit none
  !  type(huffcode),               intent(in)  :: hcode
  !  
  !  real(sp),       dimension(:), intent(out) :: x_out


  !end subroutine huffman_decode_vlen

  subroutine huffman_decode2_sp(hcode, x_in, x_out)
    implicit none
    type(huffcode),               intent(in)  :: hcode
    byte,           dimension(:), intent(in)  :: x_in
    real(sp),       dimension(:), intent(out) :: x_out

    integer(i4b) :: i, j, k, n, nb, ich, l,nc,node, curr, prev
    integer(i4b), allocatable, dimension(:) :: buf

    n  = size(x_out)
   
    write(*,*) x_in(1:100), hcode%nodemax, hcode%symbols(1:100)
 
    i = 2 ! Byte counter
    j = 7 ! Bit counter
    curr = 0
    do k = 1, n
       node=hcode%nodemax
       do while (node > hcode%nch)
          if (btest(x_in(i),j)) then 
             node=hcode%iright(node) 
          else
             node=hcode%left(node)
          end if
          j = j-1
          if (j == -1) then
             i = i+1
             j = 7
          end if
       end do
       x_out(k) = hcode%symbols(node)
       if (k > 1) x_out(k) = x_out(k-1) + x_out(k)
    end do

  end subroutine huffman_decode2_sp


  ! Public routines
  subroutine huffman_decode4(hcode, x_in, x_out, imod)
    implicit none
    type(huffcode),               intent(in)  :: hcode
    byte,           dimension(:), intent(in)  :: x_in
    integer(i4b),   dimension(:), intent(out) :: x_out
    integer(i4b),                 intent(in), optional :: imod

    integer(i4b) :: i, j, k, n, nb, ich, l,nc,node
    integer(i4b), allocatable, dimension(:) :: buf

!!$    if (.not. present(offset)) then
!!$      offset_ = 0
!!$    else
!!$      offset_ = offset
!!$    end if

    n  = size(x_out)
    k = 1
    node=hcode%nodemax
    do i = 2, size(x_in)  ! First byte does not contain real data
       do j = 7, 0, -1
          if (btest(x_in(i),j)) then 
             node=hcode%iright(node) 
          else
             node=hcode%left(node)
          end if
          if (node <= hcode%nch) then
             x_out(k) = hcode%symbols(node)
             if (k > 1)         x_out(k) = x_out(k-1) + x_out(k)
             if (present(imod)) x_out(k) = iand(x_out(k),imod)
             k    = k + 1
             if (k > n) return
             node = hcode%nodemax
          end if
       end do
    end do

  end subroutine huffman_decode4

  ! Public routines
  subroutine huffman_decode3(hcode, x_in, x_out, imod, offset)
    implicit none
    type(huffcode),               intent(in)  :: hcode
    byte,           dimension(:), intent(in)  :: x_in
    integer(i4b),   dimension(:), intent(out) :: x_out
    integer(i4b),                 intent(in), optional :: imod
    integer(i4b),                 intent(in), optional :: offset 
    !an offset from the start of the chunk, data from offset:offset+len is
    !returned

    integer(i4b) :: i, j, k, n, nb, ich, l,nc,node, offset_
    integer(i4b), allocatable, dimension(:) :: buf

    if (.not. present(offset)) then
      offset_ = 0
    else
      offset_ = offset
    end if

    n  = size(x_out)
    allocate(buf(offset_ + n))

!!$    nb = 8       ! First byte does not contain real data
!!$    do i = 1, n
!!$       node=hcode%nodemax
!!$       do 
!!$          nc=shifta(nb,3)+1    !nc=nb/8+1
!!$          l=iand(nb,7)         !l=mod(nb,8) 
!!$          nb=nb+1
!!$          if (btest(x_in(nc),7-l)) then 
!!$             node=hcode%iright(node) 
!!$          else
!!$             node=hcode%left(node)
!!$          end if
!!$          if (node <= hcode%nch) then
!!$             x_out(i) = hcode%symbols(node)
!!$             if (i > 1) x_out(i) = x_out(i-1) + x_out(i)
!!$             if (present(imod)) x_out(i) = iand(x_out(i),imod)
!!$             exit
!!$          end if
!!$       end do
!!$    end do

    k = 1
    node=hcode%nodemax
    do i = 2, size(x_in)  ! First byte does not contain real data
       do j = 7, 0, -1
          if (btest(x_in(i),j)) then 
             node=hcode%iright(node) 
          else
             node=hcode%left(node)
          end if
          if (node == 0) then ! case where there was a single entry ie. all 0
             buf(k) = 0
             k = k+1
             if (k > n+offset_) then
               x_out(1:n) = buf(1+offset_:n+offset_)
               return
             end if
             node = hcode%nodemax
          else if (node <= hcode%nch) then
             buf(k) = hcode%symbols(node)
             if (k > 1)         buf(k) = buf(k-1) + buf(k)
             if (present(imod)) buf(k) = iand(buf(k),imod)
             k    = k + 1
             if (k > n+offset_) then
               x_out(1:n) = buf(1+offset_:n+offset_)
               return
             end if
             node = hcode%nodemax
          end if
       end do
    end do

  end subroutine huffman_decode3


  ! Public routines
  subroutine huffman_decode(hcode, x_in, x_out)
    implicit none
    type(huffcode),               intent(in)  :: hcode
    byte,           dimension(:), intent(in)  :: x_in
    integer(i4b),   dimension(:), intent(out) :: x_out

    integer(i4b) :: i, n, nb, ich

    n  = size(x_out)
    nb = 8       ! First byte does not contain real data
    do i = 1, n
       call hufdec(ich, x_in, nb, hcode)
       x_out(i) = hcode%symbols(ich)
    end do

  end subroutine huffman_decode

  subroutine hufmak_precomp(symbols,tree,hcode)
    implicit none
    integer(i4b), dimension(:), intent(in) :: symbols, tree
    type(huffcode) :: hcode

    integer(i4b) :: n
    
    call huff_deallocate(hcode)
    n             = size(symbols)
    hcode%nch     = n
    hcode%nodemax = tree(1)
    allocate(hcode%left(2*n-1), hcode%iright(2*n-1), hcode%symbols(n))
    hcode%left    = 0
    hcode%iright  = 0
    hcode%symbols = symbols
    hcode%left(n+1:2*n-1)   = tree(2:n)
    hcode%iright(n+1:2*n-1) = tree(n+1:2*n-1)

  end subroutine hufmak_precomp

  subroutine hufmak(symbols,nfreq,hcode)
    implicit none
    integer(i4b), dimension(:), intent(in)  :: symbols,nfreq
    type(huffcode) :: hcode
    integer(i4b) :: ibit,j,k,n,node,nused,nerr,itmp(1), ilong, nlong
    integer(i4b), dimension(2*size(nfreq)-1) :: indx,iup,nprob
    iup=0; indx=0; nprob=0
    hcode%nch=size(nfreq)
    call huff_allocate(hcode,size(nfreq))
    hcode%symbols = symbols
    hcode%nfreq = nfreq
    nused=0
    nprob(1:hcode%nch)=nfreq(1:hcode%nch)
    call array_copy(pack(arth(1,1,hcode%nch), nfreq(1:hcode%nch) /= 0 ),&
         indx,nused,nerr)
    do j=nused,1,-1 
       call hufapp(j)
    end do
    k=hcode%nch
    do 
       if (nused <= 1) exit
       node=indx(1)
       indx(1)=indx(nused)
       nused=nused-1
       call hufapp(1)
       k=k+1
       nprob(k)=nprob(indx(1))+nprob(node)
       hcode%left(k)=node        
       hcode%iright(k)=indx(1)
       iup(indx(1))=-k 
       iup(node)=k 
       indx(1)=k
       call hufapp(1)
    end do
    hcode%nodemax=k
    iup(hcode%nodemax)=0
    do j=1,hcode%nch 
       if (nprob(j) /= 0) then
          n=0
          ibit=0
          node=iup(j)
          do
             if (node == 0) exit
             if (node < 0) then
                n=ibset(n,ibit)
                node=-node
             end if
             node=iup(node)
             ibit=ibit+1
          end do
          hcode%icode(j)=n
          hcode%ncode(j)=ibit
       end if
    end do
    itmp=maxloc(hcode%ncode(1:hcode%nch))
    ilong=itmp(1)
    nlong=hcode%ncode(ilong)
    if (nlong > bit_size(1_i4b)) then
       write(*,*) 'Huffman error: Number of possible bits for code exceeded'
       stop
    end if
  contains
    subroutine hufapp(l)
      implicit none
      integer(i4b), intent(in) :: l
      integer(i4b) :: i,j,k,n
      n=nused
      i=l
      k=indx(i)
      do
         if (i > n/2) exit
         j=i+i
         if (j < n) then
            if (nprob(indx(j)) > nprob(indx(j+1))) j=j+1
         end if
         if (nprob(k) <= nprob(indx(j))) exit
         indx(i)=indx(j)
         i=j
      end do
      indx(i)=k
    end subroutine hufapp
  end subroutine hufmak
  
  function get_bitstring(hcode, i)
    implicit none
    type(huffcode) :: hcode
    integer(i4b)   :: i
    character(len=hcode%ncode(i)) :: get_bitstring

    integer(i4b) :: j

    do j = 1, hcode%ncode(i)
       if (btest(hcode%icode(i),j-1)) then
          get_bitstring(j:j) = '1'
       else
          get_bitstring(j:j) = '0'
       end if
    end do

  end function get_bitstring

  ! routines
  subroutine huff_allocate(hcode,mc)
    implicit none
    type(huffcode) :: hcode
    integer(i4b) :: mc
    integer(i4b) :: mq
    mq=2*mc-1
    allocate(hcode%icode(mq),hcode%ncode(mq),hcode%left(mq),hcode%iright(mq),hcode%symbols(mc),hcode%nfreq(mc))
    hcode%icode(:)=0
    hcode%ncode(:)=0
  end subroutine huff_allocate

  subroutine huff_deallocate(hcode)
    implicit none
    type(huffcode) :: hcode
    if (allocated(hcode%iright))  deallocate(hcode%iright)
    if (allocated(hcode%left))    deallocate(hcode%left)
    if (allocated(hcode%symbols)) deallocate(hcode%symbols)
    if (allocated(hcode%icode))   deallocate(hcode%icode)
    if (allocated(hcode%ncode))   deallocate(hcode%ncode)
    if (allocated(hcode%nfreq))   deallocate(hcode%nfreq)
  end subroutine huff_deallocate


!!$  subroutine hufenc(ich,codep,nb,hcode)
!!$    implicit none
!!$    integer(i4b), intent(in) :: ich
!!$    integer(i4b), intent(inout) :: nb
!!$    character(1), allocatable, dimension(:) :: codep
!!$    type(huffcode) :: hcode
!!$    integer(i4b) :: k,l,n,nc,ntmp
!!$    do n=hcode%ncode(k),1,-1 
!!$       nc=nb/8+1 
!!$       if (nc > size(codep)) codep=>reallocate(codep,2*size(codep))
!!$       l=mod(nb,8)
!!$       if (l == 0) codep(nc)=char(0)
!!$       if (btest(hcode%icode(k),n-1)) then 
!!$          ntmp=ibset(ichar(codep(nc)),l)
!!$          codep(nc)=char(ntmp)
!!$       end if
!!$       nb=nb+1
!!$    end do
!!$  end subroutine hufenc
  
  subroutine hufdec(ich,code,nb,hcode)
    implicit none
    integer(i4b), intent(out) :: ich
    integer(i4b), intent(inout) :: nb
    byte, dimension(:), intent(in) :: code
    type(huffcode) :: hcode
    integer(i4b) :: l,nc,node
    node=hcode%nodemax
    do 
       nc=nb/8+1
       if (nc > size(code)) then 
          ich=hcode%nch 
          return
       end if
       l=mod(nb,8) 
       nb=nb+1
       if (btest(code(nc),7-l)) then 
          node=hcode%iright(node) 
       else
          node=hcode%left(node)
       end if
       if (node <= hcode%nch) then 
          ich=node
          return
       end if
    end do
  end subroutine hufdec

  subroutine array_copy(src,dest,n_copied,n_not_copied)
    implicit none
    integer(i4b), dimension(:), intent(in)  :: src
    integer(i4b), dimension(:), intent(out) :: dest
    integer(i4b),               intent(out) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  end subroutine array_copy

  function arth(first,increment,n)
    integer(i4b), intent(in) :: first,increment
    integer(i4b), intent(in) :: n
    integer(i4b), dimension(n) :: arth
    integer(i4b) :: k
    if (n > 0) arth(1)=first
    do k=2,n
       arth(k)=arth(k-1)+increment
    end do
  end function arth

end module comm_huffman_mod
