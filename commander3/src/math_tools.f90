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
module math_tools
  use healpix_types
  USE OMP_LIB


  interface invert_matrix
     module procedure invert_matrix_dpc, invert_matrix_dp, invert_matrix_sp
  end interface

  interface invert_matrix_with_mask
     module procedure invert_matrix_with_mask_dpc, invert_matrix_with_mask_dp
  end interface

  interface convert_fract2sigma
     module procedure convert_fract2sigma_sp, convert_fract2sigma_dp
  end interface

contains

  subroutine invert_matrix_dpc(matrix)
    implicit none

    complex(dpc), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    complex(dpc), allocatable, dimension(:)   :: work

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(n))

    call ZGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'ZGETRF: LU factorization failed. Info = ', info
       stop
    else

       call ZGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'ZGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_dpc


  subroutine invert_matrix_dp(matrix, cholesky, status, ln_det)
    implicit none

    real(dp), dimension(1:,1:), intent(inout)         :: matrix
    logical(lgt),               intent(in),  optional :: cholesky
    integer(i4b),               intent(out), optional :: status
    real(dp),                   intent(out), optional :: ln_det

    integer(i4b)     :: i, j, n, lda, info, lwork
    logical(lgt)     :: use_cholesky
    character(len=1) :: uplo
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(dp),     allocatable, dimension(:)   :: work

    if(present(status)) status = 0
    use_cholesky = .false.; if (present(cholesky)) use_cholesky = cholesky
    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    uplo  = 'l'
    allocate(ipiv(n))
    allocate(work(n))

    if (use_cholesky) then
       call DPOTRF(uplo, n, matrix, lda, info)
       if (present(ln_det)) then
          ln_det = 0.d0
          do i = 1, n
             if (matrix(i,i) > 0.d0) then
                ln_det = ln_det + 2.d0*log(matrix(i,i))
             else
                ln_det = -1.d30
                exit
             end if
          end do
       end if
    else
       call DGETRF(n, n, matrix, lda, ipiv, info)
    end if
    if (info /= 0) then
       if(present(status)) then
          status = info
          return
       end if
       write(*,*) 'DGETRF: Factorization failed. Info = ', info
       stop
    else

       if (use_cholesky) then
          call DPOTRI(uplo, n, matrix, lda, info)
       else
          call DGETRI(n, matrix, lda, ipiv, work, lwork, info)
       end if

       if (info /= 0) then
          if(present(status)) then
             status = info
             return
          end if
          write(*,*) 'DGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    if (use_cholesky) then
       do i = 1, n
          do j = i+1, n
             matrix(i,j) = matrix(j,i)
          end do
       end do
    end if

    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_dp

!!$  subroutine invert_matrix_dp(matrix)
!!$    implicit none
!!$
!!$    real(dp), dimension(1:,1:), intent(inout) :: matrix
!!$
!!$    integer(i4b)     :: i, j, n, lda, info, lwork
!!$    integer(i4b), allocatable, dimension(:)   :: ipiv
!!$    real(dp),     allocatable, dimension(:)   :: work
!!$
!!$    n     = size(matrix(1,:))
!!$    lda   = n
!!$!    lwork = n
!!$    lwork = n*64
!!$    allocate(ipiv(n))
!!$!    allocate(work(n))
!!$    allocate(work(lwork))
!!$
!!$    call DGETRF(n, n, matrix, lda, ipiv, info)
!!$    if (info /= 0) then
!!$       write(*,*) 'DGETRF: LU factorization failed. Info = ', info
!!$       open(58,file='matrix_inv_failed.dat')
!!$       do i = 1, n
!!$          do j = 1, n
!!$             write(58,*) i, j, matrix(i,j)
!!$          end do
!!$       end do
!!$       close(58)
!!$       stop
!!$    else
!!$
!!$       call DGETRI(n, matrix, lda, ipiv, work, lwork, info)
!!$
!!$       if (info /= 0) then
!!$          write(*,*) 'DGETRI: Inversion failed. Info = ', info
!!$          stop
!!$       end if
!!$
!!$    end if
!!$
!!$    deallocate(work)
!!$    deallocate(ipiv)
!!$
!!$  end subroutine invert_matrix_dp

  subroutine invert_matrix_sp(matrix)
    implicit none

    real(sp), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(sp),     allocatable, dimension(:)   :: work

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(n))

    call SGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'SGETRF: LU factorization failed. Info = ', info
       stop
    else

       call SGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'SGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_sp

  subroutine compute_pseudo_inverse(A, Across)
    implicit none
    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: Across

    integer(i4b) :: i, m, n, lwork, info
    real(dp)     :: threshold
    real(dp),     allocatable, dimension(:)   :: work, S
    real(dp),     allocatable, dimension(:,:) :: Atmp, U, Vt, Sigma
    
    threshold = 1d-12
    m = size(A,1)
    n = size(A,2)
    lwork = 2*max(1,3*min(m,n)+max(m,n),5*min(m,n))

    ! Compute SVD factorization
    allocate(Atmp(m,n), work(lwork), S(min(m,n)), U(m,m), Vt(n,n), Sigma(n,m))
    Atmp = A
    call DGESVD('A', 'A', m, n, Atmp, m, S, U, m, VT, n, work, lwork, info)
    if (info /= 0) then
       do i = 1, m
          write(*,*) real(A(i,:),sp)
       end do
       write(*,*) 'Error DGESVD = ', info
       stop
    end if

    ! Invert S
    Sigma = 0.d0
    do i = min(m,n), 1, -1
       if (abs(S(i)) > threshold*abs(S(1))) then
          !write(*,*) S(i), threshold, S(1)
          Sigma(i,i) = 1.d0/S(i)
       else
          Sigma(i,i) = 0.d0
       end if
    end do

!!$    ! Multiply S^+ with Vt
!!$    do i = 1, min(m,n)
!!$       Vt(i,:) = S(i)*Vt(i,:)
!!$    end do
!!$
!!$    Across = matmul(transpose(Vt), transpose(U))

    ! Multiply transposr(Vt) with S^+
!!$    Across = 0.d0
!!$    do i = 1, min(m,n)
!!$       do j = 1, min(m,n)
!!$          Across(i,j) = Vt(j,i) * S(j)
!!$       end do
!!$    end do

    Sigma = matmul(transpose(Vt), Sigma)
    Across = matmul(Sigma, transpose(U))

    deallocate(Atmp, work, S, U, Vt, Sigma)

  end subroutine compute_pseudo_inverse

  subroutine compute_pseudo_inverse2(A, Across)
    implicit none
    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: Across

    integer(i4b) :: i, m, n, lwork, info
    logical(lgt), allocatable, dimension(:)   :: mask
    real(dp),     allocatable, dimension(:)   :: tau, work
    real(dp),     allocatable, dimension(:,:) :: Atmp, R
    
    m = size(A,1)
    n = size(A,2)
    lwork = n

    ! Compute QR factorization
    allocate(Atmp(m,n), tau(m), work(lwork), R(n,n), mask(n))
    Atmp = A
    call DGEQRF(m, n, Atmp, m, tau, work, lwork, info)
    if (info /= 0) write(*,*) 'Error DGEQRF = ', info

    ! Extract R
    R = 0.d0
    do i = 1, n
       R(1:i,i) = Atmp(1:i,i)
    end do

    ! Find non-zero elements
    mask = .true.
    do i = 1, n
       if (R(i,i) == 0.d0) then
          mask(i) = .false.
          R(i,:)  = 0.d0
          R(:,i)  = 0.d0
          R(i,i)  = 1.d0
       end if
    end do

    ! Solve for y = R*Across
    Across = transpose(A)
    call DTRTRS('U', 'T', 'N', n, m, R, n, Across, n, info)
    if (info /= 0) write(*,*) 'Error DTRTRS 1 = ', info

    ! Solve for Across
    call DTRTRS('U', 'N', 'N', n, m, R, n, Across, n, info)
    if (info /= 0) write(*,*) 'Error DTRTRS 2 = ', info

    ! Nullify masked rows and columns
    do i = 1, n
       if (.not. mask(i)) then
          Across(i,:) = 0.d0
          Across(:,i) = 0.d0
       end if
    end do

    deallocate(Atmp, tau, work, R, mask)

  end subroutine compute_pseudo_inverse2


  subroutine invert_matrix_with_mask_dpc(matrix)
    implicit none

    complex(dpc), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: i, n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    complex(dpc), allocatable, dimension(:)   :: work
    logical(lgt), allocatable, dimension(:)   :: mask

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(n))
    allocate(mask(n))

    mask = .true.
    do i = 1, n
       if (abs(matrix(i,i)) <= 0.d0) then
          mask(i)     = .false.
          matrix(i,i) = cmplx(1.d0,0.d0,dp)
       end if
    end do

    call ZGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'ZGETRF: LU factorization failed. Info = ', info
       stop
    else

       call ZGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'ZGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    do i = 1, n
       if (.not. mask(i)) then
          matrix(i,i) = cmplx(0.d0,0.d0,dp)
       end if
    end do

    deallocate(mask)
    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_with_mask_dpc
  

  subroutine invert_matrix_with_mask_dp(matrix)
    implicit none

    real(dp), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: i, n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(dp),     allocatable, dimension(:)   :: work
    logical(lgt), allocatable, dimension(:)   :: mask

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(n))
    allocate(mask(n))

    mask = .true.
    do i = 1, n
       if (abs(matrix(i,i)) <= 0.d0) then
          mask(i)     = .false.
          matrix(i,i) = 1.d0
       end if
    end do
    
    call DGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'DGETRF: LU factorization failed. Info = ', info
       stop
    else

       call DGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'DGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    do i = 1, n
       if (.not. mask(i)) then
          matrix(i,i) = 0.d0
       end if
    end do

    deallocate(work)
    deallocate(ipiv)
    deallocate(mask)

  end subroutine invert_matrix_with_mask_dp

  subroutine invert_matrix_with_mask_sp(matrix)
    implicit none

    real(sp), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: i, n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(sp),     allocatable, dimension(:)   :: work
    logical(lgt), allocatable, dimension(:)   :: mask

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(n))
    allocate(mask(n))

    mask = .true.
    do i = 1, n
       if (abs(matrix(i,i)) <= 0.0) then
          mask(i)     = .false.
          matrix(i,i) = 1.0
       end if
    end do

    call SGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'SGETRF: LU factorization failed. Info = ', info
       stop
    else

       call SGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'SGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    do i = 1, n
       if (.not. mask(i)) then
          matrix(i,i) = 0.0
       end if
    end do

    deallocate(work)
    deallocate(ipiv)
    deallocate(mask)

  end subroutine invert_matrix_with_mask_sp



  subroutine get_eigen_decomposition(matrix, eigenvals, eigenvectors)
    
    real(dp), dimension(1:,1:), intent(in)            :: matrix
    real(dp), dimension(1:),    intent(out)           :: eigenvals
    real(dp), dimension(1:,1:), intent(out)           :: eigenvectors

    integer(i8b)     :: i, n, liwork, lwork, lda, ldb, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:)   :: work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    job    = 'v'
    uplo   = 'l'
    n      = size(eigenvals)
    lda    = n
    ldb    = n
    liwork = 5*n + 3
    lwork  = 2*n**2 + 6*n + 1
    info   = 0
      
    ! Perform eigenvalue decomposition
    allocate(work(lwork))
    allocate(iwork(liwork))
    eigenvectors = matrix
    call dsyevd(job, uplo, n, eigenvectors, lda, eigenvals, work, lwork, iwork, liwork, info)
    if (info /= 0) then
       write(*,*) 'get_eigen_decomposition -- dsyevd info = ', info
       do i = 1, n
          write(*,*) real(matrix(i,:),sp)
       end do
       stop
    end if

    deallocate(work)
    deallocate(iwork)

  end subroutine get_eigen_decomposition


  subroutine get_eigenvalues(A, eigenvals)
    
    real(dp), dimension(1:,1:), intent(in)            :: A
    real(dp), dimension(1:),    intent(out)           :: eigenvals

    integer(i4b)     :: n, liwork, lwork, lda, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:,:) :: A_copy
    real(dp),     allocatable, dimension(:)   :: work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    n      = size(eigenvals)

    if (n == 1) then

       eigenvals(1) = A(1,1)

    else if (n == 2) then

       eigenvals(1) = 0.5d0*(A(1,1)+A(2,2) + sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))
       eigenvals(2) = 0.5d0*(A(1,1)+A(2,2) - sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))

    else

       job    = 'n'
       uplo   = 'l'
       lda    = n
       liwork = 1
       lwork  = 2*n+1

       ! Perform eigenvalue decomposition
       allocate(work(lwork), iwork(liwork), A_copy(n,n))
       A_copy = A
       call dsyevd(job, uplo, n, A_copy, lda, eigenvals, work, lwork, iwork, liwork, info)
       if (info /= 0) write(*,*) 'get_eigenvalues -- dsyevd info = ', info
       
       deallocate(work, iwork, A_copy)

    end if

  end subroutine get_eigenvalues

  subroutine eigen_decomp(matrix, eigenvals, eigenvectors, status)
    implicit none
    real(dp)              :: matrix(:,:), eigenvals(:), eigenvectors(:,:)
    real(dp), allocatable :: work(:), iwork(:)
    integer(i4b)          :: n, status
    n = size(matrix,1)
    allocate(work(2*n**2+6*n+1), iwork(5*n+3))
    eigenvectors = matrix
    call dsyevd('v', 'l', n, eigenvectors, n, eigenvals, work, size(work), &
      & iwork, size(iwork), status)
    deallocate(work, iwork)
  end subroutine eigen_decomp

  subroutine compute_hermitian_root(A, pow, trunc)
    implicit none

    real(dp),                   intent(in)    :: pow
    real(dp), dimension(1:,1:), intent(inout) :: A
    logical(lgt), intent(in), optional :: trunc

    integer(i8b)     :: i, n, liwork, lwork, lda, ldb, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:,:) :: V
    real(dp),     allocatable, dimension(:)   :: W, work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    job    = 'v'
    uplo   = 'l'
    n      = size(A(1,:))
    lda    = n
    ldb    = n
    liwork = 5*n + 3
    lwork  = 2*n**2 + 6*n + 1
      
    ! Perform eigenvalue decomposition
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(V(n,n))
    allocate(W(n))
    V = A
    call dsyevd(job, uplo, n, V, lda, W, work, lwork, iwork, liwork, info)

    if (present(trunc)) then
       where (W < 0.d0)
          W = 0.d0
       end where
    else 
       if (any(W <= 0.d0)) then
!       write(*,*) 'W = ', W
          A(1,1) = -1.d30
          deallocate(V)
          deallocate(W)
          deallocate(work)
          deallocate(iwork)
          return
       end if
    end if

    ! Re-compose matrix
    do i = 1, n
       A(i,:) = W(i)**pow * V(:,i)
    end do
    A = matmul(V, A)
    
    deallocate(V)
    deallocate(W)
    deallocate(work)
    deallocate(iwork)

  end subroutine compute_hermitian_root


  subroutine compute_hermitian_root_with_mask(A, pow, A2, pow2)
    implicit none

    real(dp),                   intent(in)    :: pow
    real(dp), dimension(1:,1:), intent(inout) :: A
    real(dp),                   intent(in),    optional :: pow2
    real(dp), dimension(1:,1:), intent(inout), optional :: A2

    integer(i8b)     :: i, n, liwork, lwork, lda, ldb, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:,:) :: V, WVt
    real(dp),     allocatable, dimension(:)   :: W, work
    integer(i4b), allocatable, dimension(:)   :: iwork
    logical(lgt), allocatable, dimension(:)   :: mask

    job    = 'v'
    uplo   = 'l'
    n      = size(A(1,:))
    lda    = n
    ldb    = n
    liwork = 5*n + 3
    lwork  = 2*n**2 + 6*n + 1

    allocate(mask(n))
    do i = 1, n
       mask(i) = (A(i,i) > 0.d0) 
       if (.not. mask(i)) then
          A(i,:) = 0.d0
          A(:,i) = 0.d0
          A(i,i) = 1.d0
       end if
    end do

    ! Perform eigenvalue decomposition
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(V(n,n), WVt(n,n))
    allocate(W(n))
    V = A
    call dsyevd(job, uplo, n, V, lda, W, work, lwork, iwork, liwork, info)

    if (any(W <= 0.d0)) then
!       write(*,*) 'W = ', W
       A(1,1) = -1.d30
       return
    end if

    ! Re-compose matrix
    do i = 1, n
       WVt(i,:) = W(i)**pow * V(:,i)
    end do
    call dgemm('N','N',n,n,n,1.d0,V,n,WVt,n,0.d0,A,n)
    !A = matmul(V, A)
    
    ! Nullify masked elements
    do i = 1, n
       if (.not. mask(i)) then
          A(i,:) = 0.d0
          A(:,i) = 0.d0
       end if
    end do

    if (present(A2)) then
       ! Re-compose matrix
       do i = 1, n
          WVt(i,:) = W(i)**pow2 * V(:,i)
       end do
       call dgemm('N','N',n,n,n,1.d0,V,n,WVt,n,0.d0,A2,n)
       !A2 = matmul(V, A2)

       ! Nullify masked elements
       do i = 1, n
          if (.not. mask(i)) then
             A2(i,:) = 0.d0
             A2(:,i) = 0.d0
          end if
       end do
    end if

    deallocate(V, WVt)
    deallocate(W)
    deallocate(work)
    deallocate(iwork)
    deallocate(mask)

  end subroutine compute_hermitian_root_with_mask
  

  !------------------------------------------------------------------
  ! Subroutines for inverting a matrix. Based on 
  ! Bo Einarssons F90-manual.
  !------------------------------------------------------------------

!!$  subroutine solve_system(A, X, B)
!!$    implicit none
!!$
!!$    complex(dpc), dimension (:, :)               :: A
!!$    complex(dpc), dimension (:)                  :: X
!!$    complex(dpc), dimension (:)                  :: B
!!$
!!$    complex(dpc), dimension(size(B), size(B)+1)  :: m
!!$    integer, dimension (1)                       :: max_loc
!!$    complex(dpc), dimension(size(B)+1)           :: temp_row
!!$    integer                                      :: N, K, I 
!!$    N = size (B)
!!$    m (1:N, 1:N) = A
!!$    m (1:N, N+1) = B 
!!$    
!!$    do K = 1, N - 1
!!$
!!$       max_loc = maxloc (abs (m (K:N, K)))
!!$       if ( max_loc(1) /= 1 ) then
!!$          temp_row (K:N+1 ) = m (K, K:N+1)
!!$          m (K, K:N+1)= m (K-1+max_loc(1), K:N+1)
!!$          m (K-1+max_loc(1), K:N+1) = temp_row( K:N+1)
!!$       end if
!!$
!!$       temp_row (K+1:N) = m (K+1:N, K) / m (K, K)
!!$       do I = K+1, N
!!$          m (I, K+1:N+1) = m (I, K+1:N+1) - &
!!$               temp_row (I) * m (K, K+1:N+1)
!!$       end do
!!$       m (K+1:N, K) = cmplx(0.d0,0.d0)
!!$
!!$    end do 
!!$
!!$    do K = N, 1, -1
!!$       X (K) = ( m (K, N+1) - &
!!$            sum (m (K, K+1:N) * X (K+1:N)) ) / m (K, K)
!!$    end do
!!$
!!$  end subroutine 


  subroutine invert_singular_matrix(matrix, threshold)
    implicit none

    real(dp), dimension(1:,1:), intent(inout)         :: matrix
    real(dp),                   intent(in)            :: threshold
    
    real(dp), allocatable, dimension(:)               :: eigenvals
    real(dp), allocatable, dimension(:,:)             :: eigenvectors, matrix2
    integer(i4b)                                      :: myid, i, n    
    real(dp)                                          :: maxeigenval  

    myid = 1000
    n = size(matrix(1,:))
    allocate(matrix2(n,n))
    allocate(eigenvals(n))
    allocate(eigenvectors(n,n))
    call get_eigen_decomposition(matrix, eigenvals, eigenvectors)
    maxeigenval = maxval(eigenvals)
    do i = 1, n
       if (eigenvals(i) == 0.d0) then 
          cycle
       else if (eigenvals(i) > threshold * maxeigenval .or. threshold ==0.d0) then
          eigenvals(i) = 1.d0/eigenvals(i)
       else 
          eigenvals(i) = 0.d0
       end if
    end do

    matrix = transpose(eigenvectors)
    do i = 1, n
       matrix(i,:) = matrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, matrix, n, 0.d0, matrix2, n)
    matrix = matrix2

    deallocate(eigenvals)
    deallocate(eigenvectors)
    deallocate(matrix2)

  end subroutine invert_singular_matrix


  !------------------------------------------------------------------
  ! Subroutines for inverting a matrix. Based on 
  ! Bo Einarssons F90-manual.
  !------------------------------------------------------------------

  subroutine solve_system_real(A, X, B)
    
    real(dp), dimension(1:,1:), intent(in)  :: A
    real(dp), dimension(1:),    intent(in)  :: B
    real(dp), dimension(1:),    intent(out) :: X

    integer(i4b) :: N, nrhs, lda, ldb, info
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(dp),     allocatable, dimension(:,:) :: b_int, A_int

!    real(dp), dimension(size(B), size(B)+1)  :: m
!    integer, dimension (1)                       :: max_loc
!    real(dp), dimension(size(B)+1)           :: temp_row
!    integer                                      :: K, I 

    N    = size(X)
    nrhs = 1
    lda  = n
    ldb  = n

    allocate(A_int(N,N))
    allocate(b_int(N,1))
    allocate(ipiv(N))

    A_int      = A
    b_int(:,1) = B

    call dgesv(N, nrhs, A_int, lda, ipiv, b_int, ldb, info)
    if (info /= 0) then
       write(*,*) 'Error in solution of real system. Info = ', info
       stop
    end if

    X = b_int(:,1)

    deallocate(ipiv)
    deallocate(A_int)
    deallocate(b_int)

!!$    N = size (B)
!!$    m (1:N, 1:N) = A
!!$    m (1:N, N+1) = B 
!!$    
!!$    do K = 1, N - 1
!!$
!!$       max_loc = maxloc (abs (m (K:N, K)))
!!$       if ( max_loc(1) /= 1 ) then
!!$          temp_row (K:N+1 ) = m (K, K:N+1)
!!$          m (K, K:N+1)= m (K-1+max_loc(1), K:N+1)
!!$          m (K-1+max_loc(1), K:N+1) = temp_row( K:N+1)
!!$       end if
!!$
!!$       temp_row (K+1:N) = m (K+1:N, K) / m (K, K)
!!$       do I = K+1, N
!!$          m (I, K+1:N+1) = m (I, K+1:N+1) - &
!!$               temp_row (I) * m (K, K+1:N+1)
!!$       end do
!!$       m (K+1:N, K) = 0.d0
!!$
!!$    end do 
!!$
!!$    do K = N, 1, -1
!!$       X (K) = ( m (K, N+1) - &
!!$            sum (m (K, K+1:N) * X (K+1:N)) ) / m (K, K)
!!$    end do
!!$
!!$    write(*,*) 'X = ', X
!!$    write(*,*)
!!$    write(*,*) 'hei!'
!!$!    call mpi_finalize(i)
!!$!    stop

  end subroutine 


  
  !  THIS routine returns the value of normalised P_lm(theta) such that
  !  2.*PI*Integral_-1^+1 dx P_lm(x)^2 = 1., where x = cos(theta)
  !  
  !  modified P_lm generating routine from HEALPix, by K.M.G. 4, Sept. 2000
  
  !=======================================================================
  subroutine comp_normalised_Plm(nlmax, m, theta, plm)
    !=======================================================================
    !nlmax (integer) = maximum l
    !theta in radians (double precision)
    !lambda (double precision) = modulus of the complex spherical harmonics
    !  contains lambda(l,m) with l,m in [0,nlmax]
    !  contains lambda(l-m) with l-m in [0,nlmax-m]
    !=======================================================================
    IMPLICIT none
    !
    REAL(dp), PARAMETER:: max_dp  = HUGE(1.0d0)
    REAL(DP), PARAMETER :: PI = 3.141592653589793238462643383279502884197d0
    
    INTEGER(I4B), INTENT(IN)  :: m
    INTEGER(I4B), INTENT(IN)  :: nlmax
    REAL(DP),     INTENT(IN)  :: theta
    
    
    REAL(DP) :: lambda(0:nlmax)
    REAL(DP),     DIMENSION(0:nlmax), INTENT(OUT) :: plm
    
    INTEGER(I4B) :: nmmax
    INTEGER(I4B) l, mm               !, m ...  alm related
    
    REAL(DP) sq4pi_inv
    REAL(DP) cth, sth
    REAL(DP) a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) f2m, fm2, fl2
    
    Character(LEN=7), PARAMETER :: code = 'ALM2MAP'
    
    REAL(DP), PARAMETER :: bignorm = 1.d-20*max_dp
    !=======================================================================
    
    
    !      write(*,*)'   PI   =    ',PI
    
    nmmax = nlmax
    
    LAMBDA = 0.0d0
    
    !     --------------------------------------------
    sq4pi_inv = 1.D0 / SQRT(4.D0*PI)
    
    cth = COS(theta)
    sth = SIN(theta)
    
    plm=0.d0
    
    !      write(*,*)cth,sth
    
    !        lambda_mm tends to go down when m increases (risk of underflow)
    !        lambda_lm tends to go up   when l increases (risk of overflow)
    
    lam_mm = sq4pi_inv * bignorm ! lambda_0_0 * big number --> to avoid underflow
    
    !      do m = 0, nmmax
    
    fm2 = DFLOAT(m) **2
    
    !           ---------- l = m ----------
    par_lm = 1.d0  ! = (-1)^(l+m)
    if (m .ge. 1) then ! lambda_0_0 for m>0
       do mm = 1, m
          f2m = 2.d0 * mm
          lam_mm = - lam_mm * sth * DSQRT( (f2m+1.d0)/ f2m )
       enddo
    endif
    
    lam_lm = lam_mm / bignorm ! actual lambda_mm
    
    LAMBDA(M) = LAM_LM
    
    
    !           ---------- l > m ----------
    lam_0 = 0.d0
    lam_1 = 1.d0 / bignorm    ! small number --> to avoid overflow
    fl2 = DFLOAT(m+1) **2
    a_rec = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
    lam_2 = cth * lam_1 * a_rec
    do l = m+1, nmmax
       par_lm = - par_lm  ! = (-1)^(l+m)
       
       lam_lm = lam_2 * lam_mm ! actual lambda_lm (small and big numbers cancel out)
       
       !            lambda(l,m) = lam_lm
       !            lambda(l-m) = lam_lm
       
       LAMBDA(L) = LAM_LM
       
       lam_0 = lam_1 / a_rec
       lam_1 = lam_2
       fl2 = DFLOAT(l+1) **2
       a_rec = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       lam_2 = (cth * lam_1 - lam_0) * a_rec
    enddo
    
    !      enddo
    
    
    plm = lambda
    
    return
  end subroutine comp_normalised_Plm


  subroutine convert_sigma2frac(frac, sigma)
    implicit none

    real(sp), intent(in)  :: sigma
    real(sp), intent(out) :: frac

    real(dp) :: int_frac, int_sigma

    int_sigma = real(sigma,dp)

    int_frac = corr_erf(int_sigma / sqrt(2.d0))

    frac = (1. + real(int_frac,sp)) / 2.

  end subroutine convert_sigma2frac

  ! Solve the equation erf(x/sqrt(2)) - y = 0 using bisection
  subroutine convert_fract2sigma_dp(sigma, fract)
    implicit none

    real(dp), intent(in)  :: fract
    real(dp), intent(out) :: sigma

    integer(i4b)   :: maxit, j
    real(dp)       :: fr, fl, fm, xl, xr, xm

    maxit = 40

    if (fract > 0.999999426d0) then
       sigma = 5.d0
    else if (fract < 1d-7) then
       sigma = 0.d0
    else

       xl = 0.
       xr = 10.d0
       xm = (xl+xr)/2.d0

       fl = -fract
       fr = corr_erf(xr)-fract

       do j = 1, maxit
          if (abs(xl-xr) < 1d-7) exit

          xm = (xl+xr)/2.d0
          fm = corr_erf(xm)-fract

          if (fm*fl < 0.d0) then
             xr = xm
             fr = fm
          else
             xl = xm
             fl = fm
          end if

       end do

       if (j == maxit) then
          write(*,*) 'ERROR: Too many iterations in the fract2sigma search'
       end if

       sigma = sqrt(2.d0) * xm

    end if

  end subroutine convert_fract2sigma_dp

  ! Solve the equation erf(x/sqrt(2)) - y = 0 using bisection
  subroutine convert_fract2sigma_sp(sigma, fract)
    implicit none

    real(sp), intent(in)  :: fract
    real(sp), intent(out) :: sigma

    integer(i4b)   :: maxit, j
    real(dp)       :: fr, fl, fm, xl, xr, xm

    maxit = 40

    if (fract > 0.999999426) then
       sigma = 5.
    else if (fract < 1e-7) then
       sigma = 0.
    else

       xl = 0.
       xr = 10.
       xm = (xl+xr)/2.

       fl = -fract
       fr = corr_erf(xr)-fract

       do j = 1, maxit
          if (abs(xl-xr) < 1e-7) exit

          xm = (xl+xr)/2.
          fm = corr_erf(xm)-fract

          if (fm*fl < 0.0) then
             xr = xm
             fr = fm
          else
             xl = xm
             fl = fm
          end if

       end do

       if (j == maxit) then
          write(*,*) 'ERROR: Too many iterations in the fract2sigma search'
       end if

       sigma = sqrt(2.0) * real(xm,sp)

    end if

  end subroutine convert_fract2sigma_sp

  ! Computes the error function. Borrowed from Numerical Recipes
  real(dp) function corr_erf(x)
    implicit none

    real(dp), intent(in) :: x

    real(dp) :: ans, z, t

    z = abs(x)
    t=1.d0/(1.d0+0.5d0*z)
    ans = t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(0.37409196d0+t*(0.09678418d0+&
         & t*(-0.18628806d0+t*(0.27886807d0+t*(-1.13520398d0+t*(1.48851587d0+&
         & t*(-0.82215223d0+t*0.17087277d0)))))))))
    if (x >= 0.d0) then
       corr_erf = 1.d0-ans
    else
       corr_erf = ans - 1.d0
    end if

  end function corr_erf

  subroutine compute_inv_Cl_and_determ(cls, inv_Cl, log_det_Cl)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: cls
    real(dp), dimension(1:,1:), intent(out) :: inv_Cl
    real(dp),                   intent(out) :: log_det_Cl

    integer(i4b) :: p

    p = size(cls(1,:))

    if (p == 1) then

       inv_Cl(1,1) = 1.d0 / cls(1,1)
       log_det_Cl  = log(cls(1,1))

    else if (p == 2) then

       inv_Cl(1,1) =  cls(2,2)
       inv_Cl(1,2) = -cls(1,2)
       inv_Cl(2,1) = -cls(2,1)
       inv_Cl(2,2) =  cls(1,1)
       inv_Cl      = inv_Cl / (cls(1,1)*cls(2,2) - cls(2,1)**2)

       log_det_Cl  = log(cls(1,1)*cls(2,2) - cls(2,1)**2)

    else if (p == 3) then

       ! Cholesky decompose Cl
       inv_Cl = 0.d0
       
       inv_Cl(1,1) = sqrt(cls(1,1))
       inv_Cl(2,1) = cls(1,2) / inv_Cl(1,1)
       inv_Cl(3,1) = cls(1,3) / inv_Cl(1,1)
       inv_Cl(2,2) = sqrt(cls(2,2) - inv_Cl(2,1)**2)
!       write(*,*) real(cls(2,3),sp), real(inv_Cl(2,1),sp), real(inv_Cl(3,1),sp), real(inv_Cl(2,2),sp)
       inv_Cl(3,2) = (cls(2,3) - inv_Cl(2,1)*inv_Cl(3,1)) / inv_Cl(2,2)
       inv_Cl(3,3) = sqrt(cls(3,3) - inv_Cl(3,1)**2 - inv_Cl(3,2)**2)
       
       log_det_Cl = 2.d0 * (log(inv_Cl(1,1)) + log(inv_Cl(2,2)) + log(inv_Cl(3,3)))
       
       inv_Cl(1,1) = 1.d0 / inv_Cl(1,1)
       inv_Cl(2,2) = 1.d0 / inv_Cl(2,2)
       inv_Cl(3,3) = 1.d0 / inv_Cl(3,3)
       
       inv_Cl(3,1) = inv_Cl(3,2)*inv_Cl(2,1) * inv_Cl(1,1)*inv_Cl(2,2)*inv_Cl(3,3) - &
            & inv_Cl(3,1) * inv_Cl(1,1) * inv_Cl(3,3)
       inv_Cl(3,2) = -inv_Cl(3,2) * inv_Cl(2,2)*inv_Cl(3,3)
       inv_Cl(2,1) = -inv_Cl(2,1) * inv_Cl(1,1)*inv_Cl(2,2)
       
       inv_Cl(1,1) = inv_Cl(1,1)**2 + inv_Cl(2,1)**2 + inv_Cl(3,1)**2
       inv_Cl(1,2) = inv_Cl(2,1)*inv_Cl(2,2) + inv_Cl(3,1)*inv_Cl(3,2)
       inv_Cl(1,3) = inv_Cl(3,1)*inv_Cl(3,3)
       inv_Cl(2,1) = inv_Cl(2,1)*inv_Cl(2,2) + inv_Cl(3,2)*inv_Cl(3,1)
       inv_Cl(2,2) = inv_Cl(2,2)**2 + inv_Cl(3,2)**2
       inv_Cl(2,3) = inv_Cl(3,2) * inv_Cl(3,3)
       inv_Cl(3,1) = inv_Cl(3,1) * inv_Cl(3,3)
       inv_Cl(3,2) = inv_Cl(3,2) * inv_Cl(3,3)
       inv_Cl(3,3) = inv_Cl(3,3)**2

    end if

  end subroutine compute_inv_Cl_and_determ


  function log_det(cls)
    implicit none

    real(dp), dimension(1:,1:),   intent(in)  :: cls
    real(dp)                                  :: log_det

    integer(i4b) :: p

    real(dp), allocatable, dimension(:,:) :: m

    p = size(cls(1,:))

    if (p == 1) then

       log_det = log(cls(1,1))

    else if (p == 2) then

       log_det = log(cls(1,1)*cls(2,2) - cls(2,1)**2)

    else if (p == 3) then

       allocate(m(p,p))
       
       ! Cholesky decompose Cl
       m(1,1) = sqrt(cls(1,1))
       m(2,1) = cls(1,2) / m(1,1)
       m(3,1) = cls(1,3) / m(1,1)
       m(2,2) = sqrt(cls(2,2) - m(2,1)**2)
       m(3,2) = (cls(2,3) - m(2,1)*m(3,1)) / m(2,2)
       m(3,3) = sqrt(cls(3,3) - m(3,1)**2 - m(3,2)**2)
       
       log_det = 2.d0 * (log(m(1,1)) + log(m(2,2)) + log(m(3,3)))
       
       deallocate(m)

    else
       write(*,*) 'math_tools: Unsupported dimension in log_det', p
       stop
    end if
    
  end function log_det


  function trace_sigma_inv_Cl(sigma, inv_Cl)
    implicit none

    real(dp), dimension(1:,1:), intent(in) :: sigma, inv_Cl
    real(dp)                               :: trace_sigma_inv_Cl
    
    integer(i4b) :: i, j, nmaps

    nmaps = size(sigma(:,1))

    trace_sigma_inv_Cl = 0.d0
    do i = 1, nmaps
       do j = 1, nmaps
          trace_sigma_inv_Cl = trace_sigma_inv_Cl + sigma(i,j) * inv_Cl(j,i)
       end do
    end do

  end function trace_sigma_inv_Cl


!!$  real(dp) function gammln(xx)
!!$    implicit none
!!$
!!$    real(dp), intent(in) :: xx
!!$
!!$    real(dp) :: x, tmp
!!$    real(dp) :: stp = 2.5066282746310005d0
!!$    integer(i4b) :: i
!!$
!!$    real(dp), dimension(6) :: coef
!!$
!!$    coef(1) = 76.18009172947146d0
!!$    coef(2) = -86.50532032941677d0
!!$    coef(3) = 24.01409824083091d0
!!$    coef(4) = -1.231739572450155d0
!!$    coef(5) = 0.001208650973866179d0
!!$    coef(6) = -0.000005395239384953d0
!!$
!!$    x = xx
!!$    tmp = x + 5.5d0
!!$    tmp = (x+0.5d0) * log(tmp) - tmp
!!$    gammln = 0.d0
!!$    do i = 1, 6
!!$       gammln = gammln + coef(i)/(x+real(i,dp))
!!$    end do
!!$    gammln = tmp+log(stp*(gammln + 1.000000000190015d0)/x)
!!$
!!$  end function gammln


  subroutine cholesky_decompose(A, L, ierr)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: L
    integer(i4b),             intent(out), optional :: ierr

    integer(i4b) :: N, i, j, stat

    N = size(A(1,:))

    L = A
    call dpotrf( 'L', N, L, N, stat )

    if (stat /= 0) then

       if (present(ierr)) then
          ierr = stat
       else
          write(*,*) 'Cholesky decomposition failed. stat = ', stat
          stop
       end if
    end if

    if (stat == 0) then
       do i = 1, N
          do j = i+1, N
             L(i,j) = 0.d0
          end do
       end do
    end if

  end subroutine cholesky_decompose


  subroutine cholesky_decompose_single(A, ierr)
    implicit none

    real(dp), dimension(:,:), intent(inout)  :: A
    integer(i4b),             intent(out), optional :: ierr

    integer(i4b) :: N, i, j, stat

    N = size(A(1,:))

    call dpotrf( 'L', N, A, N, stat )

    if (present(ierr)) ierr = stat
    if (stat /= 0 .and. .not. present(ierr)) then
       write(*,*) 'Cholesky decomposition failed. stat = ', stat
       stop
    end if

    if (stat == 0) then
       do i = 1, N
          do j = i+1, N
             A(i,j) = 0.d0
          end do
       end do
    end if

  end subroutine cholesky_decompose_single

  subroutine cholesky_solve(L, b, x)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: L
    real(dp), dimension(:),   intent(in)  :: b
    real(dp), dimension(:),   intent(out) :: x

    integer(i4b) :: N, info
    real(dp), allocatable, dimension(:,:) :: b_int

    N = size(L(1,:))
    allocate(b_int(n,1))
    b_int(:,1) = b
    call dpotrs('L', N, 1, L, N, b_int, N, info )
    x = b_int(:,1)
    deallocate(b_int)

  end subroutine cholesky_solve


  subroutine matmul_symm(A, B, C)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: A, B
    real(dp), dimension(:,:), intent(out) :: C

    integer(i4b) :: n

    n = size(A(:,1))

    call dsymm('l', 'l', n, n, 1.d0, A, n, B, n, 0.d0, C, n )

  end subroutine matmul_symm

  subroutine matmul_gen(A, B, C)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: A, B
    real(dp), dimension(:,:), intent(out) :: C

    integer(i4b) :: n

    n = size(A(:,1))

    call dgemm('N','N',n,n,n,1.d0,A,n,B,n,0.d0,C,n)

  end subroutine matmul_gen

  function mean(array) result(res)
    implicit none
    real(dp) :: array(:), res
    res = sum(array)/size(array)
  end function mean

  function variance(array) result(res)
    implicit none
    real(dp) :: array(:), res
    res = sum((array-mean(array))**2)/(size(array)-1)
  end function


  subroutine tridag(a, b, c, r, u)
    implicit none

    real(dp), dimension(:), intent(in)  :: a, b, c, r
    real(dp), dimension(:), intent(out) :: u

    real(dp)     :: bet
    integer(i4b) :: n, j
    real(dp), dimension(size(b)) :: gam
    
    n   = size(b)
    bet = b(1)

    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1)/bet
       bet    = b(j)-a(j-1)*gam(j)
       u(j)   = (r(j)-a(j-1)*u(j-1))/bet
    end do

    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do

  end subroutine tridag



  subroutine MOV3(a1, b1, c1, a2, b2, c2) 
    implicit none

    real(dp), intent(inout) :: a1, b1, c1, a2, b2, c2

    real(dp) :: temp

    temp = a2
    a2   = a1
    a1   = temp

    temp = b2
    b2   = b1
    b1   = temp

    temp = c2
    c2   = c1
    c1   = temp

  end subroutine MOV3


  subroutine compute_covariance_matrix(data, C,  get_sqrt)
    implicit none

    real(dp),     dimension(1:,1:), intent(in)           :: data
    real(dp),     dimension(1:,1:), intent(out)          :: C
    logical(lgt),                   intent(in), optional :: get_sqrt

    integer(i4b) :: i, j, n, m
    real(dp), allocatable, dimension(:) :: mu

    n = size(data,1)
    m = size(data,2)

    allocate(mu(m))
    do i = 1, m
       mu(i) = mean(data(:,i))
    end do

    do i = 1, m
       do j = i, m
          C(i,j) = sum((data(:,i)-mu(i)) * (data(:,j)-mu(j)))
          C(j,i) = C(i,j)
       end do
    end do
    C = C / (n-1.d0)

    if (present(get_sqrt)) call compute_hermitian_root(C, 0.5d0, .true.)

    deallocate(mu)

  end subroutine compute_covariance_matrix

  subroutine fit_polynomial(x, y, a)
    implicit none

    real(dp), dimension(1:),  intent(in) :: x, y
    real(dp), dimension(0:),  intent(out) :: a

    integer(i4b) :: i, j
    real(dp), dimension(0:size(a)-1)              :: b
    real(dp), dimension(0:size(a)-1, 0:size(a)-1) :: C

    do i = 0, size(a)-1
       b(i) = sum(y * x**i)
       do j = i, size(a)-1
          C(i,j) = sum(x**(i+j))
          C(j,i) = C(i,j)
       end do
    end do
    call solve_system_real(C, a, b)

  end subroutine fit_polynomial

end module math_tools
