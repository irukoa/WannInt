module WannInt_utilities

  use WannInt_kinds, only: wp => dp
  use WannInt_definitions, only: cmplx_0, cmplx_1, cmplx_i, &
    pi

  implicit none

  private

  public :: diagonalize
  public :: inverse

  interface diagonalize
    module procedure :: diagonalize_symmetric
    module procedure :: diagonalize_hermitian
  end interface

  interface inverse
    module procedure :: inverse_real
    module procedure :: inverse_complex
  end interface

  public :: dirac_delta
  public :: deg_list
  public :: schur
  public :: SVD
  public :: expsh
  public :: logu

contains

  subroutine diagonalize_symmetric(matrix, P, D, eig)
    !===========================================!
    !                                           !
    !  Given a symmetric matrix, computes the   !
    !  elements for its diagonalization P, D:   !
    !  matrix = P*D*P^T, where                  !
    !  D_nm = delta_nm eig_n for eig_n real.    !
    !                                           !
    !===========================================!

    real(wp), intent(in)  :: matrix(:, :)
    real(wp), intent(out) :: P(size(matrix(:, 1)), size(matrix(1, :)))

    real(wp), optional, intent(out) :: D(size(matrix(:, 1)), size(matrix(1, :)))
    real(wp), optional, intent(out) :: eig(size(matrix(:, 1)))

    real(wp), allocatable :: work(:)
    real(wp)              :: lst(size(matrix(:, 1)))
    integer               :: dim, info, lwork, i
    external              :: dsyev

    character(len=1024) :: errormsg

    if (size(matrix(:, 1)) /= size(matrix(1, :))) error stop &
      "WannInt: Error #1: matrix to diagonalize is not square."
    dim = size(matrix(:, 1))

    !Initialization.
    P = matrix

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call dsyev('V', 'U', dim, P, dim, lst, work, lwork, info)
    lwork = nint(real(work(1), wp))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call dsyev('V', 'U', dim, P, dim, lst, work, lwork, info)
    deallocate (work)

    if (info /= 0) then
      write (errormsg, "(i20)") info
      errormsg = "WannInt: Error #2: subroutine dsyev failed with info = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    if (present(D)) then
      D = 0.0_wp
      do i = 1, dim
        D(i, i) = lst(i)
      enddo
    endif

    if (present(eig)) eig = lst

  end subroutine diagonalize_symmetric

  subroutine diagonalize_hermitian(matrix, P, D, eig)
    !===========================================!
    !                                           !
    !  Given a Hermitian matrix, computes the   !
    !  elements for its diagonalization P, D:   !
    !  matrix = P*D*P^dagger, where             !
    !  D_nm = delta_nm eig_n for eig_n real.    !
    !                                           !
    !===========================================!

    complex(wp), intent(in)  :: matrix(:, :)
    complex(wp), intent(out) :: P(size(matrix(:, 1)), size(matrix(1, :)))

    complex(wp), optional, intent(out) :: D(size(matrix(:, 1)), size(matrix(1, :)))
    real(wp), optional, intent(out) :: eig(size(matrix(:, 1)))

    complex(wp), allocatable :: work(:)
    real(wp)                 :: lst(size(matrix(:, 1))), &
                                rwork(3*size(matrix(:, 1)) - 2)
    integer                  :: dim, info, lwork, i
    external                 :: zheev

    character(len=1024) :: errormsg

    if (size(matrix(:, 1)) /= size(matrix(1, :))) error stop &
      "WannInt: Error #1: matrix to diagonalize is not square."
    dim = size(matrix(:, 1))

    !Initialization.
    P = matrix

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call zheev('V', 'U', dim, P, dim, lst, work, lwork, rwork, info)
    lwork = nint(real(work(1), wp))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call zheev('V', 'U', dim, P, dim, lst, work, lwork, rwork, info)
    deallocate (work)

    if (info /= 0) then
      write (errormsg, "(i20)") info
      errormsg = "WannInt: Error #2: subroutine zheev failed with info = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    if (present(D)) then
      D = cmplx_0
      do i = 1, dim
        D(i, i) = cmplx(lst(i), 0.0_wp, wp)
      enddo
    endif

    if (present(eig)) eig = lst

  end subroutine diagonalize_hermitian

  pure elemental function dirac_delta(x, smr)

    !Approximation of the Dirac delta of x.

    real(wp), intent(in) :: x, smr
    real(wp) :: arg
    real(wp) :: dirac_delta

    arg = min(200.0_wp, 0.5_wp*(x/smr)**2)
    dirac_delta = exp(-arg)/(smr*sqrt(2*pi))

  end function dirac_delta

  pure function deg_list(eig, degen_thr)
    !========================================================================!
    !Auxiliary routine to get the degree of degeneracy for a given list eig. !
    !The list is supposed to have it's elements stored in ascending order:   !
    !eig(i + 1) >= eig(i).                                                   !
    !The degeneracy threshold is given by degen_thr.                         !
    !The results are set up such that if deg(i) = N,                         !
    !then eig(i) = eig(i + 1) = ... = eig(i + N - 1)                         !
    !with the equality holding up to degen_thr.                              !
    !If i < j < i + N - 1, then deg(j) = 0.                                  !
    !If the value is nondegenerate, then deg(j) = 1.                         !
    !========================================================================!

    implicit none

    real(wp), intent(in) :: degen_thr
    real(wp), intent(in) :: eig(:)

    integer :: deg_list(size(eig))

    integer :: i, j, dim

    if (degen_thr < 0.0_wp) error stop &
      "WannInt: Error #3: 'degen_thr' must be a positive real."

    deg_list = 0
    dim = size(eig)

    do i = 1, dim
      do j = i, dim !In ascending order,
        if (abs(eig(j) - eig(i)) <= degen_thr) then
          !count number of elements equal to eig(i).
          deg_list(i) = deg_list(i) + 1
        endif
      enddo
    enddo

    !In the case of eig(j+2) - eig(j+1) < degen_thr and eig(j+1) - eig(j) < degen_thr,
    !but eig(j+2) - eig(j) > degen_thr (closely packed levels),
    do i = dim - 1, 2, -1
      if ((deg_list(i) > 1) .and. (deg_list(i - 1) > 1)) then
        deg_list(i - 1) = deg_list(i) + 1   !increase the degeneracy value of the
        !degenerate level according to the degenerate levels following it.
        deg_list(i + 1) = 0            !Set the next levels to 0.
      endif
    enddo

    do i = dim - 1, 1, -1
      !At this point, the second index of a closely packed level shall be corrected.
      if ((deg_list(i) /= 1) .and. (deg_list(i) >= deg_list(i + 1)) .and. (deg_list(i + 1) >= 1)) then
        deg_list(i + 1) = 0
      endif
    enddo

  end function deg_list

  subroutine schur(matrix, Z, S, T)
    !==============================================!
    !                                              !
    !  Given a square complex matrix, computes the !
    !  elements for its Schur decomposition Z, S:  !
    !  matrix = Z*S*Z^dagger. T are the diagonal   !
    !  elements of S.                              !
    !                                              !
    !==============================================!

    complex(wp), intent(in)  :: matrix(:, :)
    complex(wp), intent(out) :: Z(size(matrix(:, 1)), size(matrix(1, :)))

    complex(wp), optional, intent(out) :: S(size(matrix(:, 1)), size(matrix(1, :)))
    complex(wp), optional, intent(out) :: T(size(matrix(:, 1)))

    complex(wp)              :: B(size(matrix(:, 1)), size(matrix(1, :)))
    complex(wp)              :: eig(size(matrix(:, 1)))
    complex(wp), allocatable :: work(:)
    real(wp)                 :: rwork(size(matrix(:, 1)), size(matrix(1, :)))
    integer                  :: info, lwork, dim, sdim, i
    logical                  :: bwork(size(matrix(:, 1)), size(matrix(1, :))), select
    external                 :: zgees

    character(len=1024) :: errormsg

    if (size(matrix(:, 1)) /= size(matrix(1, :))) error stop &
      "WannInt: Error #1: matrix to diagonalize is not square."
    dim = size(matrix(:, 1))

    !Initialization.
    B = matrix

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call zgees('V', 'N', select, dim, B, dim, sdim, eig, Z, dim, work, lwork, rwork, bwork, info)
    lwork = nint(real(work(1), wp))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call zgees('V', 'N', select, dim, B, dim, sdim, eig, Z, dim, work, lwork, rwork, bwork, info)
    if (present(S)) S = B
    deallocate (work)

    if (info /= 0) then
      write (errormsg, "(i20)") info
      errormsg = "WannInt: Error #2: subroutine zgees failed with info = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    if (present(T)) T = eig

    if (present(S)) then
      S = cmplx_0
      do i = 1, dim
        S(i, i) = eig(i)
      enddo
    endif

  end subroutine schur

  subroutine SVD(matrix, U, V, sigma, eig)
    !==============================================!
    !                                              !
    !  Given a complex matrix, computes the        !
    !  elements for its SVD decomposition U, V:    !
    !  matrix = U*sigma*V^dagger.                  !
    !  eig are the diagonal elements of sigma.     !
    !                                              !
    !==============================================!

    complex(wp), intent(in)  :: matrix(:, :)
    complex(wp), intent(out) :: U(size(matrix(:, 1)), size(matrix(:, 1))), &
                                V(size(matrix(1, :)), size(matrix(1, :)))

    complex(wp), optional, intent(out) :: sigma(size(matrix(:, 1)), size(matrix(1, :)))
    real(wp), optional, intent(out)    :: eig(min(size(matrix(:, 1)), size(matrix(1, :))))

    complex(wp)              :: B(size(matrix(:, 1)), size(matrix(1, :)))
    real(wp)                 :: sigmaw(min(size(matrix(:, 1)), size(matrix(1, :))))
    complex(wp), allocatable :: work(:)
    integer                  :: m, n, info, lwork, i
    real(wp)                 :: rwork(5*min(size(matrix(1, :)), size(matrix(:, 1))))
    external                 :: zgesvd

    character(len=1024) :: errormsg

    m = size(matrix(:, 1)) !Number of rows of the input matrix.
    n = size(matrix(1, :)) !Number of cols of the input matrix.

    !Initialization.
    B = matrix

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call zgesvd('A', 'A', m, n, B, m, sigmaw, U, m, V, n, work, lwork, rwork, info)
    lwork = nint(real(work(1), wp))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call zgesvd('A', 'A', m, n, B, m, sigmaw, U, m, V, n, work, lwork, rwork, info)
    V = conjg(transpose(V))
    deallocate (work)

    if (info /= 0) then
      write (errormsg, "(i20)") info
      errormsg = "WannInt: Error #2: subroutine zgesvd failed with info = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    if (present(sigma)) then
      sigma = cmplx_0
      do i = 1, size(sigmaw)
        sigma(i, i) = cmplx(sigmaw(i), 0.0_wp, wp)
      enddo
    endif
    if (present(eig)) eig = sigmaw

  end subroutine SVD

  function expsh(matrix)
    !=======================================!
    !                                       !
    !  Given a skew Hermitian matrix,       !
    !  matrix^dagger = - matrix             !
    !  computes the unitary matrix          !
    !  exponential exp(matrix).             !
    !                                       !
    !=======================================!

    complex(wp), intent(in)  :: matrix(:, :)

    complex(wp) :: expsh(size(matrix(:, 1)), size(matrix(1, :)))

    complex(wp) :: P(size(matrix(:, 1)), size(matrix(1, :))), &
                   D(size(matrix(:, 1)), size(matrix(1, :)))

    integer :: i

    call diagonalize(matrix=-cmplx_i*matrix, P=P, D=D)

    do i = 1, size(matrix(:, 1))
      D(i, i) = exp(cmplx_i*D(i, i))
    enddo

    expsh = matmul(P, matmul(D, conjg(transpose(P))))

  end function expsh

  function logu(matrix)
    !=========================================!
    !                                         !
    !  Given a unitary matrix,                !
    !  matrix^dagger = matrix^-1              !
    !  computes the skew Hermitian matrix     !
    !  logarithm log(matrix).                 !
    !                                         !
    !=========================================!

    complex(wp), intent(in)  :: matrix(:, :)

    complex(wp) :: logu(size(matrix(:, 1)), size(matrix(1, :)))

    complex(wp) :: Z(size(matrix(:, 1)), size(matrix(1, :))), &
                   S(size(matrix(:, 1)), size(matrix(1, :)))

    integer :: i

    call schur(matrix=matrix, Z=Z, S=S)

    do i = 1, size(matrix(:, 1))
      S(i, i) = log(S(i, i))
    enddo

    logu = matmul(Z, matmul(S, conjg(transpose(Z))))

  end function logu

  function inverse_real(matrix) result(inv_r)
    !=========================================!
    !                                         !
    !  Given a general real matrix,           !
    !  try to compute its inverse.            !
    !                                         !
    !=========================================!

    real(wp), intent(in)  :: matrix(:, :)
    real(wp) :: inv_r(size(matrix(:, 1)), size(matrix(1, :)))

    complex(wp) :: matrix_c(size(matrix(:, 1)), size(matrix(1, :))), &
                   U(size(matrix(:, 1)), size(matrix(1, :))), &
                   V(size(matrix(:, 1)), size(matrix(1, :))), &
                   S(size(matrix(:, 1)), size(matrix(1, :)))

    integer :: i

    character(len=1024) :: errormsg

    if (size(matrix(:, 1)) /= size(matrix(1, :))) error stop &
      "WannInt: Error #1: matrix to invert is not square."

    matrix_c = cmplx(matrix, 0.0_wp, wp)

    call SVD(matrix=matrix_c, U=U, V=V, sigma=S)

    do i = 1, size(matrix(:, 1))
      if (abs(real(S(i, i), wp)) < 1.0E-8_wp) then
        write (errormsg, "(i20)") i
        errormsg = "WannInt: Error #2: the eigenvalue #"//trim(adjustl(errormsg))//" is 0, cannot invert."
        error stop trim(errormsg)
      endif
      S(i, i) = cmplx(1.0_wp/real(S(i, i), wp), 0.0_wp, wp)
    enddo

    inv_r = real(matmul(matmul(V, S), conjg(transpose(U))), wp)

  end function inverse_real

  function inverse_complex(matrix) result(inv_c)
    !=========================================!
    !                                         !
    !  Given a general real matrix,           !
    !  try to compute its inverse.            !
    !                                         !
    !=========================================!

    complex(wp), intent(in)  :: matrix(:, :)
    complex(wp) :: inv_c(size(matrix(:, 1)), size(matrix(1, :)))

    complex(wp) :: U(size(matrix(:, 1)), size(matrix(1, :))), &
                   V(size(matrix(:, 1)), size(matrix(1, :))), &
                   S(size(matrix(:, 1)), size(matrix(1, :)))

    integer :: i

    character(len=1024) :: errormsg

    if (size(matrix(:, 1)) /= size(matrix(1, :))) error stop &
      "WannInt: Error #1: matrix to invert is not square."

    call SVD(matrix=matrix, U=U, V=V, sigma=S)

    do i = 1, size(matrix(:, 1))
      if (abs(real(S(i, i), wp)) < 1.0E-8_wp) then
        write (errormsg, "(i20)") i
        errormsg = "WannInt: Error #2: the eigenvalue #"//trim(adjustl(errormsg))//" is 0, cannot invert."
        error stop trim(errormsg)
      endif
      S(i, i) = cmplx(1.0_wp/real(S(i, i), wp), 0.0_wp, wp)
    enddo

    inv_c = matmul(matmul(V, S), conjg(transpose(U)))

  end function inverse_complex

end module WannInt_utilities
