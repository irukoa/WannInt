module WannInt_utilities

  use WannInt_kinds, only: wp => dp
  use WannInt_definitions, only: cmplx_0, cmplx_1

  implicit none

  private

  public :: diagonalize

  interface diagonalize
    module procedure :: diagonalize_symmetric
    module procedure :: diagonalize_hermitian
  end interface

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

end module WannInt_utilities
