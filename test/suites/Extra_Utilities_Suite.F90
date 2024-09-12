module Extra_Utilities_Suite
  use WannInt_kinds, only: wp => dp
  use WannInt_definitions, only: cmplx_i, pi
  use WannInt_utilities, only: diagonalize, &
    dirac_delta, deg_list, &
    schur, SVD, expsh, logu
  use WannInt, only: crystal
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  real(wp) :: tol = 1.0E5_wp

  public :: Collect_Extra_Utilities_Tests

contains

  subroutine Collect_Extra_Utilities_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Test dirac delta", test_dirac_delta), &
                 new_unittest("Test degeneracies of eigenvalue lists", test_deg_list), &
                 new_unittest("Test singular value decomposition", test_SVD), &
                 new_unittest("Test matrix exp and log", test_exp_and_log)]

  end subroutine Collect_Extra_Utilities_Tests

  subroutine test_dirac_delta(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: x, smr

    x = 0.01_wp
    smr = 0.1_wp

    if (abs(dirac_delta(x, smr) - 3.9695254747701179_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

  end subroutine test_dirac_delta

  subroutine test_deg_list(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: Si

    complex(wp), allocatable :: H(:, :), R(:, :)
    real(wp), allocatable :: eig(:)
    integer, allocatable :: dg_lts(:)

    real(wp) :: bandgap

    call Si%construct(name="Si", &
                      from_file="./material_data/Si_tb.dat", &
                      fermi_energy=6.3869_wp)

    allocate (H(Si%num_bands(), Si%num_bands()), &
              R(Si%num_bands(), Si%num_bands()), &
              eig(Si%num_bands()), dg_lts(Si%num_bands()))

    H = Si%hamiltonian(kpt=[0.0_wp, 0.0_wp, 0.0_wp])

    call diagonalize(matrix=H, P=R, eig=eig)

    dg_lts = deg_list(eig, 0.0001_wp)
    if (dg_lts(2) /= 3) then
      allocate (error)
      return
    endif

    eig = 0.0_wp
    dg_lts = deg_list(eig, 0.0001_wp)
    if (dg_lts(1) /= Si%num_bands()) then
      allocate (error)
      return
    endif

    deallocate (H, R, eig, dg_lts)

  end subroutine test_deg_list

  subroutine test_SVD(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: mr(3, 4), eig(3)
    complex(wp) :: m(3, 4), U(3, 3), V(4, 4), sigma(3, 4)
    integer :: i, matches

    mr(1, :) = [1.0_wp, 2.0_wp, 3.0_wp, 0.0_wp]
    mr(2, :) = [4.0_wp, 1.0_wp, 3.0_wp, 0.0_wp]
    mr(3, :) = 2*mr(1, :) + mr(2, :)

    m = cmplx(mr, 0.0_wp, wp)

    call SVD(m, U, V, sigma, eig)

    matches = 0
    do i = 1, 3
      if (abs(eig(i)) < tol*epsilon(1.0_wp)) matches = matches + 1
    enddo
    if (matches /= 1) then
      allocate (error)
      return
    endif

  end subroutine test_SVD

  subroutine test_exp_and_log(error)
    type(error_type), allocatable, intent(out) :: error

    integer :: n, i, matches
    real(wp) :: nr, sv

    real(wp), allocatable :: rnd1(:, :), rnd2(:, :), eig(:)
    complex(wp), allocatable :: m(:, :), res(:, :), p(:, :)

    call random_seed()
    call random_number(nr)

    n = nint(2.0_wp + 28.0_wp*nr)

    allocate (rnd1(n, n), rnd2(n, n), eig(n), m(n, n), res(n, n), p(n, n))

    call random_number(rnd1)
    rnd1 = 1.0_wp + 9.0_wp*rnd1
    call random_number(rnd2)
    rnd2 = 1.0_wp + 9.0_wp*rnd2

    m = cmplx(rnd1, rnd2, wp)
    m = (m - conjg(transpose(m)))/2.0_wp
    call diagonalize(matrix=cmplx_i*m, P=p, eig=eig)
    sv = eig(1)

    res = logu(matrix=expsh(matrix=m))
    call diagonalize(matrix=cmplx_i*res, P=p, eig=eig)

    matches = 0
    do i = 1, n
      if (abs(aimag(log(exp(cmplx_i*sv))) - eig(i)) < 1.0E-5_wp) matches = matches + 1
    enddo
    if (matches < 1) then
      allocate (error)
      return
    endif

    deallocate (rnd1, rnd2, eig, m, res, p)

  end subroutine test_exp_and_log

end module Extra_Utilities_Suite
