module Basic_Suite
  use WannInt_kinds, only: wp => dp
  use WannInt_definitions, only: cmplx_0, cmplx_1
  use WannInt_utilities, only: diagonalize, inverse
  use WannInt, only: crystal
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  real(wp) :: tol = 1.0E2_wp

  public :: Collect_Basic_Tests

contains

  subroutine Collect_Basic_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Diagonalization of symmetric matrices (2x2)", test_diagonalization_sym_22), &
                 new_unittest("Diagonalization of Hermitian matrices (2x2)", test_diagonalization_herm_22), &
                 new_unittest("Diagonalization of symmetric matrices (random)", test_diagonalization_sym_rand), &
                 new_unittest("Diagonalization of Hermitian matrices (random)", test_diagonalization_herm_rand), &
                 new_unittest("Inverse of real matrices (random)", test_inv_re_rand), &
                 new_unittest("Inverse of complex matrices (random)", test_inv_c_rand), &
                 new_unittest("Initialization of system and return basic properties from input.", create_crystal_from_input), &
                 new_unittest("Initialization of system and return basic properties from file.", create_crystal_from_file)]

  end subroutine Collect_Basic_Tests

  subroutine test_diagonalization_sym_22(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: mat(2, 2), P(2, 2), D(2, 2)
    real(wp) ::  eig(2), val(2)

    mat(1, :) = [0.0_wp, 1.0_wp]
    mat(2, :) = [1.0_wp, 0.0_wp]

    call diagonalize(matrix=mat, P=P, D=D, eig=eig)

    if (abs(eig(1)) - 1.0_wp > 100*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    call random_seed()
    call random_number(mat)

    mat = mat + transpose(mat)

    call diagonalize(matrix=mat, P=P, D=D, eig=eig)

    associate (a=>mat(1, 1), c=>mat(2, 2), b=>mat(1, 2))

      val(1) = ((a + c) - sqrt((a - c)**2.0_wp + 4.0_wp*b**2.0_wp))/2.0_wp
      val(2) = ((a + c) + sqrt((a - c)**2.0_wp + 4.0_wp*b**2.0_wp))/2.0_wp

    end associate

    if ((abs(eig(1) - val(1)) > tol*epsilon(1.0_wp)) .or. &
        (abs(eig(2) - val(2)) > tol*epsilon(1.0_wp))) then
      allocate (error)
      return
    endif

  end subroutine test_diagonalization_sym_22

  subroutine test_diagonalization_herm_22(error)
    type(error_type), allocatable, intent(out) :: error

    complex(wp) :: mat(2, 2), P(2, 2), D(2, 2)
    real(wp) ::  eig(2), val(2)

    real(wp) :: t1, t2, t3, t4

    mat(1, :) = [cmplx_0, cmplx_1]
    mat(2, :) = [cmplx_1, cmplx_0]

    call diagonalize(matrix=mat, P=P, D=D, eig=eig)

    if (abs(eig(1)) - 1.0_wp > 100*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    call random_seed()
    call random_number(t1)
    call random_number(t2)
    call random_number(t3)
    call random_number(t4)

    mat(1, 1) = t1*cmplx_1
    mat(1, 2) = cmplx(t2, t3, wp)
    mat(2, 1) = conjg(mat(1, 2))
    mat(2, 2) = t4*cmplx_1

    mat = mat + transpose(conjg(mat))

    call diagonalize(matrix=mat, P=P, D=D, eig=eig)

    associate (a=>mat(1, 1), c=>mat(2, 2), b=>mat(1, 2))

      val(1) = (real(a + c, wp) - sqrt(real(a - c, wp)**2.0_wp + 4.0_wp*abs(b)**2.0_wp))/2.0_wp
      val(2) = (real(a + c, wp) + sqrt(real(a - c, wp)**2.0_wp + 4.0_wp*abs(b)**2.0_wp))/2.0_wp

    end associate

    if ((abs(eig(1) - val(1)) > tol*epsilon(1.0_wp)) .or. &
        (abs(eig(2) - val(2)) > tol*epsilon(1.0_wp))) then
      allocate (error)
      return
    endif

  end subroutine test_diagonalization_herm_22

  subroutine test_diagonalization_sym_rand(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp), allocatable :: mat(:, :), P(:, :), D(:, :)

    real(wp) :: t1
    integer :: n, i, j

    call random_seed()
    call random_number(t1)

    n = nint(50.0_wp*t1) + 1

    allocate (mat(n, n), P(n, n), D(n, n))

    do i = 1, n
      do j = 1, i
        call random_number(mat(i, j))
        mat(j, i) = mat(i, j)
      enddo
    enddo

    call diagonalize(matrix=mat, P=P, D=D)

    D = matmul(P, matmul(D, transpose(P)))

    if (abs(D(1, 1) - mat(1, 1)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    deallocate (mat, P, D)

  end subroutine test_diagonalization_sym_rand

  subroutine test_inv_re_rand(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp), allocatable :: mat(:, :), prod(:, :)
    real(wp) :: szr
    real(wp) :: offdiag
    integer :: sz, n, m

    call random_seed()
    call random_number(szr)

    sz = nint(50.0_wp*szr) + 5

    allocate (mat(sz, sz), prod(sz, sz))

    call random_number(mat)
    mat = 10.0_wp*(mat - 0.5_wp)

    !Since strictly diagonally dominant matrices are invertible,
    !(see https://en.wikipedia.org/wiki/Diagonally_dominant_matrix),
    !create a random diagonally dominant matrix.
    do n = 1, sz
      offdiag = 0.0_wp
      do m = 1, sz
        if (n == m) cycle
        offdiag = offdiag + abs(mat(n, m))
      enddo
      mat(n, n) = offdiag + 1.0_wp
    enddo

    prod = matmul(mat, inverse(mat))

    do n = 1, sz
      if (abs(prod(n, n) - 1.0_wp) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    enddo

    deallocate (mat, prod)

  end subroutine test_inv_re_rand

  subroutine test_inv_c_rand(error)
    type(error_type), allocatable, intent(out) :: error

    complex(wp), allocatable :: mat(:, :), prod(:, :)
    real(wp) :: szr, entryr, entryc
    real(wp) :: offdiag
    integer :: sz, n, m

    call random_seed()
    call random_number(szr)

    sz = nint(50.0_wp*szr) + 5

    allocate (mat(sz, sz), prod(sz, sz))

    do n = 1, sz
      do m = 1, sz
        call random_number(entryr)
        call random_number(entryc)
        entryr = 10.0_wp*(entryr - 0.5_wp)
        entryc = 10.0_wp*(entryc - 0.5_wp)
        mat(n, m) = cmplx(entryr, entryc, wp)
      enddo
    enddo

    !Since strictly diagonally dominant matrices are invertible,
    !(see https://en.wikipedia.org/wiki/Diagonally_dominant_matrix),
    !create a random diagonally dominant matrix.
    do n = 1, sz
      offdiag = 0.0_wp
      do m = 1, sz
        if (n == m) cycle
        offdiag = offdiag + abs(mat(n, m))
      enddo
      mat(n, n) = cmplx(offdiag + 1.0_wp, 0.0_wp, wp)
    enddo

    prod = matmul(mat, inverse(mat))

    do n = 1, sz
      if (abs(prod(n, n) - 1.0_wp) > tol*epsilon(1.0_wp)) then
        allocate (error)
        return
      endif
    enddo

    deallocate (mat, prod)

  end subroutine test_inv_c_rand

  subroutine test_diagonalization_herm_rand(error)
    type(error_type), allocatable, intent(out) :: error

    complex(wp), allocatable :: mat(:, :), P(:, :), D(:, :)

    real(wp) :: t1, t2
    integer :: n, i, j

    call random_seed()
    call random_number(t1)

    n = nint(50.0_wp*t1) + 1

    allocate (mat(n, n), P(n, n), D(n, n))

    do i = 1, n
      do j = 1, i
        call random_number(t1)
        if (i == j) then
          mat(i, i) = t1*cmplx_1
        else
          call random_number(t2)
          mat(i, j) = cmplx(t1, t2, wp)
          mat(j, i) = conjg(mat(i, j))
        endif
      enddo
    enddo

    call diagonalize(matrix=mat, P=P, D=D)

    D = matmul(P, matmul(D, transpose(conjg(P))))

    if (abs(abs(D(1, 1)) - abs(mat(1, 1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    deallocate (mat, P, D)

  end subroutine test_diagonalization_herm_rand

  subroutine create_crystal_from_input(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    real(wp) :: dlb(3, 3), rot(3, 3), eig(3), vlm
    integer :: rpts(5, 3)
    complex(wp) :: tunnellings(5, 2, 2), &
                   dipoles(5, 2, 2, 3)

    dlb(1, :) = [1.0_wp, 2.0_wp, 1.5_wp]
    dlb(2, :) = [3.0_wp, 0.5_wp, 5.0_wp]
    dlb(3, :) = [8.0_wp, 0.0_wp, 2.0_wp]

    rpts(1, :) = [-2, 1, 1]
    rpts(2, :) = [-1, 1, 1]
    rpts(3, :) = [-0, 1, 1]
    rpts(4, :) = [1, 1, 1]
    rpts(5, :) = [2, 1, 1]

    tunnellings(3, 1, 2) = sqrt(2.0_wp)
    dipoles(3, 1, 2, 1) = sqrt(3.0_wp)

    if (dummy%initialized()) then
      allocate (error)
      return
    endif

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=2, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=3.0_wp)

    if (dummy%name() /= "dummy") then
      allocate (error)
      return
    endif

    dlb = dummy%direct_lattice_basis()
    if (abs(dlb(1, 1) - 1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    dlb = dummy%reciprocal_lattice_basis()
    if (abs(dlb(2, 1) - (-0.39893240045584766_wp)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    dlb = dummy%metric_tensor()
    call diagonalize(matrix=dlb, P=rot, eig=eig)
    if (abs(dlb(3, 2) - (34.0_wp)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    vlm = dummy%cell_volume()
    if (abs(vlm - sqrt(product(eig))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    if (dummy%nrpts() /= 5) then
      allocate (error)
      return
    endif

    rpts = dummy%rpts()
    if (rpts(1, 1) /= -2) then
      allocate (error)
      return
    endif

    rpts(:, 1) = dummy%deg_rpts()
    if (rpts(1, 1) /= 1) then
      allocate (error)
      return
    endif

    tunnellings = dummy%get_real_space_hamiltonian_elements()
    if (abs(tunnellings(3, 1, 2) - sqrt(2.0_wp)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    dipoles = dummy%get_real_space_position_elements()
    if (abs(dipoles(3, 1, 2, 1) - sqrt(3.0_wp)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    vlm = dummy%fermi_energy()
    if (abs(vlm - 3.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

  end subroutine create_crystal_from_input

  subroutine create_crystal_from_file(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), rot(3, 3), eig(3), vlm
    integer :: rpts(5, 3)
    complex(wp) :: tunnellings(5, 2, 2), &
                   dipoles(5, 2, 2, 3)

    if (dummy%initialized()) then
      allocate (error)
      return
    endif

    call dummy%construct(name="dummy", &
                         from_file="./test/suites/dummy.dat", &
                         fermi_energy=3.0_wp)

    if (dummy%name() /= "dummy") then
      allocate (error)
      return
    endif

    dlb = dummy%direct_lattice_basis()
    if (abs(dlb(1, 1) - 1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    dlb = dummy%reciprocal_lattice_basis()
    if (abs(dlb(2, 1) - (-0.39893240045584766_wp)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    dlb = dummy%metric_tensor()
    call diagonalize(matrix=dlb, P=rot, eig=eig)
    if (abs(dlb(3, 2) - (34.0_wp)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    vlm = dummy%cell_volume()
    if (abs(vlm - sqrt(product(eig))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    if (dummy%nrpts() /= 5) then
      allocate (error)
      return
    endif

    rpts = dummy%rpts()
    if (rpts(1, 1) /= -2) then
      allocate (error)
      return
    endif

    rpts(:, 1) = dummy%deg_rpts()
    if (rpts(1, 1) /= 1) then
      allocate (error)
      return
    endif

    tunnellings = dummy%get_real_space_hamiltonian_elements()
    if (abs(tunnellings(3, 1, 2) - 1.5_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    dipoles = dummy%get_real_space_position_elements()
    if (abs(dipoles(3, 1, 2, 1) - 5.7_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    vlm = dummy%fermi_energy()
    if (abs(vlm - 3.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

  end subroutine create_crystal_from_file

end module Basic_Suite
