module Basic_Interpolation_Suite
  use MAC, only: container
  use WannInt_kinds, only: wp => dp
  use WannInt_definitions, only: cmplx_0, cmplx_1, cmplx_i, pi
  use WannInt, only: crystal
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  real(wp) :: tol = 1.0E1_wp

  public :: Collect_Basic_Interpolation_Tests

contains

  subroutine Collect_Basic_Interpolation_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Test Fourier expansion of trig. functions (fixed kpt)", test_FT_of_trig_fix), &
                 new_unittest("Test Fourier expansion of trig. functions (random kpt)", test_FT_of_trig_rand)]

  end subroutine Collect_Basic_Interpolation_Tests

  subroutine test_FT_of_trig_fix(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(2, 3)

    complex(wp) :: tunnellings(2, 2, 2), &
                   dipoles(2, 2, 2, 3)

    complex(wp) :: resH(2, 2), resA(2, 2, 3), val
    type(container) :: contained_res
    type(container), allocatable :: alloc_container_res(:)

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    rpts(1, :) = [-1, 0, 0]
    rpts(2, :) = [1, 0, 0]

    tunnellings = cmplx_0

    tunnellings(1, 1, 1) = 0.5_wp*cmplx_1
    tunnellings(2, 1, 1) = 0.5_wp*cmplx_1

    tunnellings(1, 1, 2) = 0.5_wp*cmplx_1
    tunnellings(2, 1, 2) = 0.5_wp*cmplx_1

    tunnellings(1, 2, 1) = 0.5_wp*cmplx_i
    tunnellings(2, 2, 1) = -0.5_wp*cmplx_i

    tunnellings(1, 2, 2) = 0.5_wp*cmplx_i
    tunnellings(2, 2, 2) = -0.5_wp*cmplx_i

    dipoles = cmplx_0

    dipoles(1, 1, 1, 1) = 0.5_wp*cmplx_1
    dipoles(2, 1, 1, 1) = 0.5_wp*cmplx_1

    dipoles(1, 1, 2, 1) = 0.5_wp*cmplx_1
    dipoles(2, 1, 2, 1) = 0.5_wp*cmplx_1

    dipoles(1, 2, 1, 1) = 0.5_wp*cmplx_i
    dipoles(2, 2, 1, 1) = -0.5_wp*cmplx_i

    dipoles(1, 2, 2, 1) = 0.5_wp*cmplx_i
    dipoles(2, 2, 2, 1) = -0.5_wp*cmplx_i

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=2, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    !This implements cos(2*pi*k_x) in the first row of H(k) and
    !sin(2*pi*k_x) in the second.

    k = [0.0_wp, 0.0_wp, 0.0_wp] !H(1, 1) = 1, H(2, 2) = 0
    resH = dummy%hamiltonian(kpt=k)
    if (abs(resH(1, 1) - 1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resH(2, 2) - 0.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    k = [0.5_wp, 0.0_wp, 0.0_wp] !H(1, 1) = -1, H(2, 2) = 0
    resH = dummy%hamiltonian(kpt=k)
    if (abs(resH(1, 1) + cmplx_1*1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resH(2, 2)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    k = [0.25_wp, 0.0_wp, 0.0_wp] !H(1, 1) = 0, H(2, 2) = 1
    resH = dummy%hamiltonian(kpt=k)
    if (abs(resH(1, 1)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resH(2, 2) - cmplx_1*1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    k = [0.125_wp, 0.0_wp, 0.0_wp] !H(1, 1) = H(2, 2) = sqrt(2)/2.
    resH = dummy%hamiltonian(kpt=k)
    if (abs(resH(1, 1) - cmplx_1*sqrt(2.0_wp)/2.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resH(2, 2) - cmplx_1*sqrt(2.0_wp)/2.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    !This implements the 5rd derivative of cos(x) and sin(x) evaluated at pi/2.
    contained_res = dummy%hamiltonian(kpt=[0.25_wp, 0.0_wp, 0.0_wp], derivative=5)
    !H(1, 1)''''' = -sin(pi/2).
    val = contained_res%cdp_storage(contained_res%ind([1, 1, 1, 1, 1, 1, 1]))
    if (abs(val + cmplx_1) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    !H(2, 2)''''' = cos(pi/2).
    val = contained_res%cdp_storage(contained_res%ind([2, 2, 1, 1, 1, 1, 1]))
    if (abs(val) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    alloc_container_res = dummy%hamiltonian(kpt=[0.25_wp, 0.0_wp, 0.0_wp], derivative=5, all=.true.)
    !Same thing.
    val = alloc_container_res(6)%cdp_storage(alloc_container_res(6)%ind([1, 1, 1, 1, 1, 1, 1]))
    if (abs(val + cmplx_1) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    !H(2, 2)''''' = cos(pi/2).
    val = alloc_container_res(6)%cdp_storage(alloc_container_res(6)%ind([2, 2, 1, 1, 1, 1, 1]))
    if (abs(val) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    !Repeat same tests for A.

    k = [0.0_wp, 0.0_wp, 0.0_wp]
    resA = dummy%berry_connection(kpt=k)
    if (abs(resA(1, 1, 1) - 1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resA(2, 2, 1) - 0.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    k = [0.5_wp, 0.0_wp, 0.0_wp]
    resA = dummy%berry_connection(kpt=k)
    if (abs(resA(1, 1, 1) + cmplx_1*1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resA(2, 2, 1)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    k = [0.25_wp, 0.0_wp, 0.0_wp]
    resA = dummy%berry_connection(kpt=k)
    if (abs(resA(1, 1, 1)) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resA(2, 2, 1) - cmplx_1*1.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    k = [0.125_wp, 0.0_wp, 0.0_wp]
    resA = dummy%berry_connection(kpt=k)
    if (abs(resA(1, 1, 1) - cmplx_1*sqrt(2.0_wp)/2.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resA(2, 2, 1) - cmplx_1*sqrt(2.0_wp)/2.0_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    !This implements the 5rd derivative of cos(x) and sin(x) evaluated at pi/2.
    contained_res = dummy%berry_connection(kpt=[0.25_wp, 0.0_wp, 0.0_wp], derivative=5)

    val = contained_res%cdp_storage(contained_res%ind([1, 1, 1, 1, 1, 1, 1, 1]))
    if (abs(val + cmplx_1) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    val = contained_res%cdp_storage(contained_res%ind([2, 2, 1, 1, 1, 1, 1, 1]))
    if (abs(val) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    alloc_container_res = dummy%berry_connection(kpt=[0.25_wp, 0.0_wp, 0.0_wp], derivative=5, all=.true.)

    val = alloc_container_res(6)%cdp_storage(alloc_container_res(6)%ind([1, 1, 1, 1, 1, 1, 1, 1]))
    if (abs(val + cmplx_1) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    val = alloc_container_res(6)%cdp_storage(alloc_container_res(6)%ind([2, 2, 1, 1, 1, 1, 1, 1]))
    if (abs(val) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

  end subroutine test_FT_of_trig_fix

  subroutine test_FT_of_trig_rand(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(2, 3)

    complex(wp) :: tunnellings(2, 2, 2), &
                   dipoles(2, 2, 2, 3)

    complex(wp) :: resH(2, 2), resA(2, 2, 3), val
    type(container) :: contained_res

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    rpts(1, :) = [-1, 0, 0]
    rpts(2, :) = [1, 0, 0]

    tunnellings = cmplx_0

    tunnellings(1, 1, 1) = 0.5_wp*cmplx_1
    tunnellings(2, 1, 1) = 0.5_wp*cmplx_1

    tunnellings(1, 1, 2) = 0.5_wp*cmplx_1
    tunnellings(2, 1, 2) = 0.5_wp*cmplx_1

    tunnellings(1, 2, 1) = 0.5_wp*cmplx_i
    tunnellings(2, 2, 1) = -0.5_wp*cmplx_i

    tunnellings(1, 2, 2) = 0.5_wp*cmplx_i
    tunnellings(2, 2, 2) = -0.5_wp*cmplx_i

    dipoles = cmplx_0

    dipoles(1, 1, 1, 1) = 0.5_wp*cmplx_1
    dipoles(2, 1, 1, 1) = 0.5_wp*cmplx_1

    dipoles(1, 1, 2, 1) = 0.5_wp*cmplx_1
    dipoles(2, 1, 2, 1) = 0.5_wp*cmplx_1

    dipoles(1, 2, 1, 1) = 0.5_wp*cmplx_i
    dipoles(2, 2, 1, 1) = -0.5_wp*cmplx_i

    dipoles(1, 2, 2, 1) = 0.5_wp*cmplx_i
    dipoles(2, 2, 2, 1) = -0.5_wp*cmplx_i

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=2, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    call random_seed()
    call random_number(k)

    k = k - 0.5_wp

    resH = dummy%hamiltonian(kpt=k)
    if (abs(resH(1, 1) - cmplx_1*cos(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resH(2, 2) - cmplx_1*sin(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    contained_res = dummy%hamiltonian(kpt=k, derivative=2)
    val = contained_res%cdp_storage(contained_res%ind([1, 1, 1, 1]))
    if (abs(val + cmplx_1*cos(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    val = contained_res%cdp_storage(contained_res%ind([2, 2, 1, 1]))
    if (abs(val + cmplx_1*sin(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    resA = dummy%berry_connection(kpt=k)
    if (abs(resA(1, 1, 1) - cmplx_1*cos(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif
    if (abs(resA(2, 2, 1) - cmplx_1*sin(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    contained_res = dummy%berry_connection(kpt=k, derivative=2)
    val = contained_res%cdp_storage(contained_res%ind([1, 1, 1, 1, 1]))
    if (abs(val + cmplx_1*cos(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    val = contained_res%cdp_storage(contained_res%ind([2, 2, 1, 1, 1]))
    if (abs(val + cmplx_1*sin(2*pi*k(1))) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

  end subroutine test_FT_of_trig_rand

end module Basic_Interpolation_Suite
