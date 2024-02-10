module Material_Property_Interpolation_Suite
  use WannInt_kinds, only: wp => dp
  use WannInt_utilities, only: diagonalize
  use WannInt, only: crystal
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  real(wp) :: tol = 1.0E11_wp

  public :: Collect_Material_Property_Interpolation_Tests

contains

  subroutine Collect_Material_Property_Interpolation_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Test band structure of Si (Gamma point)", test_Si_gamma), &
                 new_unittest("Test band structure of GaAs (S point)", test_GaAs_S)]

  end subroutine Collect_Material_Property_Interpolation_Tests

  subroutine test_Si_gamma(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: Si

    complex(wp), allocatable :: H(:, :), R(:, :)
    real(wp), allocatable :: eig(:)

    real(wp) :: bandgap

    call Si%construct(name="Si", &
                      from_file="./material_data/Si_tb.dat", &
                      fermi_energy=6.3869_wp)

    allocate (H(Si%num_bands(), Si%num_bands()), &
              R(Si%num_bands(), Si%num_bands()), &
              eig(Si%num_bands()))

    H = Si%hamiltonian(kpt=[0.0_wp, 0.0_wp, 0.0_wp])

    call diagonalize(matrix=H, P=R, eig=eig)

    bandgap = eig(5) - eig(4)

    if (abs(bandgap - 2.5732061981822474_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    deallocate (H, R, eig)

  end subroutine test_Si_gamma

  subroutine test_GaAs_S(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: GaAs

    complex(wp), allocatable :: H(:, :), R(:, :)
    real(wp), allocatable :: eig(:)
    real(wp) :: bandgap

    call GaAs%construct(name="GaAs", &
                        from_file="./material_data/GaAs_tb.dat", &
                        fermi_energy=7.7414_wp)

    allocate (H(GaAs%num_bands(), GaAs%num_bands()), &
              R(GaAs%num_bands(), GaAs%num_bands()), &
              eig(GaAs%num_bands()))

    H = GaAs%hamiltonian(kpt=[0.5_wp, 0.5_wp, 0.0_wp])

    call diagonalize(matrix=H, P=R, eig=eig)

    bandgap = eig(5) - eig(4)

    if (abs(bandgap - 4.1408237848598901_wp) > tol*epsilon(1.0_wp)) then
      allocate (error)
      return
    endif

    deallocate (H, R, eig)

  end subroutine test_GaAs_S

end module Material_Property_Interpolation_Suite
