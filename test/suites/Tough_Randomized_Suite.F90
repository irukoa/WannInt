module Tough_Randomized_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MAC, only: container
  use WannInt_kinds, only: wp => dp
  use WannInt, only: crystal
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  public :: Collect_Tough_Randomized_Tests

contains

  subroutine Collect_Tough_Randomized_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Stress randomized (timeout after 1 min.)", test_tough_rnd)]

  end subroutine Collect_Tough_Randomized_Tests

  subroutine test_tough_rnd(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    call execute_command_line("./app/runtime_error_detection/script.sh", &
                              wait=.true., exitstat=i)
    write (error_unit, "(A, i4, A)") "Exitstat: ", i, "."
    if (i == 1) allocate (error)
    if (allocated(error)) return
  end subroutine test_tough_rnd

end module Tough_Randomized_Suite
