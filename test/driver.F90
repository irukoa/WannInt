program driver
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
  & select_suite, run_selected, get_argument
  !Tests:
  use Basic_Suite, only: Collect_Basic_Tests
  use Extra_Utilities_Suite, only: Collect_Extra_Utilities_Tests
  use Basic_Interpolation_Suite, only: Collect_Basic_Interpolation_Tests
  use Material_Property_Interpolation_Suite, only: Collect_Material_Property_Interpolation_Tests
  use Tough_Randomized_Suite, only: Collect_Tough_Randomized_Tests

  implicit none
  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [new_testsuite("Basic Tests", Collect_Basic_Tests), &
                new_testsuite("Extra Utilities Tests", Collect_Extra_Utilities_Tests), &
                new_testsuite("Basic Interpolation Tests", Collect_Basic_Interpolation_Tests), &
                new_testsuite("Material Property Interpolation Tests", Collect_Material_Property_Interpolation_Tests), &
                new_testsuite("Stress and Randomized Tests", Collect_Tough_Randomized_Tests)]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
    is = select_suite(testsuites, suite_name)
    if (is > 0 .and. is <= size(testsuites)) then
      if (allocated(test_name)) then
        write (error_unit, fmt) "Suite:", testsuites(is)%name
        call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
        if (stat < 0) then
          error stop 1
        end if
      else
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end if
    else
      write (error_unit, fmt) "Available testsuites"
      do is = 1, size(testsuites)
        write (error_unit, fmt) "-", testsuites(is)%name
      end do
      error stop 1
    end if
  else
    do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  end if

  if (stat > 0) then
    write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

end program driver
