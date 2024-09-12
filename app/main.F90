program main
  use iso_fortran_env, only: output_unit
  use WannInt_definitions
  use WannInt_utilities
  use WannInt

  implicit none

  character(len=10) :: ver = "1.0.1"

  write (unit=output_unit, fmt="(A)") ""
  write (unit=output_unit, fmt="(A)") " ____    __    ____  ___      .__   __. .__   __.  __  .__   __. .___________. "
  write (unit=output_unit, fmt="(A)") " \   \  /  \  /   / /   \     |  \ |  | |  \ |  | |  | |  \ |  | |           | "
  write (unit=output_unit, fmt="(A)") "  \   \/    \/   / /  ^  \    |   \|  | |   \|  | |  | |   \|  | `---|  |----` "
  write (unit=output_unit, fmt="(A)") "   \            / /  /_\  \   |  . `  | |  . `  | |  | |  . `  |     |  |      "
  write (unit=output_unit, fmt="(A)") "    \    /\    / /  _____  \  |  |\   | |  |\   | |  | |  |\   |     |  |      "
  write (unit=output_unit, fmt="(A)") "     \__/  \__/ /__/     \__\ |__| \__| |__| \__| |__| |__| \__|     |__|      "
  write (unit=output_unit, fmt="(A)") ""
  write (unit=output_unit, fmt="(A)") " Wannier Interpolation (WannInt) v"//trim(adjustl(ver))//" built."

end program main
