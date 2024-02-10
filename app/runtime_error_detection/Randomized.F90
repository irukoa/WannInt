program Randomized
  use iso_fortran_env, only: error_unit
  use MAC, only: container
  use WannInt_kinds, only: wp => dp
  use WannInt, only: crystal

  implicit none

  integer, parameter :: minbnd = 20, maxbnd = 300, &
                        minnr = 500, maxnr = 3000, &
                        minder = 1, maxder = 10

  real(wp) :: rnd, &
              rnbnd, rnr, rnder

  integer :: nbnd, nr, nder

  call random_seed()
  call random_number(rnd)

  if (nint(rnd) == 0) then
    write (error_unit, fmt="(A)") "Testing Hamiltonian:"
    call random_number(rnd)
    rnd = 3.0_wp*rnd - 0.001_wp
    select case (int(rnd))
    case (0)
      write (error_unit, fmt="(A)") "Branch 1:"
      call random_number(rnbnd)
      nbnd = minbnd + nint(real(maxbnd - minbnd, wp)*rnbnd)
      call random_number(rnr)
      nr = minnr + nint(real(maxnr - minnr, wp)*rnr)
      call random_number(rnder)
      nder = minder + nint(real(maxder - minder, wp)*rnder)
      write (error_unit, fmt="(A, i3, A, i4, A, i2, A)") "Parameters: nbnd = ", nbnd, ", nr = ", nr, ", nder = ", nder, "."
      call random_exec_h_b1(nbnd, nr, nder)
    case (1)
      write (error_unit, fmt="(A)") "Branch 2:"
      call random_number(rnbnd)
      nbnd = minbnd + nint(real(maxbnd - minbnd, wp)*rnbnd)
      call random_number(rnr)
      nr = minnr + nint(real(maxnr - minnr, wp)*rnr)
      call random_number(rnder)
      nder = minder + nint(real(maxder - minder, wp)*rnder)
      write (error_unit, fmt="(A, i3, A, i4, A, i2, A)") "Parameters: nbnd = ", nbnd, ", nr = ", nr, ", nder = ", nder, "."
      call random_exec_h_b2(nbnd, nr, nder)
    case (2)
      write (error_unit, fmt="(A)") "Branch 3:"
      call random_number(rnbnd)
      nbnd = minbnd + nint(real(maxbnd - minbnd, wp)*rnbnd)
      call random_number(rnr)
      nr = minnr + nint(real(maxnr - minnr, wp)*rnr)
      write (error_unit, fmt="(A, i3, A, i4, A)") "Parameters: nbnd = ", nbnd, ", nr = ", nr, "."
      call random_exec_h_b3(nbnd, nr)
    end select
  else
    write (error_unit, fmt="(A)") "Testing Berry connection:"
    call random_number(rnd)
    rnd = 3.0_wp*rnd - 0.001_wp
    select case (int(rnd))
    case (0)
      write (error_unit, fmt="(A)") "Branch 1:"
      call random_number(rnbnd)
      nbnd = minbnd + nint(real(maxbnd - minbnd, wp)*rnbnd)
      call random_number(rnr)
      nr = minnr + nint(real(maxnr - minnr, wp)*rnr)
      call random_number(rnder)
      nder = minder + nint(real(maxder - minder, wp)*rnder)
      write (error_unit, fmt="(A, i3, A, i4, A, i2, A)") "Parameters: nbnd = ", nbnd, ", nr = ", nr, ", nder = ", nder, "."
      call random_exec_a_b1(nbnd, nr, nder)
    case (1)
      write (error_unit, fmt="(A)") "Branch 2:"
      call random_number(rnbnd)
      nbnd = minbnd + nint(real(maxbnd - minbnd, wp)*rnbnd)
      call random_number(rnr)
      nr = minnr + nint(real(maxnr - minnr, wp)*rnr)
      call random_number(rnder)
      nder = minder + nint(real(maxder - minder, wp)*rnder)
      write (error_unit, fmt="(A, i3, A, i4, A, i2, A)") "Parameters: nbnd = ", nbnd, ", nr = ", nr, ", nder = ", nder, "."
      call random_exec_a_b2(nbnd, nr, nder)
    case (2)
      write (error_unit, fmt="(A)") "Branch 3:"
      call random_number(rnbnd)
      nbnd = minbnd + nint(real(maxbnd - minbnd, wp)*rnbnd)
      call random_number(rnr)
      nr = minnr + nint(real(maxnr - minnr, wp)*rnr)
      write (error_unit, fmt="(A, i3, A, i4, A)") "Parameters: nbnd = ", nbnd, ", nr = ", nr, "."
      call random_exec_a_b3(nbnd, nr)
    end select
  endif

contains

  subroutine random_exec_h_b1(nbnd, nr, nder)

    integer, intent(in) :: nbnd, nr, nder

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(nr, 3)

    complex(wp) :: tunnellings(nr, nbnd, nbnd), &
                   dipoles(nr, nbnd, nbnd, 3)

    type(container), allocatable :: alloc_container_res(:)

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=nbnd, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    write (error_unit, fmt="(A)") "All variables allocated."
    alloc_container_res = dummy%hamiltonian(kpt=[0.0_wp, 0.0_wp, 0.0_wp], derivative=nder, all=.true.)

  end subroutine random_exec_h_b1

  subroutine random_exec_h_b2(nbnd, nr, nder)

    integer, intent(in) :: nbnd, nr, nder

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(nr, 3)

    complex(wp) :: tunnellings(nr, nbnd, nbnd), &
                   dipoles(nr, nbnd, nbnd, 3)

    type(container) :: contained_res

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=nbnd, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    write (error_unit, fmt="(A)") "All variables allocated."
    contained_res = dummy%hamiltonian(kpt=[0.0_wp, 0.0_wp, 0.0_wp], derivative=nder)

  end subroutine random_exec_h_b2

  subroutine random_exec_h_b3(nbnd, nr)

    integer, intent(in) :: nbnd, nr

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(nr, 3)

    complex(wp) :: tunnellings(nr, nbnd, nbnd), &
                   dipoles(nr, nbnd, nbnd, 3)

    complex(wp) :: resH(nbnd, nbnd)

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=nbnd, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    write (error_unit, fmt="(A)") "All variables allocated."
    resH = dummy%hamiltonian(kpt=[0.0_wp, 0.0_wp, 0.0_wp])

  end subroutine random_exec_h_b3

  subroutine random_exec_a_b1(nbnd, nr, nder)

    integer, intent(in) :: nbnd, nr, nder

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(nr, 3)

    complex(wp) :: tunnellings(nr, nbnd, nbnd), &
                   dipoles(nr, nbnd, nbnd, 3)

    type(container), allocatable :: alloc_container_res(:)

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=nbnd, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    write (error_unit, fmt="(A)") "All variables allocated."
    alloc_container_res = dummy%berry_connection(kpt=[0.0_wp, 0.0_wp, 0.0_wp], derivative=nder, all=.true.)

  end subroutine random_exec_a_b1

  subroutine random_exec_a_b2(nbnd, nr, nder)

    integer, intent(in) :: nbnd, nr, nder

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(nr, 3)

    complex(wp) :: tunnellings(nr, nbnd, nbnd), &
                   dipoles(nr, nbnd, nbnd, 3)

    type(container) :: contained_res

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=nbnd, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    write (error_unit, fmt="(A)") "All variables allocated."
    contained_res = dummy%berry_connection(kpt=[0.0_wp, 0.0_wp, 0.0_wp], derivative=nder)

  end subroutine random_exec_a_b2

  subroutine random_exec_a_b3(nbnd, nr)

    integer, intent(in) :: nbnd, nr

    type(crystal) :: dummy

    real(wp) :: dlb(3, 3), k(3)
    integer :: rpts(nr, 3)

    complex(wp) :: tunnellings(nr, nbnd, nbnd), &
                   dipoles(nr, nbnd, nbnd, 3)

    complex(wp) :: resA(nbnd, nbnd, 3)

    dlb(1, :) = [1.0_wp, 0.0_wp, 0.0_wp]
    dlb(2, :) = [0.0_wp, 1.0_wp, 0.0_wp]
    dlb(3, :) = [0.0_wp, 0.0_wp, 1.0_wp]

    call dummy%construct(name="dummy", &
                         direct_lattice_basis=dlb, &
                         num_bands=nbnd, &
                         R_points=rpts, &
                         tunnellings=tunnellings, &
                         dipoles=dipoles, &
                         fermi_energy=0.0_wp)

    write (error_unit, fmt="(A)") "All variables allocated."
    resA = dummy%berry_connection(kpt=[0.0_wp, 0.0_wp, 0.0_wp])

  end subroutine random_exec_a_b3

end program Randomized
