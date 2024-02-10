program Band_Structure

  use WannInt_kinds, only: wp => dp
  use WannInt_utilities, only: diagonalize
  use WannInt, only: crystal

  !In this example we calculate the band structure of Si in the G - X path.

  type(crystal) :: Si

  integer, parameter :: nk = 100

  complex(wp), allocatable :: H(:, :), R(:, :)
  real(wp), allocatable :: eig(:, :)
  real(wp) :: k(3)
  integer :: out
  integer :: ik, ibnd

  call Si%construct(name="Si", &
                    from_file="./material_data/Si_tb.dat", &
                    fermi_energy=6.3869_wp)

  allocate (H(Si%num_bands(), Si%num_bands()), &
            R(Si%num_bands(), Si%num_bands()), &
            eig(nk, Si%num_bands()))

!$OMP   PARALLEL DO PRIVATE(ik, k, H, R) SHARED(eig)
  do ik = 1, nk

    k = [0.5_wp, 0.0_wp, 0.5_wp]*real(ik - 1, wp)/real(nk - 1, wp)

    H = Si%hamiltonian(kpt=k)
    call diagonalize(matrix=H, P=R, eig=eig(ik, :))

  enddo

  open (newunit=out, action="write", file=trim(Si%name())//"-bands.dat")
  do ibnd = 1, Si%num_bands()
    do ik = 1, nk
      write (unit=out, fmt=*) ik, eig(ik, ibnd)
    enddo
    write (unit=out, fmt=*) ""
  enddo
  close (unit=out)

  deallocate (H, R, eig)

end program Band_Structure
