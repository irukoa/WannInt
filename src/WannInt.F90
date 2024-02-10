#include "macros.inc"
module WannInt

  use OMP_LIB

  use WannInt_kinds, only: wp => dp
  use WannInt_definitions, only: cmplx_0, cmplx_i, pi
  use WannInt_utilities, only: diagonalize
  use MAC, only: container_specifier, container

  implicit none

  private

  public :: diagonalize

  public :: container_specifier, container

  type, public :: crystal
    private
    character(len=120) :: nm
    !Crystallographyc data.
    real(wp) :: direct_l_b(3, 3), &
                reciprocal_l_b(3, 3)
    real(wp) :: m_tensor(3, 3)
    real(wp) :: c_volume
    !Electronic structure data.
    integer                  :: bands
    integer, allocatable     :: R_points(:, :)
    integer                  :: num_R_points
    integer, allocatable     :: deg_R_points(:)
    complex(wp), allocatable :: real_space_hamiltonian_elements(:, :, :)
    complex(wp), allocatable :: real_space_position_elements(:, :, :, :)
    real(wp) :: e_fermi = 0.0_wp
    !
    logical :: is_initialized = .false.
  contains
    private
    procedure, pass(self) :: construct_from_input
    procedure, pass(self) :: construct_from_file
    generic, public :: construct => construct_from_input, construct_from_file
    procedure, public, pass(self) :: name
    procedure, public, pass(self) :: direct_lattice_basis
    procedure, public, pass(self) :: reciprocal_lattice_basis
    procedure, public, pass(self) :: metric_tensor
    procedure, public, pass(self) :: cell_volume
    procedure, public, pass(self) :: num_bands
    procedure, public, pass(self) :: nrpts
    procedure, public, pass(self) :: rpts
    procedure, public, pass(self) :: deg_rpts
    procedure, public, pass(self) :: get_real_space_hamiltonian_elements
    procedure, public, pass(self) :: get_real_space_position_elements
    procedure, public, pass(self) :: fermi_energy
    procedure, public, pass(self) :: initialized
    procedure, pass(self) :: h_branch_1, h_branch_2, h_branch_3, &
      a_branch_1, a_branch_2, a_branch_3
    generic, public :: hamiltonian => h_branch_1, h_branch_2, h_branch_3
    generic, public :: berry_connection => a_branch_1, a_branch_2, a_branch_3
  end type

contains

  subroutine construct_from_input(self, name, &
                                  direct_lattice_basis, &
                                  num_bands, &
                                  R_points, &
                                  tunnellings, &
                                  dipoles, &
                                  fermi_energy)

    class(crystal), intent(out) :: self

    character(len=*), intent(in) :: name

    real(wp), intent(in) :: direct_lattice_basis(3, 3)
    integer, intent(in) :: num_bands
    integer, intent(in) :: R_points(:, :)
    complex(wp), intent(in) :: tunnellings(:, :, :)
    complex(wp), intent(in) :: dipoles(:, :, :, :)
    real(wp), optional, intent(in) :: fermi_energy

    integer :: i, j
    real(wp) :: P(3, 3), D(3, 3)

    character(len=1024) :: errormsg
    integer :: istat

    if (num_bands < 1) error stop &
      "WannInt: Error #3: 'num_bands' must be greater than 0."
    if (size(R_points(1, :)) /= 3) error stop &
      "WannInt: Error #3: Bravais lattice vector dimension is not 3."
    if (size(tunnellings(:, 1, 1)) /= size(R_points(:, 1))) error stop &
    "WannInt: Error #3: size of 1st index of tunnellings does not equal the &
    &number of points in 'R_points'."
    if (size(tunnellings(1, :, 1)) /= num_bands) error stop &
      "WannInt: Error #3: size of 2nd index of tunnellings is not 'num_bands'."
    if (size(tunnellings(1, 1, :)) /= num_bands) error stop &
      "WannInt: Error #3: size of 3rd index of tunnellings is not 'num_bands'."
    if (size(dipoles(:, 1, 1, 1)) /= size(R_points(:, 1))) error stop &
    "WannInt: Error #3: size of 1st index of dipoles does not equal the &
    &number of points in 'R_points'."
    if (size(dipoles(1, :, 1, 1)) /= num_bands) error stop &
      "WannInt: Error #3: size of 2nd index of dipoles is not 'num_bands'."
    if (size(dipoles(1, 1, :, 1)) /= num_bands) error stop &
      "WannInt: Error #3: size of 3rd index of dipoles is not 'num_bands'."
    if (size(dipoles(1, 1, 1, :)) /= 3) error stop &
      "WannInt: Error #3: size of 4th index of dipoles is not 3."

    self%nm = trim(name)

    !Crystallographyc data block.

    self%direct_l_b = direct_lattice_basis

    do i = 1, 3
      do j = 1, 3
        self%m_tensor(i, j) = dot_product(self%direct_l_b(i, :), self%direct_l_b(j, :))
      enddo
    enddo

    !Get diagonalization of metric tensor.
    call diagonalize(matrix=self%m_tensor, &
                     P=P, D=D)

    !Invert metric tensor's eigenvalues. Also get cell volume from eigenvalues.
    self%c_volume = 1.0_wp
    do i = 1, 3
      if (abs(D(i, i)) < 1.0E-6_wp) error stop &
        "WannInt: Error #4: an eigenvalue of the metric tensor is zero."
      self%c_volume = self%c_volume*D(i, i)
      D(i, i) = 1.0_wp/D(i, i)
    enddo
    self%c_volume = sqrt(self%c_volume)

    !Get inverse metric tensor.
    D = 2*pi*matmul(matmul(P, D), transpose(P))
    !Get reciprocal lattice.
    self%reciprocal_l_b = matmul(D, self%direct_l_b)

    !Electronic structure data block.

    self%bands = num_bands

    self%num_R_points = size(R_points(:, 1))

    allocate (self%R_points(self%num_R_points, 3), &
              self%deg_R_points(self%num_R_points), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "WannInt: Error #5: failure allocating R_points and deg_R_points. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%R_points = R_points
    self%deg_R_points = 1

    allocate (self%real_space_hamiltonian_elements(self%num_R_points, self%bands, self%bands), &
              self%real_space_position_elements(self%num_R_points, self%bands, self%bands, 3), &
              stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "WannInt: Error #5: failure allocating real_space_hamiltonian_elements &
                  &and real_space_position_elements. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%real_space_hamiltonian_elements = tunnellings
    self%real_space_position_elements = dipoles

    if (present(fermi_energy)) self%e_fermi = fermi_energy

    self%is_initialized = .true.

  end subroutine construct_from_input

  subroutine construct_from_file(self, name, &
                                 from_file, &
                                 fermi_energy)

    class(crystal), intent(out) :: self

    character(len=*), intent(in) :: name

    character(len=*), intent(in) :: from_file
    real(wp), optional, intent(in) :: fermi_energy

    integer :: stdin
    integer :: i, j, irpts, &
               division, remainder
    real(wp) :: P(3, 3), D(3, 3)

    character(len=1024) :: errormsg
    integer :: istat

    integer :: dummy1, dummy2, dummy3
    real(wp), allocatable :: dummyR(:)

    self%nm = trim(name)

    !Will read Wannier90 formatted file.
    open (newunit=stdin, action="read", file=trim(from_file), status="unknown")

    read (unit=stdin, fmt=*)

    !Crystallographyc data block.

    do i = 1, 3
      read (unit=stdin, fmt=*) (self%direct_l_b(i, j), j=1, 3)
    enddo

    do i = 1, 3
      do j = 1, 3
        self%m_tensor(i, j) = dot_product(self%direct_l_b(i, :), self%direct_l_b(j, :))
      enddo
    enddo

    !Get diagonalization of metric tensor.
    call diagonalize(matrix=self%m_tensor, &
                     P=P, D=D)

    !Invert metric tensor's eigenvalues. Also get cell volume from eigenvalues.
    self%c_volume = 1.0_wp
    do i = 1, 3
      if (abs(D(i, i)) < 1.0E-6_wp) error stop &
        "WannInt: Error #4: an eigenvalue of the metric tensor is zero."
      self%c_volume = self%c_volume*D(i, i)
      D(i, i) = 1.0_wp/D(i, i)
    enddo
    self%c_volume = sqrt(self%c_volume)

    !Get inverse metric tensor.
    D = 2*pi*matmul(matmul(P, D), transpose(P))
    !Get reciprocal lattice.
    self%reciprocal_l_b = matmul(D, self%direct_l_b)

    !Electronic structure data block.

    read (unit=stdin, fmt=*) self%bands

    read (unit=stdin, fmt=*) self%num_R_points

    allocate (self%R_points(self%num_R_points, 3), &
              self%deg_R_points(self%num_R_points), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "WannInt: Error #5: failure allocating R_points and deg_R_points. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    division = self%num_R_points/15
    remainder = modulo(self%num_R_points, 15)

    if (remainder == 0) then
      do i = 1, division
        read (unit=stdin, fmt=*) (self%deg_R_points(15*(i - 1) + j), j=1, 15)
      enddo
    else
      do i = 1, division
        read (unit=stdin, fmt=*) (self%deg_R_points(15*(i - 1) + j), j=1, 15)
      enddo
      read (unit=stdin, fmt=*) (self%deg_R_points(15*(i - 1) + j), j=1, remainder)
    endif
    read (unit=stdin, fmt=*)

    allocate (self%real_space_hamiltonian_elements(self%num_R_points, self%bands, self%bands), &
              self%real_space_position_elements(self%num_R_points, self%bands, self%bands, 3), &
              stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "WannInt: Error #5: failure allocating real_space_hamiltonian_elements &
              &and real_space_position_elements. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    !Read Hamiltonian.
    allocate (dummyR(2), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "WannInt: Error #5: failure allocating dummyR. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%real_space_hamiltonian_elements = cmplx_0
    do irpts = 1, self%num_R_points
      read (unit=stdin, fmt=*) (self%R_points(irpts, i), i=1, 3)
      do i = 1, self%bands
        do j = 1, self%bands
          read (unit=stdin, fmt=*) dummy1, dummy2, dummyR(1), dummyR(2)
          !As pointed out in the W90 source code v3.1.0, get_oper.F90, ln 147, addition is required instead of equality.
          self%real_space_hamiltonian_elements(irpts, i, j) = &
            self%real_space_hamiltonian_elements(irpts, i, j) + cmplx(dummyR(1), dummyR(2), wp)
        enddo
      enddo
      read (unit=stdin, fmt=*)
    enddo
    deallocate (dummyR)

    !Read position.
    allocate (dummyR(6), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "WannInt: Error #5: failure allocating dummyR. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%real_space_position_elements = cmplx_0
    do irpts = 1, self%num_R_points
      read (unit=stdin, fmt=*) dummy1, dummy2, dummy3
      do i = 1, self%bands
        do j = 1, self%bands
          read (unit=stdin, fmt=*) dummy1, dummy2, dummyR(1), dummyR(2), dummyR(3), dummyR(4), dummyR(5), dummyR(6)
          !As pointed out in the W90 source code v3.1.0, get_oper.F90, ln 147, addition is required instead of equality.
          self%real_space_position_elements(irpts, i, j, 1) = &
            self%real_space_position_elements(irpts, i, j, 1) + cmplx(dummyR(1), dummyR(2), wp)
          self%real_space_position_elements(irpts, i, j, 2) = &
            self%real_space_position_elements(irpts, i, j, 2) + cmplx(dummyR(3), dummyR(4), wp)
          self%real_space_position_elements(irpts, i, j, 3) = &
            self%real_space_position_elements(irpts, i, j, 3) + cmplx(dummyR(5), dummyR(6), wp)
        enddo
      enddo
      if (irpts == self%num_R_points) cycle
      read (unit=stdin, fmt=*)
    enddo
    deallocate (dummyR)

    close (unit=stdin)

    if (present(fermi_energy)) self%e_fermi = fermi_energy

    self%is_initialized = .true.

  end subroutine construct_from_file

  pure elemental function name(self)
    class(crystal), intent(in) :: self
    character(len=120) :: name
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    name = trim(self%nm)
  end function name

  pure function direct_lattice_basis(self)
    class(crystal), intent(in) :: self
    real(wp) :: direct_lattice_basis(3, 3)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    direct_lattice_basis = self%direct_l_b
  end function direct_lattice_basis

  pure function reciprocal_lattice_basis(self)
    class(crystal), intent(in) :: self
    real(wp) :: reciprocal_lattice_basis(3, 3)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    reciprocal_lattice_basis = self%reciprocal_l_b
  end function reciprocal_lattice_basis

  pure function metric_tensor(self)
    class(crystal), intent(in) :: self
    real(wp) :: metric_tensor(3, 3)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    metric_tensor = self%m_tensor
  end function metric_tensor

  pure function cell_volume(self)
    class(crystal), intent(in) :: self
    real(wp) :: cell_volume
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    cell_volume = self%c_volume
  end function cell_volume

  pure integer function num_bands(self)
    class(crystal), intent(in) :: self
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    num_bands = self%bands
  end function num_bands

  pure integer function nrpts(self)
    class(crystal), intent(in) :: self
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    nrpts = self%num_R_points
  end function nrpts

  pure function rpts(self)
    class(crystal), intent(in) :: self
    integer :: rpts(self%num_R_points, 3)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    rpts = self%R_points
  end function rpts

  pure function deg_rpts(self)
    class(crystal), intent(in) :: self
    integer :: deg_rpts(self%num_R_points)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    deg_rpts = self%deg_R_points
  end function deg_rpts

  pure function get_real_space_hamiltonian_elements(self)
    class(crystal), intent(in) :: self
    complex(wp) :: get_real_space_hamiltonian_elements(self%num_R_points, self%bands, self%bands)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    get_real_space_hamiltonian_elements = self%real_space_hamiltonian_elements
  end function get_real_space_hamiltonian_elements

  pure function get_real_space_position_elements(self)
    class(crystal), intent(in) :: self
    complex(wp) :: get_real_space_position_elements(self%num_R_points, self%bands, self%bands, 3)
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    get_real_space_position_elements = self%real_space_position_elements
  end function get_real_space_position_elements

  pure function fermi_energy(self)
    class(crystal), intent(in) :: self
    real(wp) :: fermi_energy
    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."
    fermi_energy = self%e_fermi
  end function fermi_energy

  pure elemental logical function initialized(self)
    class(crystal), intent(in) :: self
    initialized = self%is_initialized
  end function initialized

  function h_branch_1(self, kpt, derivative, all)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: kpt(3)
    integer, intent(in) :: derivative
    logical, intent(in) :: all

    type(container), allocatable :: h_branch_1(:)

    integer :: i

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."

    if (derivative < 0) error stop &
    "WannInt: Error #7: calculating Hamiltonian: invalid value for requested derivative &
    &while requesting all derivatives."

    if (all) then
      allocate (h_branch_1(derivative + 1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "WannInt: Error #5: failure allocating h_branch_1. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      do i = 1, derivative + 1
        h_branch_1(i) = hamiltonian_der_driver(self, kpt, i - 1)
      enddo
    else
      allocate (h_branch_1(1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "WannInt: Error #5: failure allocating h_branch_1. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      h_branch_1(1) = hamiltonian_der_driver(self, kpt, derivative)
    endif

  end function h_branch_1

  function h_branch_2(self, kpt, derivative)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: kpt(3)
    integer, intent(in) :: derivative

    type(container) :: h_branch_2

    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."

    if (derivative < 0) error stop &
      "WannInt: Error #7: calculating Hamiltonian: invalid value for requested derivative."

    h_branch_2 = hamiltonian_der_driver(self, kpt, derivative)

  end function h_branch_2

  function h_branch_3(self, kpt)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: kpt(3)

    complex(wp) :: h_branch_3(self%bands, self%bands)
    type(container) :: h

    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."

    h = hamiltonian_der_driver(self, kpt, 0)

    h_branch_3 = reshape(h%cdp_storage, [self%bands, self%bands])

  end function h_branch_3

  function hamiltonian_der_driver(self, k, der)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: k(3)
    integer, intent(in) :: der !In calling procedure: ensure that der \in [0, inf].
    type(container) :: hamiltonian_der_driver

    integer :: tensor_shape(2 + der), &
               array_layout_handle(2 + der), &
               memory_layout_handle
    integer :: i, irpts

    complex(wp) :: sum

    real(wp) :: kdotr, vec(3), prod_r

    !Define the shape of the output array.
    !First two are band indices.
    !Each requested derivative adds a further spatial index in
    !the range [1, 3].
    tensor_shape(1) = self%bands
    tensor_shape(2) = self%bands
    if (der > 0) then
      do i = 3, size(tensor_shape)
        tensor_shape(i) = 3
      enddo
    endif

    call hamiltonian_der_driver%construct(container_type="complex_dp", &
                                          dimension_specifier=tensor_shape)

    !Loop though band indices and requested derivative indices.
    !_OMPTGT_(PARALLEL DO &)
    !_OMPTGT_(PRIVATE (memory_layout_handle, array_layout_handle, sum, kdotr, vec, prod_r))
    do memory_layout_handle = 1, hamiltonian_der_driver%size()

      array_layout_handle = hamiltonian_der_driver%ind(memory_layout_handle)
      sum = cmplx_0

      !_OMPTGT_(SIMD &)
      !_OMPTGT_(PRIVATE (kdotr, vec, prod_r) &)
      !_OMPTGT_(REDUCTION (+: sum))
      do irpts = 1, self%num_R_points

        !Compute factor appearing in the exponential.
        !(k is in coords relative to recip. lattice vectors).
        kdotr = 2.0_wp*pi*dot_product(self%R_points(irpts, :), k)

        !Compute Bravais lattice vector for label irpts.
        vec = 0.0_wp
        do i = 1, 3
          vec = vec + self%R_points(irpts, i)*self%direct_l_b(i, :)
        enddo

        prod_r = 1.0_wp
        if (der > 0) then
          do i = 1, der
            prod_r = prod_r*vec(array_layout_handle(i + 2))
          enddo
        endif

        sum = sum + & !Compute sum.
              ((cmplx_i)**real(der, wp))*(prod_r)*exp(cmplx_i*kdotr)* &
              self%real_space_hamiltonian_elements(irpts, array_layout_handle(1), array_layout_handle(2)) &
              /real(self%deg_R_points(irpts), wp)

      enddo
      !Store to result.
      call hamiltonian_der_driver%set(val=sum, at=memory_layout_handle)
    enddo

  end function hamiltonian_der_driver

  function a_branch_1(self, kpt, derivative, all)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: kpt(3)
    integer, intent(in) :: derivative
    logical, intent(in) :: all

    type(container), allocatable :: a_branch_1(:)

    integer :: i

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."

    if (derivative < 0) error stop &
    "WannInt: Error #7: calculating Berry connection: invalid value for requested derivative &
    &while requesting all derivatives."

    if (all) then
      allocate (a_branch_1(derivative + 1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "WannInt: Error #5: failure allocating a_branch_1. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      do i = 1, derivative + 1
        a_branch_1(i) = berry_der_driver(self, kpt, i - 1)
      enddo
    else
      allocate (a_branch_1(1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "WannInt: Error #5: failure allocating a_branch_1. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      a_branch_1(1) = berry_der_driver(self, kpt, derivative)
    endif

  end function a_branch_1

  function a_branch_2(self, kpt, derivative)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: kpt(3)
    integer, intent(in) :: derivative

    type(container) :: a_branch_2

    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."

    if (derivative < 0) error stop &
      "WannInt: Error #7: calculating Berry connection: invalid value for requested derivative."

    a_branch_2 = berry_der_driver(self, kpt, derivative)

  end function a_branch_2

  function a_branch_3(self, kpt)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: kpt(3)

    complex(wp) :: a_branch_3(self%bands, self%bands, 3)
    type(container) :: a

    if (.not. (self%is_initialized)) error stop &
      "WannInt: Error #6: crystal is not initialized."

    a = berry_der_driver(self, kpt, 0)

    a_branch_3 = reshape(a%cdp_storage, [self%bands, self%bands, 3])

  end function a_branch_3

  function berry_der_driver(self, k, der)
    class(crystal), intent(in) :: self
    real(wp), intent(in) :: k(3)
    integer, intent(in) :: der
    type(container) :: berry_der_driver

    integer :: tensor_shape(3 + der), &
               array_layout_handle(3 + der), &
               memory_layout_handle
    integer :: i, irpts

    complex(wp) :: sum

    real(wp) :: kdotr, vec(3), prod_r

    !Define the shape of the output array.
    !First two are band indices and third is spatial index.
    !Each requested derivative adds a further spatial index in
    !the range [1, 3].
    tensor_shape(1) = self%bands
    tensor_shape(2) = self%bands
    tensor_shape(3) = 3
    if (der > 0) then
      do i = 4, size(tensor_shape)
        tensor_shape(i) = 3
      enddo
    endif

    call berry_der_driver%construct(container_type="complex_dp", &
                                    dimension_specifier=tensor_shape)

    !Loop though band indices, spatial index and requested derivative indices.
    !_OMPTGT_(PARALLEL DO &)
    !_OMPTGT_(PRIVATE (memory_layout_handle, array_layout_handle, sum, kdotr, vec, prod_r))
    do memory_layout_handle = 1, berry_der_driver%size()

      array_layout_handle = berry_der_driver%ind(memory_layout_handle)
      sum = cmplx_0

      !_OMPTGT_(SIMD &)
      !_OMPTGT_(PRIVATE (kdotr, vec, prod_r) &)
      !_OMPTGT_(REDUCTION (+: sum))
      do irpts = 1, self%num_R_points

        !Compute factor appearing in the exponential.
        !(k is in coords relative to recip. lattice vectors).
        kdotr = 2.0_wp*pi*dot_product(self%R_points(irpts, :), k)

        !Compute Bravais lattice vector for label irpts.
        vec = 0.0_wp
        do i = 1, 3
          vec = vec + self%R_points(irpts, i)*self%direct_l_b(i, :)
        enddo

        prod_r = 1.0_wp
        if (der > 0) then
          do i = 1, der
            prod_r = prod_r*vec(array_layout_handle(i + 3))
          enddo
        endif

        sum = sum + & !Compute sum.
              ((cmplx_i)**real(der, wp))*(prod_r)*exp(cmplx_i*kdotr)* &
              self%real_space_position_elements(irpts, array_layout_handle(1), array_layout_handle(2), array_layout_handle(3)) &
              /real(self%deg_R_points(irpts), wp)

      enddo
      !Store to result.
      call berry_der_driver%set(val=sum, at=memory_layout_handle)
    enddo

  end function berry_der_driver

end module WannInt
