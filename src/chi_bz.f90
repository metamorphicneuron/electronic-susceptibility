!===============================================================================
! chi_bz.f90  ->  bin/chi_bz.x
!
! Puts chi0(q) from lindhard_chi.x (chi_<E_F>_<T>.dat, or a k-resolved
! chi_k<i>_<E_F>_<T>.dat) into the true geometry of reciprocal space, together
! with the Brillouin-zone boundaries, so that chi0 can be plotted as an overlay
! on the Brillouin zone and compared with experiments: q in 1/Angstrom (STM-FFT,
! ARPES, electron diffraction) and in reciprocal-lattice units H K L (X-ray
! diffuse scattering, inelastic scattering).
! Serial; run it in the folder of the chi file, as often as needed:
!
!   bin/chi_bz.x bz.in            (or  bin/chi_bz.x < bz.in)
!
! Units: 1/Angstrom with the 2 pi included (|b| = 2 pi / a). qx, qy are the
! cartesian axes when the plane is a cartesian plane (normal along x, y or z),
! e.g. qx, qy for a 2D material in the xy plane; otherwise qx runs along the
! first in-plane reciprocal axis and qy perpendicular to it. The headers give
! both axes, the reciprocal axes, the largest values with their q, and ready
! gnuplot lines; examples/plot_bz.py draws the overlay as a PDF, asking in the
! terminal for colour scale, colour map, style and font size.
! Errors are written to standard output and standard error.
!===============================================================================

program chi_bz
  use, intrinsic :: iso_fortran_env, only: dp => real64, input_unit, output_unit, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  implicit none

  real(dp), parameter :: PI = 3.14159265358979323846_dp
  real(dp), parameter :: BOHR = 0.529177210903_dp           ! Angstrom

  type :: string_t
    character(len=:), allocatable :: s
  end type string_t

  type :: poly_t
    real(dp), allocatable :: p(:,:)
  end type poly_t

  ! input
  character(len=:), allocatable :: input_name, chi_file, nscf_out, prefix, plane, zone, zone_shape
  logical :: box_zone                                      ! zone = the cell -1/2 ... 1/2 in crystal coordinates
  real(dp) :: extent, fixc(3), avec(3,3)                   ! avec: columns a1 a2 a3 (Angstrom)
  logical :: have_a(3), fix_given(3)
  type(string_t), allocatable :: plab(:)                   ! labelled points ...
  real(dp), allocatable :: pk(:,:)                         ! ... (3, np), crystal coordinates

  ! the chi file
  integer :: nmesh(3), nt
  real(dp) :: per(3)
  real(dp), allocatable :: val(:,:)                        ! (2, nt) Re, Im on the mesh sites
  type(string_t), allocatable :: chi_head(:)               ! its '#' lines

  ! geometry
  real(dp) :: bvec(3,3)                                    ! columns b1 b2 b3, 1/Angstrom (2 pi included)
  real(dp), allocatable :: hv(:,:)                         ! reciprocal-lattice vectors /= 0 near Gamma
  real(dp) :: nan

  nan = ieee_value(nan, ieee_quiet_nan)
  write(output_unit, '(a)') 'chi_bz: chi0 of lindhard_chi.x on the Brillouin zone'
  call read_input()
  call read_chi()
  call setup_cell()
  if (count(nmesh > 1) == 0) call fatal(chi_file // ': the mesh is a single point')
  if (count(nmesh > 1) == 1) then
    call line_output()
  else
    call plane_output()
  end if

contains

  !-----------------------------------------------------------------------------
  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    write(output_unit, '(a)') 'ERROR: ' // msg
    flush(output_unit)
    write(error_unit, '(a)') 'ERROR: ' // msg
    flush(error_unit)
    stop 1
  end subroutine fatal

  !=============================================================================
  ! input
  !=============================================================================
  subroutine read_input()
    character(len=:), allocatable :: line, text, key
    type(string_t), allocatable :: tok(:)
    logical, allocatable :: quoted(:)
    integer :: u, ios, n, i, j, ntok
    logical :: started

    if (command_argument_count() >= 1) then
      call get_command_argument(1, length=n)
      allocate(character(len=n) :: input_name)
      call get_command_argument(1, input_name)
      open(newunit=u, file=input_name, status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal('cannot open the input file "' // input_name // '"')
    else
      input_name = '(standard input)'
      u = input_unit
    end if
    text = ''
    do
      call read_line(u, line, ios)
      if (ios /= 0) exit
      text = text // ' ' // strip_comment(line)
    end do
    if (u /= input_unit) close(u)
    call tokenize(text, tok, quoted)

    chi_file = ''
    nscf_out = ''
    prefix = ''
    plane = ''
    zone = 'extended'
    zone_shape = 'auto'
    extent = 1.0_dp
    fixc = 0.0_dp
    fix_given = .false.
    have_a = .false.
    avec = 0.0_dp
    allocate(plab(0), pk(3, 0))

    ntok = size(tok)
    started = .false.
    i = 1
    do while (i <= ntok)
      if (.not. quoted(i) .and. lower(tok(i)%s) == '&bz') then
        started = .true.
        exit
      end if
      i = i + 1
    end do
    if (.not. started) call fatal(input_name // ': no section &bz')
    i = i + 1
    do
      if (i > ntok) call fatal(input_name // ': the section &bz is not closed with /')
      if (.not. quoted(i) .and. tok(i)%s == '/') exit
      if (i + 1 > ntok .or. quoted(i)) call fatal(input_name // ': name = value expected at "' // tok(i)%s // '"')
      if (quoted(i + 1) .or. tok(i + 1)%s /= '=') &
        call fatal(input_name // ': name = value expected at "' // tok(i)%s // '"')
      key = tok(i)%s
      j = i + 2                                    ! the values run to the next  name =  or to /
      do while (j <= ntok)
        if (.not. quoted(j) .and. tok(j)%s == '/') exit
        if (j + 1 <= ntok) then
          if (.not. quoted(j + 1) .and. tok(j + 1)%s == '=') exit
        end if
        j = j + 1
      end do
      call set_key(key, tok(i + 2:j - 1))
      i = j
    end do

    if (len(chi_file) == 0) call fatal(input_name // ': chi_file is required (a chi_<E_F>_<T>.dat of lindhard_chi.x)')
    if (zone /= 'extended' .and. zone /= 'first') call fatal('zone = ''extended'' or ''first'', got ''' // zone // '''')
    if (zone_shape /= 'auto' .and. zone_shape /= 'box' .and. zone_shape /= 'wigner-seitz') &
      call fatal('zone_shape = ''auto'', ''box'' or ''wigner-seitz'', got ''' // zone_shape // '''')
    if (len(plane) > 0 .and. plane /= 'ab' .and. plane /= 'bc' .and. plane /= 'ac') &
      call fatal('plane = ''ab'', ''bc'' or ''ac'', got ''' // plane // '''')
    if (extent <= 0.0_dp) call fatal('extent must be > 0')
  end subroutine read_input

  !> one entry  name = values
  subroutine set_key(key, v)
    character(len=*), intent(in) :: key
    type(string_t), intent(in) :: v(:)
    character(len=:), allocatable :: k
    real(dp) :: x(3)
    integer :: i

    if (size(v) == 0) call fatal(input_name // ': no value for ' // key)
    k = lower(key)
    select case (k)
    case ('chi_file')
      chi_file = one_text(key, v)
      return
    case ('nscf_out')
      nscf_out = one_text(key, v)
      return
    case ('prefix')
      prefix = one_text(key, v)
      return
    case ('plane')
      plane = lower(one_text(key, v))
      return
    case ('zone_shape')
      zone_shape = lower(one_text(key, v))
      if (zone_shape == 'ws') zone_shape = 'wigner-seitz'
      return
    case ('zone')
      zone = lower(one_text(key, v))
      return
    case ('extent')
      if (size(v) /= 1) call fatal(input_name // ': extent takes one number')
      extent = number(key, v(1)%s)
      return
    case ('a1', 'a2', 'a3')
      if (size(v) /= 3) call fatal(input_name // ': ' // key // ' takes three numbers (Angstrom)')
      do i = 1, 3
        x(i) = number(key, v(i)%s)
      end do
      i = ichar(k(2:2)) - ichar('0')
      avec(:, i) = x
      have_a(i) = .true.
      return
    case ('k1', 'k2', 'k3')
      if (size(v) == 1) then                       ! fixed coordinate of the plane; with 3 numbers a point
        i = ichar(k(2:2)) - ichar('0')
        fixc(i) = number(key, v(1)%s)
        fix_given(i) = .true.
        return
      end if
    end select
    if (size(v) /= 3) call fatal(input_name // ': "' // key // '" is not a known name, and a point needs ' // &
      'three numbers (label = k1 k2 k3); got ' // istr(size(v)) // ' value(s)')
    do i = 1, 3
      x(i) = number(key, v(i)%s)
    end do
    call push_point(key, x)
  end subroutine set_key

  function one_text(key, v) result(s)
    character(len=*), intent(in) :: key
    type(string_t), intent(in) :: v(:)
    character(len=:), allocatable :: s
    if (size(v) /= 1) call fatal(input_name // ': ' // key // ' takes one value')
    s = v(1)%s
  end function one_text

  !> a number, or a fraction such as 1/3
  function number(key, s) result(x)
    character(len=*), intent(in) :: key, s
    real(dp) :: x, a, b
    integer :: ios, k
    k = index(s, '/')
    if (k > 1 .and. k < len(s)) then
      read(s(1:k - 1), *, iostat=ios) a
      if (ios == 0) then
        read(s(k + 1:), *, iostat=ios) b
        if (ios == 0) then
          if (abs(b) > 0.0_dp) then
            x = a / b
            return
          end if
        end if
      end if
    else if (k == 0) then
      read(s, *, iostat=ios) x
      if (ios == 0) return
    end if
    call fatal(input_name // ': ' // key // ': "' // s // '" is not a number')
    x = 0.0_dp
  end function number

  subroutine push_point(label, x)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: x(3)
    type(string_t), allocatable :: t(:)
    real(dp), allocatable :: p(:,:)
    integer :: n, i
    n = size(plab)
    allocate(t(n + 1), p(3, n + 1))
    do i = 1, n
      t(i)%s = plab(i)%s
    end do
    t(n + 1)%s = label
    p(:, 1:n) = pk
    p(:, n + 1) = x
    call move_alloc(t, plab)
    call move_alloc(p, pk)
  end subroutine push_point

  !> words, quoted texts, '=' and '/'; blanks and commas separate
  subroutine tokenize(text, tok, quoted)
    character(len=*), intent(in) :: text
    type(string_t), allocatable, intent(out) :: tok(:)
    logical, allocatable, intent(out) :: quoted(:)
    character(len=*), parameter :: STOPS = ' ,=''"' // char(9)
    integer :: i, j, n
    character :: c

    allocate(tok(0), quoted(0))
    n = len(text)
    i = 1
    do while (i <= n)
      c = text(i:i)
      if (c == ' ' .or. c == ',' .or. c == char(9)) then
        i = i + 1
      else if (c == '=') then
        call push_token(tok, quoted, '=', .false.)
        i = i + 1
      else if (c == '''' .or. c == '"') then
        j = index(text(i + 1:), c)
        if (j == 0) call fatal(input_name // ': a quote is not closed')
        call push_token(tok, quoted, text(i + 1:i + j - 1), .true.)
        i = i + j + 1
      else
        j = i
        do while (j <= n)
          if (index(STOPS, text(j:j)) > 0) exit
          j = j + 1
        end do
        if (j - i > 1 .and. text(j - 1:j - 1) == '/') then   ! e.g.  0.5/  : the value, then the end
          call push_token(tok, quoted, text(i:j - 2), .false.)
          call push_token(tok, quoted, '/', .false.)
        else
          call push_token(tok, quoted, text(i:j - 1), .false.)
        end if
        i = j
      end if
    end do
  end subroutine tokenize

  subroutine push_token(tok, quoted, s, q)
    type(string_t), allocatable, intent(inout) :: tok(:)
    logical, allocatable, intent(inout) :: quoted(:)
    character(len=*), intent(in) :: s
    logical, intent(in) :: q
    type(string_t), allocatable :: t(:)
    integer :: n, i
    n = size(tok)
    allocate(t(n + 1))
    do i = 1, n
      t(i)%s = tok(i)%s
    end do
    t(n + 1)%s = s
    call move_alloc(t, tok)
    quoted = [quoted, q]
  end subroutine push_token

  !=============================================================================
  ! the chi file of lindhard_chi.x
  !=============================================================================
  subroutine read_chi()
    character(len=:), allocatable :: line
    integer :: u, ios, ncol, k, g(3), s, nrow, nmiss
    real(dp) :: row(8), kc(3)
    logical :: have_mesh

    open(newunit=u, file=chi_file, status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('cannot open chi_file "' // chi_file // '"')
    allocate(chi_head(0))
    have_mesh = .false.
    ncol = 0
    nrow = 0
    do
      call read_line(u, line, ios)
      if (ios /= 0) exit
      line = trim(adjustl(line))
      if (len(line) == 0) cycle
      if (line(1:1) == '#') then
        call push_head(line)
        if (index(line, '# mesh ') == 1) then
          read(line(8:), *, iostat=ios) nmesh
          k = index(line, 'period')
          if (ios == 0 .and. k > 0) read(line(k + 6:), *, iostat=ios) per
          if (ios /= 0 .or. any(nmesh < 1)) call fatal(chi_file // ': cannot read the line "' // line // '"')
          if (any(abs(per) < 1.0e-12_dp)) call fatal(chi_file // ': period 0 in "' // line // '"')
          have_mesh = .true.
          nt = product(nmesh)
          allocate(val(2, nt))
          val = nan
        end if
        cycle
      end if
      if (.not. have_mesh) call fatal(chi_file // ': no line "# mesh N1 N2 N3, period ..." before the data; ' // &
                                      'is this a file of lindhard_chi.x?')
      if (ncol == 0) then
        ncol = count_tokens(line)
        if (ncol /= 5 .and. ncol /= 8) call fatal(chi_file // ': 5 columns (q1 q2 q3 Re Im) or 8 (k-resolved file) ' // &
          'expected, found ' // istr(ncol) // '. Use a file of this version of lindhard_chi.x')
      end if
      read(line, *, iostat=ios) row(1:ncol)
      if (ios /= 0) call fatal(chi_file // ': cannot read the line "' // line // '"')
      kc = row(1:3) / per * nmesh
      if (any(abs(kc - nint(kc)) > 1.0e-3_dp)) call fatal(chi_file // ': the q-point of the line "' // line // &
                                                           '" is not on the mesh')
      g = modulo(nint(kc), nmesh)
      s = site(g)
      val(:, s) = row(ncol - 1:ncol)
      nrow = nrow + 1
    end do
    close(u)
    if (nrow == 0) call fatal(chi_file // ': no data')
    nmiss = count(ieee_is_nan(val(1, :)))
    write(output_unit, '(a)') '  chi file    ' // chi_file // ': ' // istr(nrow) // ' lines, mesh ' // istr(nmesh(1)) // &
      ' x ' // istr(nmesh(2)) // ' x ' // istr(nmesh(3)) // pick(ncol == 8, ' (k-resolved file)', '')
    if (nmiss > 0) write(output_unit, '(a)') '  NOTE: ' // istr(nmiss) // ' mesh points have no value (NaN in the output)'
  end subroutine read_chi

  subroutine push_head(line)
    character(len=*), intent(in) :: line
    type(string_t), allocatable :: t(:)
    integer :: n, i
    n = size(chi_head)
    allocate(t(n + 1))
    do i = 1, n
      t(i)%s = chi_head(i)%s
    end do
    t(n + 1)%s = line
    call move_alloc(t, chi_head)
  end subroutine push_head

  pure integer function site(g)
    integer, intent(in) :: g(3)
    site = 1 + g(1) + nmesh(1) * (g(2) + nmesh(2) * g(3))
  end function site

  !=============================================================================
  ! the cell
  !=============================================================================
  subroutine setup_cell()
    real(dp) :: at_t(3,3)
    integer :: n1, n2, n3, nop

    if (len(nscf_out) > 0) then
      if (any(have_a)) call fatal('give either nscf_out or a1, a2, a3, not both')
      call read_cell(nscf_out)
    else if (.not. all(have_a)) then
      call fatal('the cell is needed: nscf_out = ''nscf.out'' (any pw.x output of this crystal), or a1, a2, a3 ' // &
                 '(cartesian, Angstrom)')
    end if
    if (abs(det3(avec)) < 1.0e-8_dp) call fatal('a1, a2, a3 are linearly dependent')
    at_t = transpose(avec)
    bvec = 2.0_dp * PI * inv3(at_t)                       ! a_i . b_j = 2 pi delta_ij

    allocate(hv(3, 0))
    do n1 = -3, 3
      do n2 = -3, 3
        do n3 = -3, 3
          if (n1 == 0 .and. n2 == 0 .and. n3 == 0) cycle
          hv = reshape([reshape(hv, [size(hv)]), matmul(bvec, real([n1, n2, n3], dp))], [3, size(hv, 2) + 1])
        end do
      end do
    end do
    write(output_unit, '(a)') '  cell        ' // pick(len(nscf_out) > 0, 'from ' // nscf_out, 'a1, a2, a3 of the input')
    write(output_unit, '(a)') '              a1 = (' // vtext(avec(:, 1)) // ')  a2 = (' // vtext(avec(:, 2)) // &
      ')  a3 = (' // vtext(avec(:, 3)) // ') Angstrom'
    write(output_unit, '(a)') '              b1 = (' // vtext(bvec(:, 1)) // ')  b2 = (' // vtext(bvec(:, 2)) // &
      ')  b3 = (' // vtext(bvec(:, 3)) // ') 1/Angstrom'
    nop = lattice_order(avec)
    select case (zone_shape)
    case ('box')
      box_zone = .true.
    case ('wigner-seitz')
      box_zone = .false.
    case default                                          ! monoclinic and triclinic: the box
      box_zone = nop <= 4
    end select
    write(output_unit, '(a)') '  lattice     ' // istr(nop) // ' point operations (' // lattice_name(nop) // '); zone: ' // &
      pick(box_zone, 'the cell -1/2 ... 1/2 in crystal coordinates (high-symmetry points with coordinates 0 and 1/2 ' // &
      'on its border)', 'Wigner-Seitz cell') // pick(zone_shape == 'auto', '', ' (zone_shape)')
  end subroutine setup_cell

  !> the point operations of the lattice: integer matrices (entries -1, 0, 1) that keep the metric.
  !> 2 triclinic, 4 monoclinic, 8 orthorhombic, 12 rhombohedral, 16 tetragonal, 24 hexagonal, 48 cubic
  integer function lattice_order(a)
    real(dp), intent(in) :: a(3,3)
    real(dp) :: gm(3,3), m(3,3), tol
    integer :: idx, k, j, e(9)
    gm = matmul(transpose(a), a)
    tol = 1.0e-4_dp * maxval(abs(gm))
    lattice_order = 0
    do idx = 0, 3**9 - 1
      k = idx
      do j = 1, 9
        e(j) = mod(k, 3) - 1
        k = k / 3
      end do
      m = reshape(real(e, dp), [3, 3])
      if (maxval(abs(matmul(transpose(m), matmul(gm, m)) - gm)) < tol) lattice_order = lattice_order + 1
    end do
  end function lattice_order

  function lattice_name(nop) result(s)
    integer, intent(in) :: nop
    character(len=:), allocatable :: s
    select case (nop)
    case (2)
      s = 'triclinic'
    case (4)
      s = 'monoclinic'
    case (8)
      s = 'orthorhombic'
    case (12)
      s = 'rhombohedral'
    case (16)
      s = 'tetragonal'
    case (24)
      s = 'hexagonal'
    case (48)
      s = 'cubic'
    case default
      s = 'cell not reduced?'
    end select
  end function lattice_name

  !> alat and the crystal axes of a pw.x output
  subroutine read_cell(fname)
    character(len=*), intent(in) :: fname
    character(len=:), allocatable :: line
    integer :: u, ios, i, k1, k2
    real(dp) :: alat
    logical :: got_alat, got_at

    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('cannot open nscf_out "' // fname // '"')
    got_alat = .false.
    got_at = .false.
    alat = 0.0_dp
    do
      call read_line(u, line, ios)
      if (ios /= 0) exit
      if (.not. got_alat .and. index(line, 'lattice parameter (alat)') > 0) then
        k1 = index(line, '=')
        read(line(k1 + 1:), *, iostat=ios) alat
        got_alat = (ios == 0 .and. alat > 0.0_dp)
      else if (.not. got_at .and. index(line, 'crystal axes: (cart. coord. in units of alat)') > 0) then
        do i = 1, 3
          call read_line(u, line, ios)
          if (ios /= 0) exit
          k1 = index(line, '(', back=.true.)
          k2 = index(line, ')', back=.true.)
          if (k1 == 0 .or. k2 <= k1) then
            ios = 1
            exit
          end if
          read(line(k1 + 1:k2 - 1), *, iostat=ios) avec(:, i)
          if (ios /= 0) exit
        end do
        got_at = (ios == 0)
      end if
      if (got_alat .and. got_at) exit
    end do
    close(u)
    if (.not. (got_alat .and. got_at)) call fatal(fname // ': ''lattice parameter (alat)'' and the ''crystal axes'' ' // &
      'of a pw.x output not found')
    avec = avec * alat * BOHR
  end subroutine read_cell

  !=============================================================================
  ! plane (2D and 3D systems)
  !=============================================================================
  subroutine plane_output()
    character(len=3), parameter :: LETTERS = 'HKL'
    integer :: ia, ib, ic, gc, k, n, gcs(2)
    real(dp) :: c
    character(len=:), allocatable :: base

    ! ---- the plane: axes ia, ib in the plane, ic fixed
    if (count(nmesh > 1) == 2) then
      ic = minloc(nmesh, dim=1)                             ! the axis with N = 1
      if (len(plane) > 0 .and. plane /= pick(ic == 3, 'ab', pick(ic == 1, 'bc', 'ac'))) &
        write(output_unit, '(a)') '  NOTE: plane is ignored: the mesh is two-dimensional'
      plane = pick(ic == 3, 'ab', pick(ic == 1, 'bc', 'ac'))
    else if (len(plane) == 0) then
      plane = 'ab'
      write(output_unit, '(a)') '  NOTE: 3D mesh and no plane given: plane = ''ab'''
    end if
    select case (plane)
    case ('ab')
      ia = 1; ib = 2; ic = 3
    case ('bc')
      ia = 2; ib = 3; ic = 1
    case default
      ia = 1; ib = 3; ic = 2
    end select
    do k = 1, 3
      if (fix_given(k) .and. k /= ic) write(output_unit, '(a)') '  NOTE: k' // istr(k) // ' is ignored for plane ''' // &
        plane // ''''
    end do
    gc = nint(fixc(ic) * nmesh(ic))
    c = real(gc, dp) / nmesh(ic)
    if (abs(c - fixc(ic)) > 1.0e-9_dp) write(output_unit, '(a)') '  NOTE: k' // istr(ic) // ' = ' // &
      fstr(fixc(ic), '(f10.6)') // ' is not a mesh plane; the nearest one, k' // istr(ic) // ' = ' // fstr(c, '(f10.6)') // &
      ', is used'

    ! the planes k_ic = +c and k_ic = -c; one plane when they are the same (c = 0, or 1/2 = -1/2)
    gcs = [gc, -gc]
    n = merge(1, 2, modulo(2 * gc, nmesh(ic)) == 0)
    base = pick(len(prefix) > 0, prefix, stem(chi_file) // '_' // plane)
    if (n == 2) then
      write(output_unit, '(a)') '  planes      k' // istr(ic) // ' = +' // fstr(abs(c), '(f10.6)') // &
        ' and -' // fstr(abs(c), '(f10.6)') // ': written together, in one file'
      call one_plane(ia, ib, ic, gcs, base // '_' // LETTERS(ic:ic) // '+-' // plain_text(abs(c)))
    else
      call one_plane(ia, ib, ic, gcs(1:1), base)
    end if
  end subroutine plane_output

  !-----------------------------------------------------------------------------
  !> the map, the zone boundaries and the points of the plane k_ic = c (gcs = [gc]), or of the two
  !> planes +c and -c together (gcs = [gc, -gc]): one map with the values of both. The in-plane
  !> coordinates are those of H_a b_a + H_b b_b, so that the point H_a = H_b = 0 of each plane is at
  !> (0, 0), both planes share the grid, and a point (H K L) is drawn at its own H_a, H_b
  subroutine one_plane(ia, ib, ic, gcs, pfx)
    integer, intent(in) :: ia, ib, ic, gcs(:)
    character(len=*), intent(in) :: pfx
    character(len=3), parameter :: LETTERS = 'HKL'
    integer :: nps, nz, i, j, i0, i1, j0, j1, u, k, ip, npoly, na, nb, nc, kk, p
    real(dp) :: cs(2), nrm(3), e1(3), e2(3), on(2), sh(2,2), q(3), x, y, hkl(3), v(2,2), g0(3), st(2), mab(2,2)
    real(dp) :: xlo, xhi, ylo, yhi, scale
    real(dp) :: best(2,2), bestq(2,2), bestx(2,2,2), besth(3,2,2), far(2,2), farq(2,2), farx(2,2,2), farh(3,2,2)
    real(dp), allocatable :: poly(:,:)
    type(poly_t) :: zone1(2)
    logical :: near(2)
    character(len=2) :: names(2)
    character(len=40) :: ptag(2)
    character(len=:), allocatable :: fmap, fbz, fpt, offnote

    nps = size(gcs)
    do ip = 1, nps
      cs(ip) = real(gcs(ip), dp) / nmesh(ic)
      ptag(ip) = 'k' // istr(ic) // ' = ' // fstr(cs(ip), '(f10.6)')
    end do

    ! ---- in-plane cartesian axes e1, e2 and the normal
    nrm = cross(bvec(:, ia), bvec(:, ib))
    nrm = nrm / norm2(nrm)
    if (abs(abs(nrm(3)) - 1.0_dp) < 1.0e-6_dp) then
      e1 = [1.0_dp, 0.0_dp, 0.0_dp]
      e2 = [0.0_dp, 1.0_dp, 0.0_dp]
      names = ['qx', 'qy']
    else if (abs(abs(nrm(1)) - 1.0_dp) < 1.0e-6_dp) then
      e1 = [0.0_dp, 1.0_dp, 0.0_dp]
      e2 = [0.0_dp, 0.0_dp, 1.0_dp]
      names = ['qy', 'qz']
    else if (abs(abs(nrm(2)) - 1.0_dp) < 1.0e-6_dp) then
      e1 = [1.0_dp, 0.0_dp, 0.0_dp]
      e2 = [0.0_dp, 0.0_dp, 1.0_dp]
      names = ['qx', 'qz']
    else
      e1 = bvec(:, ia) / norm2(bvec(:, ia))
      e2 = cross(nrm, e1)
      names = ['q1', 'q2']
    end if
    nrm = cross(e1, e2)                                     ! right-handed with e1, e2
    mab(1, :) = [dot_product(bvec(:, ia), e1), dot_product(bvec(:, ib), e1)]   ! (x, y) = M (H_a, H_b)
    mab(2, :) = [dot_product(bvec(:, ia), e2), dot_product(bvec(:, ib), e2)]
    scale = maxval(norm2(bvec, dim=1))

    fmap = pfx // '_map.dat'
    fbz = pfx // '_bz.dat'
    fpt = pfx // '_points.dat'

    ! ---- the first zone in each plane
    do ip = 1, nps
      on(ip) = cs(ip) * dot_product(bvec(:, ic), nrm)                        ! distance from Gamma
      sh(:, ip) = cs(ip) * [dot_product(bvec(:, ic), e1), dot_product(bvec(:, ic), e2)]   ! (0 0 c), in-plane part
      if (box_zone) then
        zone1(ip)%p = box_poly(mab)
      else
        call zone_section([0.0_dp, 0.0_dp, 0.0_dp], e1, e2, nrm, on(ip), poly)
        if (size(poly, 2) < 3) call fatal('the plane ' // trim(ptag(ip)) // ' does not cut the first Brillouin zone')
        poly(1, :) = poly(1, :) - sh(1, ip)
        poly(2, :) = poly(2, :) - sh(2, ip)
        zone1(ip)%p = poly
      end if
    end do
    nz = nps                                                ! different zone outlines
    if (nps == 2) then
      if (same_poly(zone1(1)%p, zone1(2)%p, 1.0e-7_dp * scale)) nz = 1
    end if

    ! ---- the window of the map
    if (zone == 'first') then
      i0 = huge(1)
      i1 = -huge(1)
      j0 = huge(1)
      j1 = -huge(1)
      do ip = 1, nz
        do k = 1, size(zone1(ip)%p, 2)
          st = solve2(mab, zone1(ip)%p(:, k))
          i0 = min(i0, floor(st(1) * nmesh(ia) + 1.0e-6_dp))
          i1 = max(i1, ceiling(st(1) * nmesh(ia) - 1.0e-6_dp))
          j0 = min(j0, floor(st(2) * nmesh(ib) + 1.0e-6_dp))
          j1 = max(j1, ceiling(st(2) * nmesh(ib) - 1.0e-6_dp))
        end do
      end do
    else
      i1 = floor(extent * nmesh(ia) + 1.0e-6_dp)
      j1 = floor(extent * nmesh(ib) + 1.0e-6_dp)
      i0 = -i1
      j0 = -j1
    end if

    ! ---- largest values of each plane: anywhere, and away from Gamma and its equivalents
    best = -huge(1.0_dp)
    far = -huge(1.0_dp)
    bestq = huge(1.0_dp)
    farq = huge(1.0_dp)
    bestx = 0.0_dp
    farx = 0.0_dp
    besth = 0.0_dp
    farh = 0.0_dp
    xlo = huge(1.0_dp)
    xhi = -huge(1.0_dp)
    ylo = huge(1.0_dp)
    yhi = -huge(1.0_dp)
    do i = i0, i1
      do j = j0, j1
        call pair_point(i, j, ia, ib, ic, gcs, mab, x, y, hkl, v, near)
        xlo = min(xlo, x)
        xhi = max(xhi, x)
        ylo = min(ylo, y)
        yhi = max(yhi, y)
        do ip = 1, nps
          hkl(ic) = cs(ip)
          q = matmul(bvec, hkl)
          do p = 1, 2
            if (ieee_is_nan(v(p, ip))) cycle
            call keep_max(v(p, ip), norm2(q), x, y, hkl, best(p, ip), bestq(p, ip), bestx(:, p, ip), besth(:, p, ip))
            if (.not. near(ip)) call keep_max(v(p, ip), norm2(q), x, y, hkl, far(p, ip), farq(p, ip), farx(:, p, ip), &
                                              farh(:, p, ip))
          end do
        end do
      end do
    end do

    ! ---- zone boundaries; with two different outlines, one section per plane
    open(newunit=u, file=fbz, status='replace', action='write')
    write(u, '(a)') '# Brillouin-zone boundaries in the plane(s) ' // trim(ptag(1)) // pick(nps == 2, ' and ' // &
      trim(ptag(2)), '') // ' (spanned by b' // istr(ia) // ', b' // istr(ib) // '); one closed polygon per zone; ' // &
      pick(box_zone, 'zone = the cell -1/2 ... 1/2 in crystal coordinates', 'Wigner-Seitz zones')
    write(u, '(a)') '# columns: ' // names(1) // ' ' // names(2) // ' [1/Angstrom], (0, 0) = the point ' // &
      'H_a = H_b = 0 of the plane'
    npoly = 0
    do ip = 1, nz
      if (nps == 2) write(u, '(a)') '# plane ' // pick(nz == 1, 'both', signed_text(cs(ip)))
      if (zone == 'first') then
        call write_poly(u, zone1(ip)%p)
        npoly = npoly + 1
        cycle
      end if
      kk = ceiling(extent) + 2
      do nc = merge(0, -3, box_zone), merge(0, 3, box_zone)
        do nb = -kk, kk
          do na = -kk, kk
            if (box_zone) then
              poly = zone1(ip)%p
              poly(1, :) = poly(1, :) + mab(1, 1) * na + mab(1, 2) * nb
              poly(2, :) = poly(2, :) + mab(2, 1) * na + mab(2, 2) * nb
            else
              g0 = na * bvec(:, ia) + nb * bvec(:, ib) + nc * bvec(:, ic)
              call zone_section(g0, e1, e2, nrm, on(ip), poly)
              if (size(poly, 2) < 3) cycle
              poly(1, :) = poly(1, :) - sh(1, ip)
              poly(2, :) = poly(2, :) - sh(2, ip)
            end if
            if (maxval(poly(1, :)) < xlo .or. minval(poly(1, :)) > xhi .or. &
                maxval(poly(2, :)) < ylo .or. minval(poly(2, :)) > yhi) cycle
            call write_poly(u, poly)
            npoly = npoly + 1
          end do
        end do
      end do
    end do
    close(u)

    ! ---- the map
    open(newunit=u, file=fmap, status='replace', action='write')
    write(u, '(a)') '# kind: plane' // pick(nps == 2, ' pair', '')
    if (nps == 2) write(u, '(a)') '# pair: ' // LETTERS(ic:ic) // ' ' // signed_text(cs(1)) // ' ' // signed_text(cs(2))
    write(u, '(a)') '# chi0 of ' // chi_file // ' on the plane ' // trim(ptag(1)) // ' ' // hkl_name(ia, ib, ic, cs(1)) // &
      pick(nps == 2, ' and on the plane ' // trim(ptag(2)) // ' ' // hkl_name(ia, ib, ic, cs(2)) // ', on one grid', '') // &
      ', ' // pick(zone == 'first', 'first Brillouin zone (NaN outside it)', 'extended zones, k' // istr(ia) // ' and k' // &
      istr(ib) // ' from ' // fstr(-extent, '(f8.3)') // ' to ' // fstr(extent, '(f8.3)'))
    write(u, '(a)') '# zone: ' // pick(box_zone, 'the cell -1/2 ... 1/2 in crystal coordinates (monoclinic, triclinic)', &
      'Wigner-Seitz')
    write(u, '(a)') '# values of the mesh points (no interpolation), repeated with the periodicity of the reciprocal lattice'
    call cell_header(u)
    write(u, '(a)') '# in-plane coordinates of H_a b_a + H_b b_b (cartesian): ' // names(1) // ' along (' // vtext(e1) // &
      '), ' // names(2) // ' along (' // vtext(e2) // '); (0, 0) = the point H_a = H_b = 0 of the plane'
    write(u, '(a)') '# ' // istr(i1 - i0 + 1) // ' x ' // istr(j1 - j0 + 1) // ' points; k' // istr(ia) // &
      ' runs slowest, a blank line after each of its values'
    do ip = 1, nps
      if (nps == 2) write(u, '(a)') '# plane ' // trim(ptag(ip)) // ':'
      call max_header(u, '# ', best(:, ip), bestx(:, :, ip), besth(:, :, ip), far(:, ip), farx(:, :, ip), &
                      farh(:, :, ip), names)
    end do
    if (nps == 1) then
      write(u, '(a)') '# gnuplot: set view map; set size ratio -1; splot ''' // fmap // ''' u 1:2:3 w pm3d, ''' // &
        fbz // ''' u 1:2:(0) w l lc ''white'' notitle, ''' // fpt // ''' u 1:2:(0):3 w labels tc ''white'' notitle'
      write(u, '(a)') '#          (u 1:2:4 for Im chi0/w)'
    else
      write(u, '(a)') '# gnuplot: set view map; set size ratio -1; splot ''' // fmap // ''' u 1:2:(($3+$5)/2) w pm3d' // &
        '   (the average of the two planes; u 1:2:3 the first, u 1:2:5 the second)'
    end if
    call chi_header(u)
    if (nps == 1) then
      write(u, '(a)') '# columns: ' // names(1) // ' ' // names(2) // ' [1/Angstrom]  Re chi0  Im chi0/w  H K L [r.l.u.]'
    else
      write(u, '(a)') '# columns: ' // names(1) // ' ' // names(2) // ' [1/Angstrom]  Re chi0 (' // trim(ptag(1)) // &
        ')  Im chi0/w (' // trim(ptag(1)) // ')  Re chi0 (' // trim(ptag(2)) // ')  Im chi0/w (' // trim(ptag(2)) // &
        ')  H K L [r.l.u.] (of ' // trim(ptag(1)) // ')'
    end if
    do i = i0, i1
      do j = j0, j1
        call pair_point(i, j, ia, ib, ic, gcs, mab, x, y, hkl, v, near)
        if (nps == 1) then
          write(u, '(2f12.6,2es19.10e3,3f11.6)') x, y, v(:, 1), hkl
        else
          write(u, '(2f12.6,4es19.10e3,3f11.6)') x, y, v(:, 1), v(:, 2), hkl
        end if
      end do
      write(u, '(a)') ''
    end do
    close(u)

    ! ---- points, at their own H_a, H_b
    open(newunit=u, file=fpt, status='replace', action='write')
    write(u, '(a)') '# points of the input in the plane of ' // fmap // ', each at its own H_a, H_b (its third ' // &
      'coordinate is not used: e.g. Gamma at (0, 0))'
    write(u, '(a)') '# columns: ' // names(1) // ' ' // names(2) // ' [1/Angstrom]  label  H K L [r.l.u.]  ' // &
      'distance from the plane ' // trim(ptag(1)) // ' [1/Angstrom]'
    do p = 1, size(plab)
      q = matmul(bvec, pk(:, p))
      write(u, '(2f12.6,2x,a,3f11.6,f12.6)') mab(1, 1) * pk(ia, p) + mab(1, 2) * pk(ib, p), &
        mab(2, 1) * pk(ia, p) + mab(2, 2) * pk(ib, p), plab(p)%s, pk(:, p), dot_product(q, nrm) - on(1)
    end do
    close(u)

    ! ---- screen
    write(output_unit, '(a)') '  plane       ' // trim(ptag(1)) // ' ' // hkl_name(ia, ib, ic, cs(1)) // &
      pick(nps == 2, ' and ' // trim(ptag(2)) // ' ' // hkl_name(ia, ib, ic, cs(2)), '') // ', spanned by b' // &
      istr(ia) // ' and b' // istr(ib) // '; ' // names(1) // ' along (' // vtext(e1) // '), ' // names(2) // &
      ' along (' // vtext(e2) // ')'
    write(output_unit, '(a)') '  written     ' // fmap // ': ' // istr(i1 - i0 + 1) // ' x ' // istr(j1 - j0 + 1) // &
      ' points' // pick(nps == 2, ', both planes', '') // ', ' // pick(zone == 'first', 'first Brillouin zone', &
      'extended zones')
    write(output_unit, '(a)') '  written     ' // fbz // ': ' // istr(npoly) // ' zone boundaries' // &
      pick(nps == 2 .and. nz == 2, ' (the two planes have different zone sections)', '')
    if (size(plab) > 0) then
      write(output_unit, '(a)') '  written     ' // fpt // ':'
      do p = 1, size(plab)
        offnote = ''
        if (all(abs(pk(ic, p) - cs(1:nps)) > 1.0e-9_dp)) &
          offnote = '   (k' // istr(ic) // ' = ' // fstr(pk(ic, p), '(f8.4)') // ': drawn at its H_a, H_b)'
        write(output_unit, '(a)') '              ' // plab(p)%s // ' = (' // vtext(pk(:, p)) // ') r.l.u.  ->  ' // &
          names(1) // ', ' // names(2) // ' = ' // fstr(mab(1, 1) * pk(ia, p) + mab(1, 2) * pk(ib, p), '(f10.5)') // &
          ', ' // fstr(mab(2, 1) * pk(ia, p) + mab(2, 2) * pk(ib, p), '(f10.5)') // ' 1/Angstrom' // offnote
      end do
    end if
    do ip = 1, nps
      if (nps == 2) write(output_unit, '(a)') '  ' // trim(ptag(ip)) // ':'
      call max_header(output_unit, '  ', best(:, ip), bestx(:, :, ip), besth(:, :, ip), far(:, ip), farx(:, :, ip), &
                      farh(:, :, ip), names)
    end do
  end subroutine one_plane

  !> one point of the map, for each plane: indices i, j along the in-plane axes
  subroutine pair_point(i, j, ia, ib, ic, gcs, mab, x, y, hkl, v, near)
    integer, intent(in) :: i, j, ia, ib, ic, gcs(:)
    real(dp), intent(in) :: mab(2,2)
    real(dp), intent(out) :: x, y, hkl(3), v(2,2)
    logical, intent(out) :: near(2)
    integer :: g(3), da, db, ip
    real(dp) :: s, t
    s = real(i, dp) / nmesh(ia)
    t = real(j, dp) / nmesh(ib)
    x = mab(1, 1) * s + mab(1, 2) * t
    y = mab(2, 1) * s + mab(2, 2) * t
    g(ia) = modulo(i, nmesh(ia))
    g(ib) = modulo(j, nmesh(ib))
    da = min(g(ia), nmesh(ia) - g(ia))                      ! mesh steps from the nearest Gamma
    db = min(g(ib), nmesh(ib) - g(ib))
    v = nan
    near = .false.
    hkl(ia) = s
    hkl(ib) = t
    do ip = size(gcs), 1, -1                                ! ends with hkl of the first plane
      hkl(ic) = real(gcs(ip), dp) / nmesh(ic)
      g(ic) = modulo(gcs(ip), nmesh(ic))
      v(:, ip) = val(:, site(g))
      if (zone == 'first') then
        if (box_zone) then
          if (abs(s) > 0.5_dp + 1.0e-9_dp .or. abs(t) > 0.5_dp + 1.0e-9_dp) v(:, ip) = nan
        else
          if (.not. in_first_zone(matmul(bvec, hkl))) v(:, ip) = nan
        end if
      end if
      near(ip) = da <= 1 .and. db <= 1 .and. g(ic) == 0
    end do
  end subroutine pair_point

  !> the zone -1/2 ... 1/2 in crystal coordinates, in the plane coordinates (x, y) = M (H_a, H_b)
  pure function box_poly(mab) result(pl)
    real(dp), intent(in) :: mab(2,2)
    real(dp) :: pl(2,4)
    pl(:, 1) = matmul(mab, [-0.5_dp, -0.5_dp])
    pl(:, 2) = matmul(mab, [0.5_dp, -0.5_dp])
    pl(:, 3) = matmul(mab, [0.5_dp, 0.5_dp])
    pl(:, 4) = matmul(mab, [-0.5_dp, 0.5_dp])
  end function box_poly

  !> the same polygon (the same corners, in any order)
  logical function same_poly(a, b, tol)
    real(dp), intent(in) :: a(:,:), b(:,:), tol
    integer :: k
    same_poly = size(a, 2) == size(b, 2)
    if (.not. same_poly) return
    do k = 1, size(a, 2)
      if (minval(norm2(b - spread(a(:, k), 2, size(b, 2)), dim=1)) > tol) then
        same_poly = .false.
        return
      end if
    end do
  end function same_poly

  !> a signed number for a file name: +0.1667, -0.3333, up to 4 decimals
  function signed_text(x) result(s)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=24) :: buf
    integer :: n
    write(buf, '(f24.4)') abs(x)
    s = trim(adjustl(buf))
    n = len(s)
    do while (n > 1 .and. s(n:n) == '0')
      n = n - 1
    end do
    if (s(n:n) == '.') n = n - 1
    s = pick(x < 0.0_dp, '-', '+') // s(1:n)
  end function signed_text

  !> the same without the sign: 0.1667
  function plain_text(x) result(s)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    s = signed_text(abs(x))
    s = s(2:)
  end function plain_text

  !> the section of the zone around g0 with the plane  q . nrm = on, as a polygon
  !> in (q . e1, q . e2): |q - g0| <= |q - g0 - h| for every reciprocal-lattice
  !> vector h, i.e.  2 x (e1 . h) + 2 y (e2 . h) <= |h|^2 + 2 g0 . h - 2 on (nrm . h)
  subroutine zone_section(g0, e1, e2, nrm, on, poly)
    real(dp), intent(in) :: g0(3), e1(3), e2(3), nrm(3), on
    real(dp), allocatable, intent(out) :: poly(:,:)
    real(dp) :: big, cx, cy, h(3)
    integer :: k

    big = 3.0_dp * maxval(norm2(bvec, dim=1))
    cx = dot_product(g0, e1)
    cy = dot_product(g0, e2)
    poly = reshape([cx - big, cy - big, cx + big, cy - big, cx + big, cy + big, cx - big, cy + big], [2, 4])
    do k = 1, size(hv, 2)
      h = hv(:, k)
      call clip(poly, 2.0_dp * dot_product(e1, h), 2.0_dp * dot_product(e2, h), &
                dot_product(h, h) + 2.0_dp * dot_product(g0, h) - 2.0_dp * on * dot_product(nrm, h))
      if (size(poly, 2) < 3) exit
    end do
    if (size(poly, 2) >= 3) then                            ! a face touching the plane: no area
      if (abs(area(poly)) < 1.0e-9_dp * maxval(norm2(bvec, dim=1))**2) then
        deallocate(poly)
        allocate(poly(2, 0))
      end if
    end if
  end subroutine zone_section

  subroutine clip(poly, a, b, c)
    real(dp), allocatable, intent(inout) :: poly(:,:)
    real(dp), intent(in) :: a, b, c
    real(dp), allocatable :: out(:,:)
    real(dp) :: si, sj, tol
    integer :: n, i, j, k

    n = size(poly, 2)
    allocate(out(2, 2 * n + 2))
    tol = 1.0e-12_dp * (abs(c) + 1.0_dp)
    k = 0
    do i = 1, n
      j = merge(1, i + 1, i == n)
      si = a * poly(1, i) + b * poly(2, i) - c
      sj = a * poly(1, j) + b * poly(2, j) - c
      if (si <= tol) then
        k = k + 1
        out(:, k) = poly(:, i)
      end if
      if ((si < -tol .and. sj > tol) .or. (si > tol .and. sj < -tol)) then
        k = k + 1
        out(:, k) = poly(:, i) + si / (si - sj) * (poly(:, j) - poly(:, i))
      end if
    end do
    poly = out(:, 1:k)
  end subroutine clip

  pure real(dp) function area(poly)
    real(dp), intent(in) :: poly(:,:)
    integer :: i, j
    area = 0.0_dp
    do i = 1, size(poly, 2)
      j = merge(1, i + 1, i == size(poly, 2))
      area = area + 0.5_dp * (poly(1, i) * poly(2, j) - poly(1, j) * poly(2, i))
    end do
  end function area

  subroutine write_poly(u, poly)
    integer, intent(in) :: u
    real(dp), intent(in) :: poly(:,:)
    integer :: k
    do k = 1, size(poly, 2)
      write(u, '(2f12.6)') poly(:, k)
    end do
    write(u, '(2f12.6)') poly(:, 1)
    write(u, '(a)') ''
  end subroutine write_poly

  logical function in_first_zone(q)
    real(dp), intent(in) :: q(3)
    integer :: k
    real(dp) :: h2
    in_first_zone = .true.
    do k = 1, size(hv, 2)
      h2 = dot_product(hv(:, k), hv(:, k))
      if (2.0_dp * dot_product(q, hv(:, k)) > h2 * (1.0_dp + 1.0e-8_dp)) then
        in_first_zone = .false.
        return
      end if
    end do
  end function in_first_zone

  !> the largest value; among equal ones the one nearest to Gamma
  subroutine keep_max(v, qn, x, y, hkl, best, bestq, bestx, besth)
    real(dp), intent(in) :: v, qn, x, y, hkl(3)
    real(dp), intent(inout) :: best, bestq, bestx(2), besth(3)
    real(dp) :: tol
    tol = 1.0e-10_dp * max(abs(v), abs(best))
    if (v > best + tol .or. (abs(v - best) <= tol .and. qn < bestq - 1.0e-9_dp)) then
      best = v
      bestq = qn
      bestx = [x, y]
      besth = hkl
    end if
  end subroutine keep_max

  subroutine max_header(u, lead, best, bestx, besth, far, farx, farh, names)
    integer, intent(in) :: u
    character(len=*), intent(in) :: lead
    real(dp), intent(in) :: best(2), bestx(2,2), besth(3,2), far(2), farx(2,2), farh(3,2)
    character(len=2), intent(in) :: names(2)
    character(len=9), parameter :: what(2) = ['Re chi0  ', 'Im chi0/w']
    integer :: p
    do p = 1, 2
      if (best(p) > -huge(1.0_dp)) write(u, '(a)') lead // 'largest ' // trim(what(p)) // ': ' // &
        fstr(best(p), '(es13.5)') // ' at ' // names(1) // ', ' // names(2) // ' = ' // fstr(bestx(1, p), '(f9.5)') // &
        ', ' // fstr(bestx(2, p), '(f9.5)') // ' 1/Angstrom, H K L = ' // vtext(besth(:, p))
      if (far(p) > -huge(1.0_dp)) write(u, '(a)') lead // '   away from Gamma (> 1 mesh step): ' // &
        fstr(far(p), '(es13.5)') // ' at ' // names(1) // ', ' // names(2) // ' = ' // fstr(farx(1, p), '(f9.5)') // &
        ', ' // fstr(farx(2, p), '(f9.5)') // ' 1/Angstrom, H K L = ' // vtext(farh(:, p))
    end do
  end subroutine max_header

  !=============================================================================
  ! line (1D systems)
  !=============================================================================
  subroutine line_output()
    integer :: ia, i, i0, i1, u, k, m, m0, m1, p, g(3)
    real(dp) :: b(3), bl, smax, bh, h2, s, hkl(3), v(2), vmin(2), vmax(2), q(3), ql
    real(dp) :: best(2), bestq(2), besth(3,2), far(2), farq(2), farh(3,2), bestx(2,2), farx(2,2), xs(2)
    character(len=:), allocatable :: fmap, fbz, fpt

    ia = maxloc(nmesh, dim=1)
    b = bvec(:, ia)
    bl = norm2(b)
    smax = 0.5_dp                                           ! zone boundary along the line, in units of b
    do k = 1, merge(0, size(hv, 2), box_zone)
      bh = dot_product(b, hv(:, k))
      h2 = dot_product(hv(:, k), hv(:, k))
      if (bh > 1.0e-12_dp) smax = min(smax, h2 / (2.0_dp * bh))
    end do
    if (zone == 'first') then
      i1 = floor(smax * nmesh(ia) + 1.0e-6_dp)
    else
      i1 = floor(extent * nmesh(ia) + 1.0e-6_dp)
    end if
    i0 = -i1
    if (len(prefix) == 0) prefix = stem(chi_file) // '_line'
    fmap = prefix // '_map.dat'
    fbz = prefix // '_bz.dat'
    fpt = prefix // '_points.dat'

    best = -huge(1.0_dp)
    far = -huge(1.0_dp)
    bestq = huge(1.0_dp)
    farq = huge(1.0_dp)
    besth = 0.0_dp
    farh = 0.0_dp
    bestx = 0.0_dp
    farx = 0.0_dp
    vmin = huge(1.0_dp)
    vmax = -huge(1.0_dp)
    do i = i0, i1
      call line_point(i, ia, s, hkl, v, g)
      do p = 1, 2
        if (ieee_is_nan(v(p))) cycle
        vmin(p) = min(vmin(p), v(p))
        vmax(p) = max(vmax(p), v(p))
        call keep_max(v(p), abs(s * bl), s * bl, 0.0_dp, hkl, best(p), bestq(p), bestx(:, p), besth(:, p))
        if (min(g(ia), nmesh(ia) - g(ia)) > 1) &
          call keep_max(v(p), abs(s * bl), s * bl, 0.0_dp, hkl, far(p), farq(p), farx(:, p), farh(:, p))
      end do
    end do
    where (vmin > vmax)
      vmin = 0.0_dp
      vmax = 0.0_dp
    end where

    open(newunit=u, file=fmap, status='replace', action='write')
    write(u, '(a)') '# kind: line'
    write(u, '(a)') '# chi0 of ' // chi_file // ' along b' // istr(ia) // ' (q = H b' // istr(ia) // '), ' // &
      pick(zone == 'first', 'first Brillouin zone', 'H from ' // fstr(-extent, '(f8.3)') // ' to ' // fstr(extent, '(f8.3)'))
    write(u, '(a)') '# values of the mesh points (no interpolation), repeated with the periodicity of the reciprocal lattice'
    call cell_header(u)
    write(u, '(a)') '# q = signed distance from Gamma along b' // istr(ia) // ' = (' // vtext(b / bl) // &
      '); |b' // istr(ia) // '| = ' // fstr(bl, '(f10.6)') // ' 1/Angstrom; zone boundaries at q = ' // &
      fstr(smax * bl, '(f10.6)') // ' + n |b|'
    do p = 1, 2
      xs = [bestx(1, p), farx(1, p)]
      if (best(p) > -huge(1.0_dp)) write(u, '(a)') '# largest ' // trim(merge('Re chi0  ', 'Im chi0/w', p == 1)) // &
        ': ' // fstr(best(p), '(es13.5)') // ' at q = ' // fstr(xs(1), '(f9.5)') // ' 1/Angstrom, H K L = ' // &
        vtext(besth(:, p))
      if (far(p) > -huge(1.0_dp)) write(u, '(a)') '#    away from Gamma (> 1 mesh step): ' // fstr(far(p), '(es13.5)') // &
        ' at q = ' // fstr(xs(2), '(f9.5)') // ' 1/Angstrom, H K L = ' // vtext(farh(:, p))
    end do
    write(u, '(a)') '# gnuplot: plot ''' // fmap // ''' u 1:2 w lp, ''' // fbz // ''' u 1:2 w l lc ''black'' notitle' // &
      '   (u 1:3 for Im chi0/w)'
    call chi_header(u)
    write(u, '(a)') '# columns: q [1/Angstrom]  Re chi0  Im chi0/w  H K L [r.l.u.]'
    do i = i0, i1
      call line_point(i, ia, s, hkl, v, g)
      write(u, '(f12.6,2es19.10e3,3f11.6)') s * bl, v, hkl
    end do
    close(u)

    open(newunit=u, file=fbz, status='replace', action='write')
    write(u, '(a)') '# zone boundaries along the line, as vertical segments from the lowest to the highest value'
    write(u, '(a)') '# columns: q [1/Angstrom]  Re chi0  Im chi0/w   (plot u 1:2 with Re, u 1:3 with Im)'
    if (zone == 'first') then
      m0 = 0
      m1 = 0
    else
      m0 = -ceiling(extent) - 1
      m1 = ceiling(extent) + 1
    end if
    do m = m0, m1
      do k = -1, 1, 2
        ql = (m + k * smax) * bl
        if (zone == 'extended' .and. abs(ql) > extent * bl + 1.0e-9_dp) cycle
        write(u, '(f12.6,2es19.10e3)') ql, vmin
        write(u, '(f12.6,2es19.10e3)') ql, vmax
        write(u, '(a)') ''
      end do
    end do
    close(u)

    open(newunit=u, file=fpt, status='replace', action='write')
    write(u, '(a)') '# points of the input on the line of ' // fmap // '; a point off the line is projected onto it'
    write(u, '(a)') '# columns: q [1/Angstrom]  label  H K L [r.l.u.]  distance from the line [1/Angstrom]'
    do p = 1, size(plab)
      q = matmul(bvec, pk(:, p))
      ql = dot_product(q, b) / bl
      write(u, '(f12.6,2x,a,3f11.6,f12.6)') ql, plab(p)%s, pk(:, p), norm2(q - ql * b / bl)
    end do
    close(u)

    write(output_unit, '(a)') '  line        along b' // istr(ia) // ' = (' // vtext(b) // ') 1/Angstrom, |b| = ' // &
      fstr(bl, '(f10.6)') // ' 1/Angstrom; zone boundary at q = ' // fstr(smax * bl, '(f10.6)') // ' 1/Angstrom (H = ' // &
      fstr(smax, '(f8.5)') // ')'
    write(output_unit, '(a)') '  written     ' // fmap // ': ' // istr(i1 - i0 + 1) // ' points'
    write(output_unit, '(a)') '  written     ' // fbz // ', ' // fpt
    do p = 1, size(plab)
      q = matmul(bvec, pk(:, p))
      write(output_unit, '(a)') '              ' // plab(p)%s // ' = (' // vtext(pk(:, p)) // ') r.l.u.  ->  q = ' // &
        fstr(dot_product(q, b) / bl, '(f10.5)') // ' 1/Angstrom'
    end do
    do p = 1, 2
      if (best(p) > -huge(1.0_dp)) write(output_unit, '(a)') '  largest ' // &
        trim(merge('Re chi0  ', 'Im chi0/w', p == 1)) // ': ' // fstr(best(p), '(es13.5)') // ' at q = ' // &
        fstr(bestx(1, p), '(f9.5)') // ' 1/Angstrom, H K L = ' // vtext(besth(:, p))
      if (far(p) > -huge(1.0_dp)) write(output_unit, '(a)') '     away from Gamma (> 1 mesh step): ' // &
        fstr(far(p), '(es13.5)') // ' at q = ' // fstr(farx(1, p), '(f9.5)') // ' 1/Angstrom, H K L = ' // vtext(farh(:, p))
    end do
  end subroutine line_output

  subroutine line_point(i, ia, s, hkl, v, g)
    integer, intent(in) :: i, ia
    real(dp), intent(out) :: s, hkl(3), v(2)
    integer, intent(out) :: g(3)
    s = real(i, dp) / nmesh(ia)
    hkl = 0.0_dp
    hkl(ia) = s
    g = 0
    g(ia) = modulo(i, nmesh(ia))
    v = val(:, site(g))
  end subroutine line_point

  !=============================================================================
  ! headers
  !=============================================================================
  subroutine cell_header(u)
    integer, intent(in) :: u
    write(u, '(a)') '# cell ' // pick(len(nscf_out) > 0, 'from ' // nscf_out, 'of the input') // ' (Angstrom): a1 = (' // &
      vtext(avec(:, 1)) // '), a2 = (' // vtext(avec(:, 2)) // '), a3 = (' // vtext(avec(:, 3)) // ')'
    write(u, '(a)') '# reciprocal axes (1/Angstrom, 2 pi included): b1 = (' // vtext(bvec(:, 1)) // '), b2 = (' // &
      vtext(bvec(:, 2)) // '), b3 = (' // vtext(bvec(:, 3)) // ')'
  end subroutine cell_header

  !> the header of the chi file, for the record
  subroutine chi_header(u)
    integer, intent(in) :: u
    integer :: k
    write(u, '(a)') '# header of ' // chi_file // ':'
    do k = 1, size(chi_head)
      if (index(chi_head(k)%s, '# columns') == 1) cycle
      write(u, '(a)') '# | ' // chi_head(k)%s(2:)
    end do
  end subroutine chi_header

  !> e.g. (H K 0), (H 0.5 L)
  function hkl_name(ia, ib, ic, c) result(s)
    integer, intent(in) :: ia, ib, ic
    real(dp), intent(in) :: c
    character(len=:), allocatable :: s
    character(len=1), parameter :: letter(3) = ['H', 'K', 'L']
    character(len=12) :: part(3), buf
    integer :: n
    part(ia) = letter(ia)
    part(ib) = letter(ib)
    write(buf, '(f12.4)') c
    buf = adjustl(buf)
    n = len_trim(buf)
    do while (n > 1 .and. buf(n:n) == '0')
      n = n - 1
    end do
    if (buf(n:n) == '.') n = n - 1
    part(ic) = buf(1:n)
    s = '(' // trim(part(1)) // ' ' // trim(part(2)) // ' ' // trim(part(3)) // ')'
  end function hkl_name

  !=============================================================================
  ! small helpers
  !=============================================================================
  subroutine read_line(u, line, ios)
    integer, intent(in) :: u
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out) :: ios
    character(len=512) :: buf
    integer :: n
    line = ''
    do
      read(u, '(a)', advance='no', iostat=ios, size=n) buf
      line = line // buf(1:n)
      if (ios /= 0) exit
    end do
    if (is_iostat_eor(ios)) ios = 0
    if (ios /= 0 .and. len(line) > 0) ios = 0              ! last line without a newline
    n = len(line)
    if (n > 0) then
      if (line(n:n) == char(13)) line = line(1:n - 1)      ! Windows line ending
    end if
  end subroutine read_line

  function strip_comment(line) result(s)
    character(len=*), intent(in) :: line
    character(len=:), allocatable :: s
    integer :: i
    character :: q
    q = ' '
    do i = 1, len(line)
      if (q == ' ') then
        if (line(i:i) == '!') then
          s = line(1:i - 1)
          return
        end if
        if (line(i:i) == '''' .or. line(i:i) == '"') q = line(i:i)
      else if (line(i:i) == q) then
        q = ' '
      end if
    end do
    s = line
  end function strip_comment

  integer function count_tokens(line)
    character(len=*), intent(in) :: line
    integer :: i
    logical :: inword
    count_tokens = 0
    inword = .false.
    do i = 1, len(line)
      if (line(i:i) == ' ' .or. line(i:i) == char(9) .or. line(i:i) == ',') then
        inword = .false.
      else if (.not. inword) then
        inword = .true.
        count_tokens = count_tokens + 1
      end if
    end do
  end function count_tokens

  function lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i
    t = s
    do i = 1, len(s)
      if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') t(i:i) = achar(iachar(s(i:i)) + 32)
    end do
  end function lower

  !> the chi file name without its folder and .dat
  function stem(fname) result(s)
    character(len=*), intent(in) :: fname
    character(len=:), allocatable :: s
    integer :: k
    s = fname
    k = index(s, '/', back=.true.)
    if (k > 0) s = s(k + 1:)
    k = len(s)
    if (k > 4) then
      if (lower(s(k - 3:k)) == '.dat') s = s(1:k - 4)
    end if
  end function stem

  function istr(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=20) :: buf
    write(buf, '(i0)') i
    s = trim(buf)
  end function istr

  function fstr(x, fmt) result(s)
    real(dp), intent(in) :: x
    character(len=*), intent(in) :: fmt
    character(len=:), allocatable :: s
    character(len=40) :: buf
    write(buf, fmt) x
    s = trim(adjustl(buf))
  end function fstr

  function vtext(v) result(s)
    real(dp), intent(in) :: v(3)
    character(len=:), allocatable :: s
    s = fstr(v(1), '(f10.5)') // ' ' // fstr(v(2), '(f10.5)') // ' ' // fstr(v(3), '(f10.5)')
  end function vtext

  function pick(cond, a, b) result(s)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: a, b
    character(len=:), allocatable :: s
    if (cond) then
      s = a
    else
      s = b
    end if
  end function pick

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
  end function cross

  pure real(dp) function det3(a)
    real(dp), intent(in) :: a(3,3)
    det3 = dot_product(a(:, 1), cross(a(:, 2), a(:, 3)))
  end function det3

  pure function inv3(a) result(b)
    real(dp), intent(in) :: a(3,3)
    real(dp) :: b(3,3)
    b(1, :) = cross(a(:, 2), a(:, 3))
    b(2, :) = cross(a(:, 3), a(:, 1))
    b(3, :) = cross(a(:, 1), a(:, 2))
    b = b / det3(a)
  end function inv3

  !> (s, t) with  m (s, t) = r
  pure function solve2(m, r) result(st)
    real(dp), intent(in) :: m(2,2), r(2)
    real(dp) :: st(2), d
    d = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
    st = [m(2, 2) * r(1) - m(1, 2) * r(2), m(1, 1) * r(2) - m(2, 1) * r(1)] / d
  end function solve2

end program chi_bz
