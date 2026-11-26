!===============================================================================
! lindhard_chi.f90  ->  bin/lindhard_chi.x
!
! Lindhard susceptibility chi0(q) from band energies on a full
! Gamma-centred k-mesh, parallelised over q with MPI. The band energies are read
! directly from Quantum ESPRESSO output, or from a ready table. The real and the
! imaginary part are computed together and written to one file,
! chi_<E_F>_<T>.dat (e.g. chi_6.4521_300.dat). The post-processing tool
! bin/chi_bz.x (src/chi_bz.f90) puts the result on the Brillouin zone;
! examples/plot_bz.py draws it (2D maps, and 1D cross sections of them).
!
!   Re chi0(q)          = -(1/N) sum_k sum_nm [f(e_nk) - f(e_m,k+q)] / (e_nk - e_m,k+q)
!   Im chi0(q,w)/w|w->0 =  (pi/N) sum_k sum_nm delta(e_nk - E_F) delta(e_m,k+q - E_F)
!
!   * constant matrix elements; f = Fermi-Dirac at temperature T; for
!     |e_nk - e_m,k+q| < 1e-6 eV the real-part term is its exact limit
!     -df/de = f(1-f)/(kB T), so the intraband terms near q = 0 are kept
!   * delta = normalised Gaussian of width sigma; Im chi0/w / pi is the
!     Fermi-surface nesting function
!   * the sum over k runs over every line of the table (step1.dat), the +1/2
!     zone-edge lines included, N = number of lines
!   * per cell and per spin (x2 for spin-degenerate bands without spin-orbit)
!   * the q-points are the lines of the table; k+q is located on the mesh by
!     integer index arithmetic, so the cost is N^2 m^2 with no k-point search
!
! Run:   mpirun -np 8 bin/lindhard_chi.x chi.in      (or  ... < chi.in)
!
!===============================================================================

module lindhard_utils
  use, intrinsic :: iso_fortran_env, only: real64, iostat_end, iostat_eor
  implicit none
  private
  public :: dp, KB, DEG, string_t, read_line, strip_comment, count_tokens, lower, &
            parse_yes_no, site_index, site_vector, istr, fstr, pick, int_after, real_after, &
            vec_after, sort_small, sorted_copy, push_string

  integer, parameter :: dp = real64
  real(dp), parameter :: KB = 8.617333262e-5_dp   ! Boltzmann constant (eV/K)
  real(dp), parameter :: DEG = 1.0e-6_dp          ! eV: smaller energy differences use the degenerate limit

  type :: string_t
    character(len=:), allocatable :: s
  end type string_t

contains

  subroutine read_line(u, line, ios)
    integer, intent(in) :: u
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out) :: ios
    character(len=512) :: buf
    integer :: n, i

    line = ''
    do
      read(u, '(a)', advance='no', iostat=ios, size=n) buf
      if (ios > 0) return
      line = line // buf(1:n)
      if (ios == iostat_eor) exit
      if (ios == iostat_end) then
        if (len(line) == 0) return
        exit
      end if
    end do
    ios = 0
    n = len(line)
    if (n > 0) then
      if (line(n:n) == achar(13)) line = line(1:n-1)
    end if
    do i = 1, len(line)
      if (line(i:i) == achar(9)) line(i:i) = ' '
    end do
  end subroutine read_line

  function strip_comment(line) result(s)
    character(len=*), intent(in) :: line
    character(len=:), allocatable :: s
    character(len=1) :: q
    integer :: p, i

    p = 0
    q = ' '
    do i = 1, len(line)                              ! first '!' outside quotes
      if (q /= ' ') then
        if (line(i:i) == q) q = ' '
      else if (line(i:i) == '''' .or. line(i:i) == '"') then
        q = line(i:i)
      else if (line(i:i) == '!') then
        p = i
        exit
      end if
    end do
    if (p > 0) then
      s = trim(adjustl(line(1:p-1)))
    else
      s = trim(adjustl(line))
    end if
  end function strip_comment

  pure integer function count_tokens(line) result(n)
    character(len=*), intent(in) :: line
    integer :: i
    logical :: in_token

    n = 0
    in_token = .false.
    do i = 1, len(line)
      if (line(i:i) == ' ' .or. line(i:i) == ',') then
        in_token = .false.
      else if (.not. in_token) then
        n = n + 1
        in_token = .true.
      end if
    end do
  end function count_tokens

  pure function lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i

    t = s
    do i = 1, len(s)
      if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') t(i:i) = achar(iachar(s(i:i)) + 32)
    end do
  end function lower

  subroutine parse_yes_no(s, val, ok)
    character(len=*), intent(in) :: s
    logical, intent(out) :: val, ok
    character(len=:), allocatable :: t

    t = trim(adjustl(lower(s)))
    ok = .true.
    select case (t)
    case ('yes', 'y', 'true', 't', '.true.', '1')
      val = .true.
    case ('no', 'n', 'false', 'f', '.false.', '0')
      val = .false.
    case default
      val = .false.
      ok = .false.
    end select
  end subroutine parse_yes_no

  pure integer function site_index(i, n) result(s)
    integer, intent(in) :: i(3), n(3)
    s = 1 + i(1) + n(1) * (i(2) + n(2) * i(3))
  end function site_index

  pure function site_vector(s, n) result(i)
    integer, intent(in) :: s, n(3)
    integer :: i(3), t
    t = s - 1
    i(1) = mod(t, n(1))
    t = t / n(1)
    i(2) = mod(t, n(2))
    i(3) = t / n(2)
  end function site_vector

  function istr(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=32) :: b
    write(b, '(i0)') i
    s = trim(b)
  end function istr

  function fstr(x, fmt) result(s)
    real(dp), intent(in) :: x
    character(len=*), intent(in) :: fmt
    character(len=:), allocatable :: s
    character(len=64) :: b
    write(b, fmt) x
    s = trim(adjustl(b))
  end function fstr

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

  subroutine int_after(line, marker, n, ios)
    character(len=*), intent(in) :: line, marker
    integer, intent(out) :: n, ios
    integer :: p
    n = 0
    p = index(line, marker)
    ios = 1
    if (p > 0) read(line(p + len(marker):), *, iostat=ios) n
  end subroutine int_after

  subroutine real_after(line, marker, x, ios)
    character(len=*), intent(in) :: line, marker
    real(dp), intent(out) :: x
    integer, intent(out) :: ios
    integer :: p
    x = 0.0_dp
    p = index(line, marker)
    ios = 1
    if (p > 0) read(line(p + len(marker):), *, iostat=ios) x
  end subroutine real_after

  subroutine vec_after(line, marker, v, ios)
    character(len=*), intent(in) :: line, marker
    real(dp), intent(out) :: v(:)
    integer, intent(out) :: ios
    integer :: p, q
    v = 0.0_dp
    ios = 1
    p = index(line, marker)
    if (p == 0) return
    p = p + len(marker)
    q = index(line(p:), ')')
    if (q > 0) then
      read(line(p:p+q-2), *, iostat=ios) v
    else
      read(line(p:), *, iostat=ios) v
    end if
  end subroutine vec_after

  pure function sort_small(a) result(b)
    real(dp), intent(in) :: a(:)
    real(dp) :: b(size(a)), x
    integer :: i, j
    b = a
    do i = 2, size(b)
      x = b(i)
      j = i - 1
      do while (j >= 1)
        if (b(j) <= x) exit
        b(j + 1) = b(j)
        j = j - 1
      end do
      b(j + 1) = x
    end do
  end function sort_small

  subroutine push_string(arr, s)
    type(string_t), allocatable, intent(inout) :: arr(:)
    character(len=*), intent(in) :: s
    type(string_t), allocatable :: tmp(:)
    integer :: i, n
    n = 0
    if (allocated(arr)) n = size(arr)
    allocate(tmp(n + 1))
    do i = 1, n
      call move_alloc(arr(i)%s, tmp(i)%s)
    end do
    tmp(n + 1)%s = s
    call move_alloc(tmp, arr)
  end subroutine push_string

  function sorted_copy(a) result(b)
    real(dp), intent(in) :: a(:)
    real(dp), allocatable :: b(:)
    integer :: n, i, last
    real(dp) :: x
    b = a
    n = size(b)
    do i = n / 2, 1, -1
      call sift(i, n)
    end do
    do last = n, 2, -1
      x = b(1)
      b(1) = b(last)
      b(last) = x
      call sift(1, last - 1)
    end do
  contains
    subroutine sift(start, nend)
      integer, intent(in) :: start, nend
      integer :: root, child
      real(dp) :: y
      root = start
      do while (2 * root <= nend)
        child = 2 * root
        if (child < nend) then
          if (b(child) < b(child + 1)) child = child + 1
        end if
        if (b(root) >= b(child)) return
        y = b(root)
        b(root) = b(child)
        b(child) = y
        root = child
      end do
    end subroutine sift
  end function sorted_copy


end module lindhard_utils


!===============================================================================
!> The input namelist &lindhard
module lindhard_input
  use lindhard_utils, only: dp
  implicit none
  private
  public :: main_t, read_main, MAXSEL, UNSET

  integer, parameter :: MAXSEL = 1000                ! kpt_index(1..MAXSEL)
  real(dp), parameter :: UNSET = -huge(1.0_dp)

  type :: main_t
    character(len=512) :: nscf_out = ' ', bands_file = ' ', band_file = ' '
    integer :: band_first = 0, band_last = 0, spin = 0, mesh(3) = 0
    integer :: kpt_index(MAXSEL) = 0
    real(dp) :: ef = UNSET, temp = UNSET, sigma = 0.02_dp, period(3) = 1.0_dp
  end type main_t

contains

  subroutine read_main(u, st, ios, msg)
    integer, intent(in) :: u
    type(main_t), intent(out) :: st
    integer, intent(out) :: ios
    character(len=*), intent(out) :: msg
    character(len=512) :: nscf_out, bands_file, band_file
    integer :: band_first, band_last, spin, mesh(3), kpt_index(MAXSEL)
    real(dp) :: ef, temp, sigma, period(3)
    namelist /lindhard/ nscf_out, bands_file, band_file, temp, ef, sigma, band_first, band_last, &
                        spin, kpt_index, mesh, period

    nscf_out = st%nscf_out
    bands_file = st%bands_file
    band_file = st%band_file
    band_first = st%band_first
    band_last = st%band_last
    spin = st%spin
    mesh = st%mesh
    kpt_index = st%kpt_index
    ef = st%ef
    temp = st%temp
    sigma = st%sigma
    period = st%period
    msg = ''
    read(u, nml=lindhard, iostat=ios, iomsg=msg)
    if (ios /= 0) return
    st%nscf_out = nscf_out
    st%bands_file = bands_file
    st%band_file = band_file
    st%band_first = band_first
    st%band_last = band_last
    st%spin = spin
    st%mesh = mesh
    st%kpt_index = kpt_index
    st%ef = ef
    st%temp = temp
    st%sigma = sigma
    st%period = period
  end subroutine read_main


end module lindhard_input


!===============================================================================
!> Readers for pw.x output and the bands.x gnu file
module qe_reader
  use lindhard_utils
  implicit none
  private
  public :: pw_data_t, read_pw_output, read_bands_gnu

  type :: pw_data_t
    integer :: nks = 0, nbnd = 0, nstored = 0
    real(dp), allocatable :: kcry(:,:), kcart(:,:), e(:,:)    ! (3,nks), (3,nks), (nbnd,nks)
    real(dp) :: at(3,3) = 0.0_dp, ef = 0.0_dp, ef_spin(2) = 0.0_dp, header_mismatch = 0.0_dp
    logical :: have_at = .false., have_cart = .false., have_cryst = .false., have_e = .false.
    logical :: have_ef = .false., have_ef_spin = .false., lsda = .false.
  end type pw_data_t

contains

  subroutine read_pw_output(fname, spin, d, err)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: spin
    type(pw_data_t), intent(out) :: d
    character(len=:), allocatable, intent(out) :: err
    character(len=:), allocatable :: line, t
    integer :: u, ios, i, j, ik, cur_spin, nread, nv, nline, p
    real(dp) :: v(8), kh(3)
    logical :: in_bands, store

    err = ''
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      err = 'cannot open the pw.x output "' // fname // '"'
      return
    end if
    in_bands = .false.
    cur_spin = 1
    ik = 0
    nline = 0
    do
      call read_line(u, line, ios)
      if (ios /= 0) exit
      nline = nline + 1
      t = trim(adjustl(line))

      if (index(line, 'number of Kohn-Sham states=') > 0) then
        call int_after(line, '=', d%nbnd, ios)

      else if (index(line, 'crystal axes: (cart. coord. in units of alat)') > 0) then
        do i = 1, 3
          call read_line(u, line, ios)
          nline = nline + 1
          if (ios == 0) call vec_after(line, '= (', d%at(:, i), ios)
          if (ios /= 0) then
            err = fname // ', line ' // istr(nline) // ': cannot read the crystal axes'
            return
          end if
        end do
        d%have_at = .true.

      else if (index(line, 'number of k points=') > 0) then
        call int_after(line, '=', d%nks, ios)
        if (allocated(d%kcart)) deallocate(d%kcart, d%kcry)
        allocate(d%kcart(3, d%nks), d%kcry(3, d%nks))
        d%have_cart = .false.
        d%have_cryst = .false.

      else if (t == 'cart. coord. in units 2pi/alat' .or. t == 'cryst. coord.') then
        if (d%nks < 1) cycle
        do i = 1, d%nks
          call read_line(u, line, ios)
          nline = nline + 1
          if (ios == 0 .and. index(line, 'k(') == 0) ios = 1
          if (ios == 0) then
            if (t == 'cryst. coord.') then
              call vec_after(line, '= (', d%kcry(:, i), ios)
            else
              call vec_after(line, '= (', d%kcart(:, i), ios)
            end if
          end if
          if (ios /= 0) then
            err = fname // ', line ' // istr(nline) // ': the k-point list (' // t // ') ends after ' // &
                  istr(i - 1) // ' of ' // istr(d%nks) // ' k-points'
            return
          end if
        end do
        if (t == 'cryst. coord.') then
          d%have_cryst = .true.
        else
          d%have_cart = .true.
        end if

      else if (index(line, 'End of band structure calculation') > 0 .or. &
               index(line, 'End of self-consistent calculation') > 0) then
        in_bands = .true.
        cur_spin = 1
        ik = 0
        d%lsda = .false.
        d%nstored = 0
        d%header_mismatch = 0.0_dp
        if (allocated(d%e)) deallocate(d%e)
        allocate(d%e(max(d%nbnd, 1), max(d%nks, 1)))

      else if (in_bands .and. index(line, '------ SPIN UP') > 0) then
        d%lsda = .true.
        cur_spin = 1
        ik = 0

      else if (in_bands .and. index(line, '------ SPIN DOWN') > 0) then
        cur_spin = 2
        ik = 0

      else if (in_bands .and. (index(line, 'bands (ev)') > 0 .or. index(line, 'band energies (ev)') > 0)) then
        ik = ik + 1
        if (ik > d%nks) then
          err = fname // ', line ' // istr(nline) // ': more band-energy blocks than k-points (' // &
                istr(d%nks) // ')'
          return
        end if
        store = (.not. d%lsda) .or. cur_spin == spin
        p = index(line, 'k =')
        if (p > 0 .and. d%have_cart) then                ! k of the block (cartesian, 4 decimals)
          read(line(p+3:min(len(line), p+23)), '(3f7.4)', iostat=ios) kh
          if (ios == 0) d%header_mismatch = max(d%header_mismatch, maxval(abs(kh - d%kcart(:, ik))))
        end if
        nread = 0
        do while (nread < d%nbnd)
          call read_line(u, line, ios)
          nline = nline + 1
          if (ios /= 0) then
            err = fname // ': the band energies of k-point ' // istr(ik) // ' are incomplete'
            return
          end if
          if (len_trim(line) == 0) cycle
          nv = min(8, d%nbnd - nread)
          read(line, '(2x,8f9.4)', iostat=ios) (v(j), j = 1, nv)
          if (ios /= 0) then
            err = fname // ', line ' // istr(nline) // ': cannot read band energies'
            return
          end if
          if (store) d%e(nread+1:nread+nv, ik) = v(1:nv)
          nread = nread + nv
        end do
        if (store) d%nstored = d%nstored + 1

      else if (index(line, 'the Fermi energy is') > 0) then
        call real_after(line, 'the Fermi energy is', d%ef, ios)
        d%have_ef = (ios == 0)

      else if (index(line, 'the spin up/dw Fermi energies are') > 0) then
        p = index(line, 'are') + 3
        read(line(p:), *, iostat=ios) d%ef_spin
        d%have_ef_spin = (ios == 0)
      end if
    end do
    close(u)

    if (d%nks < 1) then
      err = fname // ' is not a pw.x output ("number of k points=" not found)'
      return
    end if
    if (d%nbnd < 1) then
      err = fname // ': "number of Kohn-Sham states=" not found'
      return
    end if
    d%have_e = (d%nstored == d%nks)
    if (.not. d%have_cryst .and. d%have_cart .and. d%have_at) then
      do i = 1, d%nks                                   ! crystal = k . a_i  (k in 2pi/alat, a in alat)
        d%kcry(:, i) = matmul(d%kcart(:, i), d%at)
      end do
    end if
  end subroutine read_pw_output

  !> bands.x gnu file: per band one block of 'x  e' lines, blocks separated by a blank line
  subroutine read_bands_gnu(fname, x, e, nk, nb, err)
    character(len=*), intent(in) :: fname
    real(dp), allocatable, intent(out) :: x(:), e(:,:)
    integer, intent(out) :: nk, nb
    character(len=:), allocatable, intent(out) :: err
    character(len=:), allocatable :: line
    integer :: u, ios, cur, ib, i, nline
    real(dp) :: xv, ev

    err = ''
    nk = 0
    nb = 0
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      err = 'cannot open the bands file "' // fname // '"'
      return
    end if
    cur = 0
    do                                                   ! pass 1: block sizes
      call read_line(u, line, ios)
      if (ios == 0 .and. len_trim(line) > 0) then
        cur = cur + 1
        cycle
      end if
      if (cur > 0) then
        nb = nb + 1
        if (nb == 1) nk = cur
        if (cur /= nk) then
          err = fname // ': band ' // istr(nb) // ' has ' // istr(cur) // ' k-points, band 1 has ' // istr(nk)
          return
        end if
        cur = 0
      end if
      if (ios /= 0) exit
    end do
    if (nb == 0) then
      err = fname // ' contains no band data (bands.x writes it only when max_lines in PP/src/bands.f90 ' // &
            'is larger than the number of k-points)'
      return
    end if
    allocate(x(nk), e(nb, nk))
    rewind(u)
    ib = 1
    i = 0
    nline = 0
    do                                                   ! pass 2: values
      call read_line(u, line, ios)
      if (ios /= 0) exit
      nline = nline + 1
      if (len_trim(line) == 0) then
        if (i > 0) ib = ib + 1
        i = 0
        cycle
      end if
      i = i + 1
      read(line, *, iostat=ios) xv, ev
      if (ios /= 0) then
        err = fname // ', line ' // istr(nline) // ': two numbers (x, energy) expected'
        return
      end if
      e(ib, i) = ev
      if (ib == 1) x(i) = xv
    end do
    close(u)
  end subroutine read_bands_gnu

end module qe_reader

!the main program!
!===============================================================================
program lindhard_chi
  use mpi
  use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit
  use lindhard_utils
  use lindhard_input
  use qe_reader
  implicit none

  real(dp), parameter :: PI = 3.14159265358979323846_dp

  ! MPI
  integer :: rank, nproc, ierr
  logical :: is_root

  ! settings (rank 0)
  character(len=:), allocatable :: input_name, table_file, nscf_out, bands_file, source_note
  logical :: qe_route, ef_given, temp_given, mesh_given
  integer :: nk_limit, col_first, col_last, spin, band_first, band_last
  real(dp) :: per(3)

  integer :: nk, m, nsel, nt, nmesh(3)
  real(dp) :: ef, temp, kt, sigma, dk(3)

  real(dp), allocatable :: kpt(:,:), eig(:,:)        ! (3, nk), (m, nk)
  integer, allocatable :: band_label(:)              ! (m)   band number of each energy column
  integer, allocatable :: site_of(:)                 ! (nk)  mesh site of each k-point of the table
  integer, allocatable :: point_of_site(:)           ! (nt)  k-point representing each site
  integer, allocatable :: rows_of_site(:)            ! (nt)  number of lines of the table on each site
  integer, allocatable :: sel_index(:)               ! (nsel) selected k-points (table numbering)

  real(dp), allocatable :: e(:,:), f(:,:), w(:,:)    ! (m, nt) energies, Fermi-Dirac, -df/de
  real(dp), allocatable :: gdel(:)                   ! (nt)  sum over bands of delta(e - E_F)
  integer, allocatable :: sel_site(:), sel_of_site(:)
  real(dp), allocatable :: wsite(:)                  ! (nt)  weight of each site in the sum over k
  real(dp) :: wsum                                   ! sum of the weights (N of the 1/N)
  real(dp), allocatable :: chi_re(:), con_re(:,:)    ! Re chi0(q) (nt), its k-resolved parts (nt, nsel)
  real(dp), allocatable :: chi_im(:), con_im(:,:)    ! Im chi0(q,w)/w at w -> 0, the same
  real(dp) :: t_start, t_compute

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  is_root = (rank == 0)
  t_start = MPI_Wtime()

  if (is_root) then
    write(output_unit, '(a,i0,a)') 'lindhard_chi: static Lindhard susceptibility (Fortran/MPI, ', &
                                   nproc, ' ranks)'
    call read_input()
    if (qe_route) call prepare_from_qe()
    call read_table()
    call map_to_mesh()
    call report_setup()
  end if
  call share_setup()
  call compute_chi()
  if (is_root) then
    call write_results()
    write(output_unit, '(a,a,a)') 'finished in ', fstr(MPI_Wtime() - t_start, '(f12.2)'), ' s'
  end if
  call MPI_Finalize(ierr)

contains

  !-----------------------------------------------------------------------------
  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    integer :: ie
    write(output_unit, '(a)') 'ERROR: ' // msg          ! in the log of standard output ...
    flush(output_unit)
    write(error_unit, '(a)') 'ERROR: ' // msg           ! ... and on the screen
    flush(error_unit)
    call MPI_Abort(MPI_COMM_WORLD, 1, ie)
    stop 1
  end subroutine fatal

  !-----------------------------------------------------------------------------
  !> Input from the file named on the command line, or from standard input: the
  !> namelist &lindhard, or the line-by-line chi.in of old version of the code.
  subroutine read_input()
    type(string_t), allocatable :: lines(:)
    character(len=:), allocatable :: s
    integer :: u, ios, n

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
    allocate(lines(0))
    do
      call read_line(u, s, ios)
      if (ios /= 0) exit
      s = strip_comment(s)
      if (len(s) > 0) call push_string(lines, s)
    end do
    if (u /= input_unit) close(u)
    if (size(lines) == 0) call fatal('the input ' // input_name // ' is empty')

    ! defaults
    nscf_out = ''
    bands_file = ''
    table_file = ''
    source_note = ''
    spin = 0
    band_first = 0
    band_last = 0
    nk_limit = -1
    col_first = 1
    col_last = 0
    per = 1.0_dp
    mesh_given = .false.
    nmesh = 0
    ef_given = .false.
    temp_given = .false.
    sigma = 0.02_dp
    allocate(sel_index(0))

    if (lines(1)%s(1:1) == '&') then
      call parse_namelists(lines)
    else
      call parse_legacy(lines)
    end if
  end subroutine read_input




  !-----------------------------------------------------------------------------
  !> the namelist &lindhard
  subroutine parse_namelists(lines)
    type(string_t), intent(in) :: lines(:)
    type(main_t) :: st
    character(len=512) :: msg
    integer :: su, ios, i

    open(newunit=su, status='scratch', action='readwrite', form='formatted')
    do i = 1, size(lines)
      write(su, '(a)') lines(i)%s
    end do
    rewind(su)
    call read_main(su, st, ios, msg)
    if (ios /= 0) call fatal('namelist &lindhard in ' // input_name // ': ' // trim(msg))
    close(su)

    ! ---- &lindhard
    nscf_out = trim(st%nscf_out)
    bands_file = trim(st%bands_file)
    table_file = trim(st%band_file)
    qe_route = len(nscf_out) > 0
    if (qe_route .eqv. len(table_file) > 0) call fatal('give either nscf_out (pw.x output) or band_file ' // &
                                                       '(a k1 k2 k3 e_1 ... table)')
    if (len(bands_file) > 0 .and. .not. qe_route) call fatal('bands_file is used together with nscf_out')
    temp_given = (st%temp > 0.5_dp * UNSET)               ! UNSET = -huge: not given
    if (.not. temp_given) call fatal('temp (K) is needed')
    if (st%temp <= 0.0_dp) call fatal('temp must be > 0 K')
    temp = st%temp
    sigma = st%sigma
    if (sigma <= 0.0_dp) call fatal('sigma (eV) must be > 0')
    ef_given = (st%ef > 0.5_dp * UNSET)
    ef = st%ef
    if (.not. qe_route .and. .not. ef_given) call fatal('ef (eV) must be given with band_file')
    if (st%spin < 0 .or. st%spin > 2) call fatal('spin must be 1 (up) or 2 (down)')
    spin = st%spin
    if (st%band_first < 0 .or. st%band_last < 0) call fatal('band_first and band_last must be >= 1')
    if (st%band_last > 0 .and. st%band_first > st%band_last) call fatal('band_first > band_last')
    band_first = st%band_first
    band_last = st%band_last
    if (any(st%mesh /= 0)) then
      if (any(st%mesh < 1)) call fatal('mesh: three integers >= 1')
      nmesh = st%mesh
      mesh_given = .true.
    end if
    sel_index = pack(st%kpt_index, st%kpt_index > 0)    ! kpt_index(1), kpt_index(2), ... in this order
    if (qe_route) then
      per = 1.0_dp
    else
      per = st%period
      col_first = max(1, band_first)
      col_last = band_last
    end if

  end subroutine parse_namelists



  !-----------------------------------------------------------------------------
  !> chi.in of old version, band file, nk, m, E_F, T, max_k, min_k, debug,
  !> then optional lines (N1 N2 N3, or one selected k-point index per line)
  subroutine parse_legacy(lines)
    type(string_t), intent(in) :: lines(:)
    character(len=:), allocatable :: s
    integer :: ios, nl, i, j, ntok, n, m_in
    integer, allocatable :: vals(:)
    real(dp) :: kmax(3), kmin(3)
    logical :: ok, debug_unused

    qe_route = .false.
    ef_given = .true.
    nl = size(lines)
    if (nl < 8) call fatal('the input ' // input_name // ' has ' // istr(nl) // &
                           ' of the 8 required lines (or start it with &lindhard for the namelist input)')
    table_file = lines(1)%s
    n = len(table_file)
    if (n >= 2) then
      if (index('''"', table_file(1:1)) > 0 .and. table_file(n:n) == table_file(1:1)) &
        table_file = table_file(2:n-1)
    end if
    read(lines(2)%s, *, iostat=ios) nk_limit
    if (ios /= 0 .or. nk_limit < 1) call fatal('line 2 (nk): a positive integer is needed, got "' // lines(2)%s // '"')
    read(lines(3)%s, *, iostat=ios) m_in
    if (ios /= 0 .or. m_in < 1) call fatal('line 3 (m): a positive integer is needed, got "' // lines(3)%s // '"')
    col_first = 1
    col_last = m_in
    read(lines(4)%s, *, iostat=ios) ef
    if (ios /= 0) call fatal('line 4 (E_F): a number is needed, got "' // lines(4)%s // '"')
    read(lines(5)%s, *, iostat=ios) temp
    if (ios /= 0 .or. temp <= 0.0_dp) call fatal('line 5 (T): a temperature > 0 K is needed, got "' &
                                                 // lines(5)%s // '"')
    read(lines(6)%s, *, iostat=ios) kmax
    if (ios /= 0) call fatal('line 6 (max_k1 max_k2 max_k3): three numbers are needed, got "' // lines(6)%s // '"')
    read(lines(7)%s, *, iostat=ios) kmin
    if (ios /= 0) call fatal('line 7 (min_k1 min_k2 min_k3): three numbers are needed, got "' // lines(7)%s // '"')
    per = kmax - kmin
    s = lines(8)%s
    if (index(s, ' ') > 0) s = s(1:index(s, ' ') - 1)
    call parse_yes_no(s, debug_unused, ok)                 ! line 8 (debug) is read but not used
    if (.not. ok) call fatal('line 8 (debug): yes or no, got "' // lines(8)%s // '"')

    do i = 9, nl
      ntok = count_tokens(lines(i)%s)
      if (allocated(vals)) deallocate(vals)
      allocate(vals(ntok))
      read(lines(i)%s, *, iostat=ios) vals
      if (ios /= 0) call fatal('line ' // istr(i) // ': integers expected (N1 N2 N3, or a k-point index), got "' &
                               // lines(i)%s // '"')
      if (ntok == 3) then
        if (mesh_given) call fatal('two mesh lines (N1 N2 N3) in the input')
        if (any(vals < 1)) call fatal('mesh dimensions N1 N2 N3 must be >= 1')
        nmesh = vals
        mesh_given = .true.
      else
        do j = 1, ntok
          if (vals(j) >= 1) sel_index = [sel_index, vals(j)]
        end do
      end if
    end do
  end subroutine parse_legacy

  !-----------------------------------------------------------------------------
  subroutine prepare_from_qe()
    type(pw_data_t) :: pw
    character(len=:), allocatable :: err, fmt
    real(dp), allocatable :: ebands(:,:), xg(:), eg(:,:), lo(:), hi(:)
    integer, allocatable :: bsel(:), g(:,:), src_of_site(:)
    logical, allocatable :: from_nscf(:)
    integer :: nks, nb, nkg, nbg, i, b, s, src, nfill, nedge, nrow, nmiss, ndup, u, i1, i2, i3
    integer :: glo(3), ghi(3), v(3), gg(3), below, above
    real(dp) :: diff
    logical :: edge_axis(3)

    ! ---- pw.x output
    call read_pw_output(nscf_out, spin, pw, err)
    if (len(err) > 0) call fatal(err)
    nks = pw%nks
    write(output_unit, '(a)') '  input       ' // input_name
    write(output_unit, '(a)') '  nscf output ' // nscf_out // ': ' // istr(nks) // ' k-points, ' // &
      istr(pw%nbnd) // ' bands' // pick(pw%lsda, ', spin-polarised (LSDA)', '')
    if (pw%lsda .and. spin == 0) call fatal(nscf_out // ' is spin-polarised: set spin = 1 (up) or 2 (down). ' // &
      'chi0 of the two spin channels add up; run once for each.')
    if (.not. pw%lsda .and. spin /= 0) write(output_unit, '(a)') '  NOTE: spin is ignored (not an LSDA run)'
    if (.not. (pw%have_cryst .or. pw%have_cart)) call fatal(nscf_out // ' has no k-point list: run the nscf ' // &
      'with verbosity = ''high''')
    if (.not. pw%have_cryst) then
      if (.not. pw%have_at) call fatal(nscf_out // ': no crystal k-point list and no crystal axes')
      write(output_unit, '(a)') '              k in crystal coordinates converted from the cartesian list'
    end if
    if (pw%have_cart .and. pw%header_mismatch > 6.0e-4_dp) call fatal(nscf_out // ': the k of the band-energy ' // &
      'blocks does not follow the k-point list')

    ! ---- band energies
    if (len(bands_file) > 0) then
      call read_bands_gnu(bands_file, xg, eg, nkg, nbg, err)
      if (len(err) > 0) call fatal(err)
      if (nkg /= nks) call fatal(bands_file // ' has ' // istr(nkg) // ' k-points per band, ' // nscf_out // &
        ' lists ' // istr(nks) // '. Run bands.x on this nscf' // pick(pw%lsda, ' (spin_component = spin)', '') // &
        ' with max_lines in PP/src/bands.f90 raised above the number of k-points.')
      if (pw%have_e) then
        if (nbg /= pw%nbnd) call fatal(bands_file // ' has ' // istr(nbg) // ' bands, ' // nscf_out // &
                                       ' has ' // istr(pw%nbnd))
        diff = 0.0_dp
        do i = 1, nks
          diff = max(diff, maxval(abs(sort_small(eg(:, i)) - sort_small(pw%e(:, i)))))
        end do
        if (diff > 2.0e-4_dp) call fatal('the energies in ' // bands_file // ' differ from those in ' // nscf_out // &
          ' by up to ' // fstr(diff, '(f10.4)') // ' eV: not the same calculation or k-point order' // &
          pick(pw%lsda, ', or another spin component', ''))
        write(output_unit, '(a)') '  bands file  ' // bands_file // ': ' // istr(nbg) // ' bands; matches the ' // &
          'energies in ' // nscf_out // ' (largest difference ' // fstr(diff, '(f7.4)') // ' eV)'
      else
        write(output_unit, '(a)') '  bands file  ' // bands_file // ': ' // istr(nbg) // ' bands (' // nscf_out // &
          ' has no band energies to compare with)'
      end if
      call move_alloc(eg, ebands)
      nb = nbg
    else
      if (.not. pw%have_e) call fatal('no band energies in ' // nscf_out // ' (run the nscf with verbosity = ' // &
                                      '''high''), and no bands_file given')
      ebands = pw%e(:, 1:nks)
      nb = pw%nbnd
      write(output_unit, '(a)') '  energies    from ' // nscf_out
    end if

    ! ---- Fermi energy
    if (.not. ef_given) then
      if (pw%lsda .and. pw%have_ef_spin) then
        ef = pw%ef_spin(spin)
      else if (pw%have_ef) then
        ef = pw%ef
      else
        call fatal('no Fermi energy in ' // nscf_out // ' (fixed occupations?): set ef')
      end if
      write(output_unit, '(a)') '  E_F         ' // fstr(ef, '(f12.4)') // ' eV, from ' // nscf_out
    end if

    ! ---- band selection
    allocate(lo(nb), hi(nb))
    do b = 1, nb
      lo(b) = minval(ebands(b, :))
      hi(b) = maxval(ebands(b, :))
    end do
    if (band_first > 0 .or. band_last > 0) then
      if (band_last == 0) band_last = nb
      band_first = max(1, band_first)
      if (band_last > nb) call fatal('band_last = ' // istr(band_last) // ', but there are ' // istr(nb) // ' bands')
      bsel = [(b, b = band_first, band_last)]
    else
      bsel = pack([(b, b = 1, nb)], lo < ef .and. hi > ef)
      if (size(bsel) == 0) call fatal('no band crosses E_F = ' // fstr(ef, '(f10.4)') // ' eV: set band_first ' // &
                                      'and band_last')
    end if
    below = 0
    above = 0
    do b = 1, nb
      if (any(bsel == b)) cycle
      if (hi(b) <= ef) below = b
      if (lo(b) >= ef .and. above == 0) above = b
    end do
    write(output_unit, '(a)') '  bands used  ' // band_list(bsel) // pick(band_first > 0, ' (from the input)', &
                              ' (the bands that cross E_F)')
    if (below > 0) write(output_unit, '(a)') '              nearest band below: ' // istr(below) // ', top at ' // &
                                             fstr(hi(below) - ef, '(f8.3)') // ' eV from E_F'
    if (above > 0) write(output_unit, '(a)') '              nearest band above: ' // istr(above) // ', bottom at +' // &
                                             fstr(lo(above) - ef, '(f8.3)') // ' eV from E_F'
    b = count(lo < ef .and. hi > ef)
    if (b > count(lo(bsel) < ef .and. hi(bsel) > ef)) write(output_unit, '(a)') '  NOTE: bands ' // &
      band_list(pack([(i, i = 1, nb)], lo < ef .and. hi > ef)) // ' cross E_F; not all of them are used'

    ! ---- sanity check files 1 and 2: k.dat and selected.dat
    call write_k_dat(pw%kcry)
    call write_selected_dat(ebands, bsel, xg)

    ! ---- mesh sites of the nscf k-points
    call detect_mesh(pw%kcry, per)
    dk = per / nmesh
    nt = product(nmesh)
    allocate(g(3, nks), src_of_site(nt))
    src_of_site = 0
    ndup = 0
    do i = 1, nks
      g(:, i) = grid_vector(pw%kcry(:, i), 'k-point ' // istr(i) // ' of ' // nscf_out)
      s = site_index(modulo(g(:, i), nmesh), nmesh)
      if (src_of_site(s) == 0) then
        src_of_site(s) = i
      else
        ndup = ndup + 1
      end if
    end do
    glo = minval(g, dim=2)

    from_nscf = src_of_site > 0
    nfill = 0
    nmiss = 0
    do s = 1, nt
      if (from_nscf(s)) cycle
      v = site_vector(s, nmesh)
      src = site_index(modulo(-v, nmesh), nmesh)
      if (from_nscf(src)) then
        src_of_site(s) = src_of_site(src)
        nfill = nfill + 1
      else
        nmiss = nmiss + 1
      end if
    end do
    if (nmiss > 0) call fatal('the k-points of ' // nscf_out // ' cover ' // istr(nt - nmiss - nfill) // ' of the ' // &
      istr(nt) // ' sites of the ' // istr(nmesh(1)) // ' x ' // istr(nmesh(2)) // ' x ' // istr(nmesh(3)) // &
      ' mesh, and ' // istr(nmiss) // ' sites are missing also as -k. Run the nscf on the full Gamma-centred ' // &
      'mesh (K_POINTS automatic N1 N2 N3 0 0 0, nosym = .true.).')

    ! ---- sanity file 3: step1.dat, the table the calculation reads. The closed box from -1/2 to
    !      +1/2 along every axis with an even number of points (the +1/2 planes repeat the
    !      energies of the -1/2 planes, edges and corners included), ordered like the step2.dat of
    !      the earlier workflow: k1 runs fastest, then k2, then k3, each in ascending order.
    edge_axis = nmesh > 1 .and. mod(nmesh, 2) == 0 .and. glo == -nmesh / 2
    ghi = glo + nmesh - 1
    where (edge_axis) ghi = ghi + 1
    nrow = product(ghi - glo + 1)
    fmt = '(3f12.7,' // istr(size(bsel)) // 'f10.4)'
    open(newunit=u, file='step1.dat', status='replace', action='write')
    do i3 = glo(3), ghi(3)
      do i2 = glo(2), ghi(2)
        do i1 = glo(1), ghi(1)
          gg = [i1, i2, i3]
          write(u, fmt) gg * dk, ebands(bsel, src_of_site(site_index(modulo(gg, nmesh), nmesh)))
        end do
      end do
    end do
    close(u)
    nedge = nrow - nt
    write(output_unit, '(a)') '  mesh        ' // istr(nmesh(1)) // ' x ' // istr(nmesh(2)) // ' x ' // &
      istr(nmesh(3)) // ' = ' // istr(nt) // ' sites' // pick(mesh_given, ' (from the input)', ' (detected)')
    if (ndup > 0) write(output_unit, '(a)') '  NOTE: ' // istr(ndup) // ' k-points of ' // nscf_out // &
                                            ' repeat a mesh site'
    if (nfill > 0) write(output_unit, '(a)') '              ' // istr(nfill) // ' sites restored from -k ' // &
      '(time reversal: e(-k) = e(k), as used by pw.x to reduce the mesh)'
    write(output_unit, '(a)') '  written     step1.dat: ' // istr(nrow) // ' k-points = ' // &
      istr(ghi(1) - glo(1) + 1) // ' x ' // istr(ghi(2) - glo(2) + 1) // ' x ' // istr(ghi(3) - glo(3) + 1) // &
      pick(nedge > 0, ' (' // istr(nt) // ' mesh points + ' // istr(nedge) // ' on the +1/2 zone edges)', '') // &
      ', k1 fastest'
    write(output_unit, '(a)') '              ' // box_text(glo, ghi)

    table_file = 'step1.dat'
    nk_limit = -1
    col_first = 1
    col_last = size(bsel)
    band_label = bsel
    source_note = nscf_out // pick(len(bands_file) > 0, ' (k-points) + ' // bands_file // ' (energies)', '') // &
                  ', bands ' // band_list(bsel) // ' -> step1.dat'
  end subroutine prepare_from_qe


  function band_list(b) result(s)
    integer, intent(in) :: b(:)
    character(len=:), allocatable :: s
    integer :: i
    if (size(b) > 1) then
      if (all(b(2:) - b(:size(b)-1) == 1)) then
        s = istr(b(1)) // '-' // istr(b(size(b)))
        return
      end if
    end if
    s = istr(b(1))
    do i = 2, size(b)
      s = s // ',' // istr(b(i))
    end do
  end function band_list

  !> 'k1 -0.5000 .. 0.5000 (225 values), ...' for the integer coordinate ranges glo..ghi
  function box_text(glo, ghi) result(s)
    integer, intent(in) :: glo(3), ghi(3)
    character(len=:), allocatable :: s
    integer :: a
    s = ''
    do a = 1, 3
      s = s // pick(a > 1, ',  ', '') // 'k' // istr(a) // ' ' // fstr(glo(a) * dk(a), '(f8.4)') // ' .. ' // &
          fstr(ghi(a) * dk(a), '(f8.4)') // ' (' // istr(ghi(a) - glo(a) + 1) // ')'
    end do
  end function box_text

  !-----------------------------------------------------------------------------
  subroutine write_k_dat(kcry)
    real(dp), intent(in) :: kcry(:,:)
    integer :: u, i
    open(newunit=u, file='k.dat', status='replace', action='write')
    do i = 1, size(kcry, 2)
      write(u, '(3f12.7)') kcry(:, i)
    end do
    close(u)
    write(output_unit, '(a)') '  written     k.dat: ' // istr(size(kcry, 2)) // ' k-points (crystal ' // &
                              'coordinates, order of ' // nscf_out // ')'
  end subroutine write_k_dat

  !> one block per selected band: x and energy (bands.x gnu), or k-point number and energy
  subroutine write_selected_dat(ebands, bsel, xg)
    real(dp), intent(in) :: ebands(:,:)
    integer, intent(in) :: bsel(:)
    real(dp), allocatable, intent(in) :: xg(:)
    integer :: u, i, j
    open(newunit=u, file='selected.dat', status='replace', action='write')
    do j = 1, size(bsel)
      do i = 1, size(ebands, 2)
        if (allocated(xg)) then
          write(u, '(2f10.4)') xg(i), ebands(bsel(j), i)
        else
          write(u, '(i10,f10.4)') i, ebands(bsel(j), i)
        end if
      end do
      write(u, '(a)') ''
    end do
    close(u)
    write(output_unit, '(a)') '  written     selected.dat: bands ' // band_list(bsel) // ', one block per band (' // &
      pick(allocated(xg), 'x and energy as in ' // bands_file, 'k-point number and energy') // ')'
  end subroutine write_selected_dat

  !-----------------------------------------------------------------------------
  !> the k-point table (step1.dat, or band_file): all lines unless nk_limit > 0,
  !> energy columns col_first..col_last (col_last = 0: up to the last column)
  subroutine read_table()
    character(len=:), allocatable :: s
    real(dp), allocatable :: tmp(:)
    integer :: u, ios, i, nline, ncol, b

    open(newunit=u, file=table_file, status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('cannot open the band file "' // table_file // '"')
    nk = 0
    ncol = -1
    do                                                   ! pass 1: k-points and columns
      call read_line(u, s, ios)
      if (ios /= 0) exit
      s = trim(adjustl(s))
      if (len(s) == 0) cycle
      if (s(1:1) == '#') cycle
      if (ncol < 0) ncol = count_tokens(s) - 3
      nk = nk + 1
      if (nk == nk_limit) exit
    end do
    if (nk_limit > 0 .and. nk < nk_limit) call fatal(table_file // ' has ' // istr(nk) // ' k-points, but nk = ' // &
                                                     istr(nk_limit))
    if (nk == 0 .or. ncol < 1) call fatal(table_file // ': no lines of the form k1 k2 k3 e_1 ...')
    if (col_last == 0) col_last = ncol
    if (col_last > ncol) call fatal(table_file // ' has ' // istr(ncol) // ' energy columns, ' // istr(col_last) // &
                                    ' requested')
    if (col_first > col_last) call fatal('band_first > number of energy columns')
    m = col_last - col_first + 1
    if (.not. allocated(band_label)) band_label = [(b, b = col_first, col_last)]

    rewind(u)
    allocate(kpt(3, nk), eig(m, nk), tmp(col_last))
    i = 0
    nline = 0
    do while (i < nk)
      call read_line(u, s, ios)
      nline = nline + 1
      s = trim(adjustl(s))
      if (len(s) == 0) cycle
      if (s(1:1) == '#') cycle
      i = i + 1
      read(s, *, iostat=ios) kpt(:, i), tmp
      if (ios /= 0) call fatal(table_file // ', line ' // istr(nline) // ': k1 k2 k3 and at least ' // &
                               istr(col_last) // ' energies expected')
      eig(:, i) = tmp(col_first:col_last)
    end do
    close(u)
  end subroutine read_table

  !-----------------------------------------------------------------------------
  !> Mesh dimensions (unless given): along each axis the smallest N for which every
  !> k_a is a multiple of period/N (to 0.1 of a step). This holds for k printed with
  !> few decimals and for meshes that pw.x reduced by time reversal.
  subroutine detect_mesh(k, period)
    real(dp), intent(in) :: k(:,:)
    real(dp), intent(inout) :: period(3)
    real(dp), allocatable :: v(:), x(:)
    real(dp) :: span, gap
    integer :: a, i, nd, n, nmax
    logical :: flat

    do a = 1, 3
      span = maxval(k(a, :)) - minval(k(a, :))
      flat = span < 1.0e-6_dp
      if (period(a) <= 0.0_dp) then
        if (.not. flat) call fatal('the period (max_k - min_k) must be > 0 along axis ' // istr(a))
        period(a) = 1.0_dp                               ! single plane: any period works
      end if
      if (mesh_given) cycle
      if (flat) then
        nmesh(a) = 1
        cycle
      end if
      v = sorted_copy(k(a, :))                           ! distinct values and the smallest gap
      nd = 1
      gap = huge(1.0_dp)
      do i = 2, size(v)
        if (v(i) - v(nd) > 1.0e-6_dp) then
          gap = min(gap, v(i) - v(nd))
          nd = nd + 1
          v(nd) = v(i)
        end if
      end do
      nmax = nint(1.5_dp * period(a) / gap) + 2
      nmesh(a) = nint(period(a) / gap)
      do n = 1, nmax
        x = v(1:nd) * n / period(a)
        if (all(abs(x - nint(x)) < 0.1_dp)) then
          nmesh(a) = n
          exit
        end if
      end do
    end do
  end subroutine detect_mesh

  !> integer coordinates g = k/dk counted from Gamma (stops if k is not a mesh point)
  function grid_vector(kv, label) result(g)
    real(dp), intent(in) :: kv(3)
    character(len=*), intent(in) :: label
    integer :: g(3), a
    real(dp) :: x
    do a = 1, 3
      x = kv(a) / dk(a)
      g(a) = nint(x)
      if (abs(x - g(a)) > 0.25_dp) call fatal(label // ' (' // fstr(kv(1), '(f12.6)') // ' ' // &
        fstr(kv(2), '(f12.6)') // ' ' // fstr(kv(3), '(f12.6)') // ') is not a point of the Gamma-centred ' // &
        istr(nmesh(1)) // ' x ' // istr(nmesh(2)) // ' x ' // istr(nmesh(3)) // ' mesh. Check the mesh, the ' // &
        'period and that k is in crystal coordinates; shifted Monkhorst-Pack meshes are not supported.')
    end do
  end function grid_vector

  !-----------------------------------------------------------------------------
  !> Mesh site of every k-point of the table. Along each axis the site index is
  !> g mod N, so k+q is site (g_k + g_q) mod N for any layout of the table.
  subroutine map_to_mesh()
    real(dp) :: mismatch
    integer :: i, s, ndup

    call detect_mesh(kpt, per)
    dk = per / nmesh
    nt = product(nmesh)
    allocate(site_of(nk), point_of_site(nt), rows_of_site(nt))
    point_of_site = 0
    rows_of_site = 0
    ndup = 0
    mismatch = 0.0_dp
    do i = 1, nk
      s = site_index(modulo(grid_vector(kpt(:, i), 'k-point ' // istr(i) // ' of ' // table_file), nmesh), nmesh)
      site_of(i) = s
      rows_of_site(s) = rows_of_site(s) + 1
      if (point_of_site(s) == 0) then
        point_of_site(s) = i
      else                                             ! periodic copy (e.g. the +1/2 plane) or repeat
        ndup = ndup + 1
        mismatch = max(mismatch, maxval(abs(eig(:, i) - eig(:, point_of_site(s)))))
      end if
    end do
    i = count(point_of_site == 0)
    if (i > 0) call fatal('the mesh is incomplete: ' // istr(i) // ' of the ' // istr(nt) // ' sites of the ' // &
      istr(nmesh(1)) // ' x ' // istr(nmesh(2)) // ' x ' // istr(nmesh(3)) // ' mesh have no k-point. ' // &
      'Run pw.x on the full mesh (nosym = .true., noinv = .true.) or unfold the irreducible wedge first' // &
      pick(mesh_given, ', or check N1 N2 N3', '') // '.')

    if (.not. qe_route) write(output_unit, '(a)') '  input       ' // input_name
    write(output_unit, '(a)') '  table       ' // table_file // ': ' // istr(nk) // ' k-points, ' // &
                              istr(m) // ' bands used'
    if (.not. qe_route) write(output_unit, '(a)') '  mesh        ' // istr(nmesh(1)) // ' x ' // istr(nmesh(2)) // &
      ' x ' // istr(nmesh(3)) // ' = ' // istr(nt) // ' distinct k-points' // pick(mesh_given, ' (from the input)', &
      ' (detected)')
    if (ndup > 0) then
      write(output_unit, '(a)') '  copies      ' // istr(ndup) // ' k-points are periodic copies of others (largest ' // &
                                'energy difference ' // fstr(mismatch, '(es9.2)') // ' eV)'
      if (mismatch > 1.0e-3_dp) write(output_unit, '(a)') '  WARNING: periodic copies have different ' // &
        'energies. Is max_k - min_k one reciprocal-lattice period (1 in crystal coordinates)?'
    end if
  end subroutine map_to_mesh

  !-----------------------------------------------------------------------------
  subroutine report_setup()
    integer :: b, j, i
    real(dp) :: lo, hi
    logical :: crosses

    kt = KB * temp
    write(output_unit, '(a)') '  E_F         ' // fstr(ef, '(f12.6)') // ' eV'
    write(output_unit, '(a)') '  real part   T = ' // fstr(temp, '(f12.3)') // ' K (kB*T = ' // &
      fstr(kt, '(f12.7)') // ' eV), Fermi-Dirac occupations'
    write(output_unit, '(a)') '  imag. part  Im chi0(q,w)/w at w -> 0, Dirac delta = Gaussian of width ' // &
      'sigma = ' // fstr(sigma, '(f10.4)') // ' eV'
    write(output_unit, '(a)') '  weights     each of the ' // istr(nk) // ' k-points of ' // table_file // &
      ' has weight 1, N = ' // istr(nk)
    crosses = .false.
    do b = 1, m
      lo = minval(eig(b, :)) - ef
      hi = maxval(eig(b, :)) - ef
      write(output_unit, '(a,i4,a,f9.3,a,f9.3,a,a)') '  band', band_label(b), ':', lo, ' ..', hi, &
        ' eV around E_F', pick(lo < 0.0_dp .and. hi > 0.0_dp, '  (crosses E_F)', '')
      crosses = crosses .or. (lo < 0.0_dp .and. hi > 0.0_dp)
    end do
    if (.not. crosses) write(output_unit, '(a)') '  NOTE: no band crosses E_F (check E_F and the bands used)'
    call report_energy_steps()

    ! selected k-points: one per mesh site
    nsel = 0
    allocate(sel_site(size(sel_index)))
    do j = 1, size(sel_index)
      i = sel_index(j)
      if (i > nk) call fatal('selected k-point ' // istr(i) // ' is beyond the ' // istr(nk) // ' k-points of ' // &
                             table_file)
      if (any(sel_site(1:nsel) == site_of(i))) then
        write(output_unit, '(a)') '  k-point ' // istr(i) // ' is the same mesh site as an earlier selected ' // &
                                  'k-point: skipped'
        cycle
      end if
      nsel = nsel + 1
      sel_site(nsel) = site_of(i)
      sel_index(nsel) = i
      write(output_unit, '(a,3f10.5,a)') '  selected    k-point ' // istr(i) // ': k = (', kpt(:, i), ')'
    end do
    sel_index = sel_index(1:nsel)
    sel_site = sel_site(1:nsel)
  end subroutine report_setup

  !-----------------------------------------------------------------------------
  subroutine report_energy_steps()
    real(dp), allocatable :: d(:), srt(:)
    character(len=:), allocatable :: s
    integer :: a, site, b, nd, v(3), nb
    real(dp) :: e1

    s = ''
    do a = 1, 3
      if (nmesh(a) < 2) cycle
      allocate(d(nt * m))
      nd = 0
      do site = 1, nt
        v = site_vector(site, nmesh)
        v(a) = modulo(v(a) + 1, nmesh(a))
        nb = point_of_site(site_index(v, nmesh))
        do b = 1, m
          e1 = eig(b, point_of_site(site))
          if (abs(e1 - ef) > 0.1_dp) cycle
          nd = nd + 1
          d(nd) = abs(eig(b, nb) - e1)
        end do
      end do
      if (nd > 0) then
        srt = sorted_copy(d(1:nd))
        s = s // pick(len(s) > 0, ',  ', '') // 'k' // istr(a) // ' ' // fstr(1000.0_dp * srt((nd + 1) / 2), '(f8.1)') // &
            ' meV'
      end if
      deallocate(d)
    end do
    if (len(s) > 0) then
      write(output_unit, '(a)') '  energy step between neighbouring k-points near E_F (median): ' // s
      write(output_unit, '(a)') '              kB*T and sigma much smaller than these steps give sampling noise'
    end if
  end subroutine report_energy_steps


  !-----------------------------------------------------------------------------
  subroutine share_setup()
    integer :: s, j

    call MPI_Bcast(nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nmesh, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nsel, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ef, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(temp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(sigma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    allocate(e(m, nt), wsite(nt))
    if (is_root) then
      do s = 1, nt
        e(:, s) = eig(:, point_of_site(s))
      end do
      wsite = real(rows_of_site, dp)                     ! every line of the table has weight 1
    else
      allocate(sel_site(nsel))
    end if
    call MPI_Bcast(e, m * nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(wsite, nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(sel_site, nsel, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    wsum = sum(wsite)

    kt = KB * temp
    allocate(f(m, nt), w(m, nt), gdel(nt))
    f = 0.5_dp * (1.0_dp - tanh(0.5_dp * (e - ef) / kt))   ! Fermi-Dirac without overflow
    w = f * (1.0_dp - f) / kt                              ! -df/de
    gdel = sum(exp(-0.5_dp * ((e - ef) / sigma)**2), dim=1) / (sqrt(2.0_dp * PI) * sigma)   ! sum over bands
    allocate(sel_of_site(nt))
    sel_of_site = 0
    do j = 1, nsel
      sel_of_site(sel_site(j)) = j
    end do
  end subroutine share_setup

  !-----------------------------------------------------------------------------
  !> sum over band pairs of the Lindhard term for k (site ik) and k+q (site jk)
  pure real(dp) function pair_sum(ik, jk) result(acc)
    integer, intent(in) :: ik, jk
    integer :: a, b
    real(dp) :: e1, f1, w1, de

    acc = 0.0_dp
    do a = 1, m
      e1 = e(a, ik)
      f1 = f(a, ik)
      w1 = w(a, ik)
      do b = 1, m
        de = e1 - e(b, jk)
        if (abs(de) > DEG) then
          acc = acc - (f1 - f(b, jk)) / de
        else
          acc = acc + 0.5_dp * (w1 + w(b, jk))           ! limit -df/de of the degenerate pair
        end if
      end do
    end do
  end function pair_sum

  !-----------------------------------------------------------------------------
  !>   Re chi0(q)          = -(1/N) sum_k w_k sum_nm [f(e_nk) - f(e_m,k+q)] / (e_nk - e_m,k+q)
  !>   Im chi0(q,w)/w|w->0 =  (pi/N) sum_k w_k sum_nm delta(e_nk - E_F) delta(e_m,k+q - E_F)
  !> (the second is pi times the Fermi-surface nesting function)
  subroutine compute_chi()
    real(dp), allocatable :: re_loc(:), im_loc(:), cre_loc(:,:), cim_loc(:,:)
    real(dp) :: inv_n, acc_re, acc_im, term, tim, t0, now
    integer :: iq, sq(3), i1, i2, i3, j1, j2, j3, ik, jk, ik0, jk0, nmine, ndone, next_pct, js

    allocate(re_loc(nt), im_loc(nt), cre_loc(nt, nsel), cim_loc(nt, nsel))
    re_loc = 0.0_dp
    im_loc = 0.0_dp
    cre_loc = 0.0_dp
    cim_loc = 0.0_dp
    inv_n = 1.0_dp / wsum
    nmine = 0
    if (rank < nt) nmine = (nt - 1 - rank) / nproc + 1
    if (is_root) then
      write(output_unit, '(a)') '  computing Re and Im chi0 at ' // istr(nt) // ' q-points on ' // istr(nproc) // &
        ' ranks (' // istr(nmine) // ' q on rank 0)'
      flush(output_unit)
    end if

    t0 = MPI_Wtime()
    ndone = 0
    next_pct = 10
    do iq = rank + 1, nt, nproc
      sq = site_vector(iq, nmesh)                        ! k+q = site (g_k + g_q) mod N
      acc_re = 0.0_dp
      acc_im = 0.0_dp
      do i3 = 0, nmesh(3) - 1
        j3 = i3 + sq(3)
        if (j3 >= nmesh(3)) j3 = j3 - nmesh(3)
        do i2 = 0, nmesh(2) - 1
          j2 = i2 + sq(2)
          if (j2 >= nmesh(2)) j2 = j2 - nmesh(2)
          ik0 = 1 + nmesh(1) * (i2 + nmesh(2) * i3)
          jk0 = 1 + nmesh(1) * (j2 + nmesh(2) * j3)
          do i1 = 0, nmesh(1) - 1
            j1 = i1 + sq(1)
            if (j1 >= nmesh(1)) j1 = j1 - nmesh(1)
            ik = ik0 + i1
            jk = jk0 + j1
            term = pair_sum(ik, jk)                      ! real part
            acc_re = acc_re + wsite(ik) * term
            tim = gdel(ik) * gdel(jk)                    ! imaginary part
            acc_im = acc_im + wsite(ik) * tim
            js = sel_of_site(ik)
            if (js > 0) then                             ! one k-point
              cre_loc(iq, js) = term * inv_n
              cim_loc(iq, js) = PI * tim * inv_n
            end if
          end do
        end do
      end do
      re_loc(iq) = acc_re * inv_n
      im_loc(iq) = PI * acc_im * inv_n

      ndone = ndone + 1
      if (is_root .and. 100 * ndone >= next_pct * nmine) then
        now = MPI_Wtime() - t0
        write(output_unit, '(a,i3,a,a,a,a,a)') '    ', next_pct, '%   elapsed ', fstr(now, '(f10.1)'), &
          ' s   remaining ', fstr(now * real(nmine - ndone, dp) / real(ndone, dp), '(f10.1)'), ' s'
        flush(output_unit)
        do while (100 * ndone >= next_pct * nmine)
          next_pct = next_pct + 10
        end do
      end if
    end do
    t_compute = MPI_Wtime() - t0

    allocate(chi_re(nt), chi_im(nt), con_re(nt, nsel), con_im(nt, nsel))
    call MPI_Reduce(re_loc, chi_re, nt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(im_loc, chi_im, nt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (nsel > 0) then
      call MPI_Reduce(cre_loc, con_re, nt * nsel, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(cim_loc, con_im, nt * nsel, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if

    if (is_root) write(output_unit, '(a)') '  chi0 done in ' // fstr(t_compute, '(f12.2)') // ' s'
  end subroutine compute_chi

  !-----------------------------------------------------------------------------
  !> chi_<E_F>_<T>.dat: Re and Im chi0 at every q-point (the lines of the table,
  !> in their order); chi_k<index>_<E_F>_<T>.dat: the contribution of each
  !> selected k-point
  subroutine write_results()
    character(len=:), allocatable :: fname, tag
    integer :: u, i, j, s, ik, gk(3), gq(3), glow(3), gkq(3)
    real(dp) :: kq(3)

    tag = '_' // name_number(ef) // '_' // name_number(temp) // '.dat'
    fname = 'chi' // tag
    open(newunit=u, file=fname, status='replace', action='write')
    write(u, '(a)') '# chi0(q): real and imaginary part at every q-point'
    call write_header(u)
    write(u, '(a)') '# columns: q1 q2 q3 (crystal coordinates, the k-points of ' // table_file // ')  Re chi0  Im chi0/w'
    do i = 1, nk
      write(u, '(3(f11.6,1x),2es19.10e3)') kpt(:, i), chi_re(site_of(i)), chi_im(site_of(i))
    end do
    close(u)
    write(output_unit, '(a)') '  written     ' // fname // ' (q1 q2 q3  Re chi0  Im chi0/w, ' // istr(nk) // ' q-points)'

    ! k+q is written inside the k box [min_k, min_k + period)
    glow = ceiling(minval(kpt, dim=2) / dk - 1.0e-6_dp)
    do j = 1, nsel
      ik = sel_site(j)
      gk = site_vector(ik, nmesh)
      fname = 'chi_k' // istr(sel_index(j)) // tag
      open(newunit=u, file=fname, status='replace', action='write')
      write(u, '(a,3(1x,f10.6),a)') '# contribution of k-point ' // istr(sel_index(j)) // ', k =', &
        kpt(:, sel_index(j)), ', to chi0(q); the sum over all k-points of the table gives chi' // tag
      call write_header(u)
      write(u, '(a)') '# columns: q1 q2 q3  (k+q)1 (k+q)2 (k+q)3  Re chi0_k(q)  Im chi0_k(q)/w'
      do i = 1, nk                                       ! one line per q of the table, same order
        s = site_of(i)
        gq = site_vector(s, nmesh)
        gkq = glow + modulo(gk + gq - glow, nmesh)
        kq = gkq * dk
        write(u, '(6(f11.6,1x),2es19.10e3)') kpt(:, i), kq, con_re(s, j), con_im(s, j)
      end do
      close(u)
      write(output_unit, '(a)') '  written     ' // fname
    end do
  end subroutine write_results

  !-----------------------------------------------------------------------------
  !> a number for a file name: up to 4 decimals, trailing zeros dropped
  !> (6.4521 -> 6.4521, 300.0 -> 300, 0.50 -> 0.5)
  function name_number(x) result(s)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=40) :: buf
    integer :: n
    write(buf, '(f40.4)') x
    s = trim(adjustl(buf))
    n = len(s)
    do while (n > 1 .and. s(n:n) == '0')
      n = n - 1
    end do
    if (s(n:n) == '.') n = n - 1
    s = s(1:n)
    if (s == '-0') s = '0'
  end function name_number

  !-----------------------------------------------------------------------------
  subroutine write_header(unit)
    integer, intent(in) :: unit
    write(unit, '(a)') '# Re chi0(q, w=0)     = -(1/N) sum_k sum_nm [f(e_nk) - f(e_m,k+q)] / (e_nk - e_m,k+q)' // &
      '   [states/eV], f = Fermi-Dirac at T'
    write(unit, '(a)') '# Im chi0(q,w)/w|w->0 =  (pi/N) sum_k sum_nm delta(e_nk - E_F) delta(e_m,k+q - E_F)' // &
      '   [1/eV^2], delta = Gaussian of width sigma; nesting function xi(q) = Im chi0/w / pi'
    write(unit, '(a)') '# band file: ' // table_file // pick(len(source_note) > 0, ' (' // source_note // ')', '')
    write(unit, '(a)') '# E_F = ' // fstr(ef, '(f12.6)') // ' eV, T = ' // fstr(temp, '(f12.3)') // ' K (kB*T = ' // &
      fstr(kt, '(f12.7)') // ' eV), sigma = ' // fstr(sigma, '(f10.4)') // ' eV, bands ' // band_list(band_label)
    write(unit, '(a,3(1x,f9.6))') '# mesh ' // istr(nmesh(1)) // ' ' // istr(nmesh(2)) // ' ' // &
      istr(nmesh(3)) // ', period', per
    write(unit, '(a)') '# k-point weights: every line of ' // table_file // ' is one k-point of weight 1, N = ' // &
      istr(nk)
    write(unit, '(a)') '# chi0 > 0 convention, constant matrix elements, per cell and per spin'
    write(unit, '(a)') '# lindhard_chi.f90 (Fortran/MPI, ' // istr(nproc) // ' ranks)'
  end subroutine write_header

end program lindhard_chi
