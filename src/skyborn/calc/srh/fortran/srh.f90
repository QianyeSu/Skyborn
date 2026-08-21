module srh_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
  use iso_c_binding, only: c_double
  implicit none
  private
  public :: srh_profile_impl, srh_grid_impl

contains

  !--------------------------------------------------------------------------
  ! Linear interpolation of val(z) on the sorted, ascending abscissa h(:).
  ! Points outside [h(1), h(m)] return NaN, matching metpy's interpolate_1d
  ! with its default fill_value=np.nan.
  real(c_double) function interp_linear(z, m, h, val) result(v)
    real(c_double), intent(in) :: z
    integer, intent(in) :: m
    real(c_double), intent(in) :: h(m), val(m)
    real(c_double) :: w
    integer :: i

    v = ieee_value(0.0d0, ieee_quiet_nan)
    if (m < 2) return
    if (z < h(1) .or. z > h(m)) return
    if (z == h(1)) then
      v = val(1)
      return
    end if
    do i = 1, m - 1
      if (h(i) <= z .and. z <= h(i + 1)) then
        if (abs(h(i + 1) - h(i)) < 1.0d-12) then
          v = val(i)
        else
          w = (z - h(i)) / (h(i + 1) - h(i))
          v = val(i) + w * (val(i + 1) - val(i))
        end if
        return
      end if
    end do
  end function interp_linear

  !--------------------------------------------------------------------------
  ! Insert (z_bound, u_bound, v_bound) into the ordered layer arrays,
  ! keeping ascending order by height. Caller checks membership first.
  subroutine insert_bound(n, h_lay, u_lay, v_lay, z_bound, u_bound, v_bound)
    integer, intent(inout) :: n
    real(c_double), intent(inout) :: h_lay(:), u_lay(:), v_lay(:)
    real(c_double), intent(in) :: z_bound, u_bound, v_bound
    integer :: i, pos

    pos = n + 1
    do i = 1, n
      if (h_lay(i) > z_bound) then
        pos = i
        exit
      end if
    end do
    do i = n, pos, -1
      h_lay(i + 1) = h_lay(i)
      u_lay(i + 1) = u_lay(i)
      v_lay(i + 1) = v_lay(i)
    end do
    h_lay(pos) = z_bound
    u_lay(pos) = u_bound
    v_lay(pos) = v_bound
    n = n + 1
  end subroutine insert_bound

  !--------------------------------------------------------------------------
  ! Sort h ascending together with u and v (simple insertion sort).
  subroutine sort3(m, h, u, v)
    integer, intent(in) :: m
    real(c_double), intent(inout) :: h(m), u(m), v(m)
    real(c_double) :: th, tu, tv
    integer :: i, j
    do i = 2, m
      th = h(i)
      tu = u(i)
      tv = v(i)
      j = i - 1
      do while (j >= 1 .and. h(j) > th)
        h(j + 1) = h(j)
        u(j + 1) = u(j)
        v(j + 1) = v(j)
        j = j - 1
      end do
      h(j + 1) = th
      u(j + 1) = tu
      v(j + 1) = tv
    end do
  end subroutine sort3

  !--------------------------------------------------------------------------
  ! Single-profile storm-relative helicity.
  ! Mirrors metpy: AGL heights, layer [bottom, bottom+depth] with the
  ! boundary heights linearly interpolated (NaN when out of range), then
  ! int_layers = (ur_{n+1}*vr_n - ur_n*vr_{n+1}) summed with sign split.
  subroutine srh_profile_impl(nlev, height_m, u_ms, v_ms, depth_m, bottom_m, &
                              storm_u_ms, storm_v_ms, srh_pos, srh_neg, srh_total)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: height_m(nlev)
    real(c_double), intent(in) :: u_ms(nlev)
    real(c_double), intent(in) :: v_ms(nlev)
    real(c_double), intent(in) :: depth_m, bottom_m
    real(c_double), intent(in) :: storm_u_ms, storm_v_ms
    real(c_double), intent(out) :: srh_pos, srh_neg, srh_total

    real(c_double), allocatable :: h(:), u(:), v(:)
    real(c_double), allocatable :: h_lay(:), u_lay(:), v_lay(:)
    real(c_double) :: hmin, top_m, u_b, v_b, ur, vr, term
    integer :: m, k, n_lay, has_bottom, has_top

    srh_pos = 0.0d0
    srh_neg = 0.0d0
    srh_total = 0.0d0

    allocate(h(nlev), u(nlev), v(nlev))
    m = 0
    do k = 1, nlev
      if (ieee_is_finite(height_m(k)) .and. ieee_is_finite(u_ms(k)) .and. &
          ieee_is_finite(v_ms(k))) then
        m = m + 1
        h(m) = height_m(k)
        u(m) = u_ms(k)
        v(m) = v_ms(k)
      end if
    end do
    if (m < 2) then
      deallocate(h, u, v)
      return
    end if

    ! Above ground level
    hmin = minval(h(1:m))
    h(1:m) = h(1:m) - hmin

    ! Ascending sort
    call sort3(m, h, u, v)

    top_m = bottom_m + depth_m

    ! Collect data levels within [bottom, top] (ascending)
    allocate(h_lay(2 * m + 2), u_lay(2 * m + 2), v_lay(2 * m + 2))
    n_lay = 0
    do k = 1, m
      if (h(k) >= bottom_m .and. h(k) <= top_m) then
        n_lay = n_lay + 1
        h_lay(n_lay) = h(k)
      end if
    end do

    ! Insert boundary heights with interpolated winds if not exactly present
    has_bottom = 0
    has_top = 0
    do k = 1, n_lay
      if (h_lay(k) == bottom_m) has_bottom = 1
      if (h_lay(k) == top_m) has_top = 1
    end do
    if (has_bottom == 0) then
      u_b = interp_linear(bottom_m, m, h, u)
      call insert_bound(n_lay, h_lay, u_lay, v_lay, bottom_m, u_b, &
                        interp_linear(bottom_m, m, h, v))
    end if
    if (has_top == 0) then
      u_b = interp_linear(top_m, m, h, u)
      call insert_bound(n_lay, h_lay, u_lay, v_lay, top_m, u_b, &
                        interp_linear(top_m, m, h, v))
    end if

    ! Note: u_lay for data levels is filled lazily below via interpolation
    ! from the raw profile (identical to metpy interpolate_1d).
    do k = 1, n_lay
      u_lay(k) = interp_linear(h_lay(k), m, h, u) - storm_u_ms
      v_lay(k) = interp_linear(h_lay(k), m, h, v) - storm_v_ms
    end do

    ! Sum the cross terms, skipping NaN entries (metpy's masking semantics)
    do k = 1, n_lay - 1
      ur = u_lay(k + 1) * v_lay(k) - u_lay(k) * v_lay(k + 1)
      if (ieee_is_finite(ur)) then
        if (ur > 0.0d0) then
          srh_pos = srh_pos + ur
        else if (ur < 0.0d0) then
          srh_neg = srh_neg + ur
        end if
        srh_total = srh_total + ur
      end if
    end do

    deallocate(h, u, v, h_lay, u_lay, v_lay)
  end subroutine srh_profile_impl

  !--------------------------------------------------------------------------
  ! Grid version: depth, bottom, and storm motion are scalars.
  subroutine srh_grid_impl(nlev, nlat, nlon, h3, u3, v3, depth_m, bottom_m, &
                           storm_u_ms, storm_v_ms, pos2, neg2, tot2)
    integer, intent(in) :: nlev, nlat, nlon
    real(c_double), intent(in) :: h3(nlev, nlat, nlon)
    real(c_double), intent(in) :: u3(nlev, nlat, nlon)
    real(c_double), intent(in) :: v3(nlev, nlat, nlon)
    real(c_double), intent(in) :: depth_m, bottom_m
    real(c_double), intent(in) :: storm_u_ms, storm_v_ms
    real(c_double), intent(out) :: pos2(nlat, nlon)
    real(c_double), intent(out) :: neg2(nlat, nlon)
    real(c_double), intent(out) :: tot2(nlat, nlon)
    integer :: j, i
    do i = 1, nlon
      do j = 1, nlat
        call srh_profile_impl(nlev, h3(:, j, i), u3(:, j, i), v3(:, j, i), &
                              depth_m, bottom_m, storm_u_ms, storm_v_ms, &
                              pos2(j, i), neg2(j, i), tot2(j, i))
      end do
    end do
  end subroutine srh_grid_impl

end module srh_mod
