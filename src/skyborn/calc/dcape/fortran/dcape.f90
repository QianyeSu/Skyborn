module dcape_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
  use iso_c_binding, only: c_double
  implicit none
  private

  public :: dcape_profile_impl, dcape_grid_impl

  real(c_double), parameter :: eps = 0.622d0
  real(c_double), parameter :: rd = 287.05d0
  real(c_double), parameter :: cp = 1004.0d0
  real(c_double), parameter :: lv = 2.5d6
  real(c_double), parameter :: kappa = rd / cp

contains

  pure real(c_double) function sat_vapor_pressure_hpa(tc) result(es)
    real(c_double), intent(in) :: tc
    es = 6.112d0 * exp(17.67d0 * tc / (tc + 243.5d0))
  end function sat_vapor_pressure_hpa

  pure real(c_double) function sat_vapor_pressure_pa_from_tk(tk) result(es_pa)
    real(c_double), intent(in) :: tk
    es_pa = sat_vapor_pressure_hpa(tk - 273.15d0) * 100.0d0
  end function sat_vapor_pressure_pa_from_tk

  pure real(c_double) function mixing_ratio_from_e_pa(p_pa, e_pa) result(r)
    real(c_double), intent(in) :: p_pa, e_pa
    real(c_double) :: e
    e = max(e_pa, 1.0d-8)
    e = min(e, 0.99d0 * p_pa)
    r = eps * e / (p_pa - e)
  end function mixing_ratio_from_e_pa

  pure real(c_double) function virtual_temp_k(tk, mixing_ratio) result(tv)
    real(c_double), intent(in) :: tk, mixing_ratio
    tv = tk * (mixing_ratio + eps) / (eps * (1.0d0 + mixing_ratio))
  end function virtual_temp_k

  pure real(c_double) function thetae_bolton_like_metpy(p_hpa, tc, tdc) result(thetae)
    real(c_double), intent(in) :: p_hpa, tc, tdc
    real(c_double) :: t, td, p_pa, e_pa, r, t_l, th_l

    t = tc + 273.15d0
    td = max(tdc + 273.15d0, 173.15d0)
    p_pa = p_hpa * 100.0d0
    e_pa = sat_vapor_pressure_hpa(tdc) * 100.0d0
    r = mixing_ratio_from_e_pa(p_pa, e_pa)

    t_l = 56.0d0 + 1.0d0 / (1.0d0 / (td - 56.0d0) + log(t / td) / 800.0d0)
    th_l = t * (100000.0d0 / (p_pa - e_pa)) ** kappa * (t / t_l) ** (0.28d0 * r)
    thetae = th_l * exp(r * (1.0d0 + 0.448d0 * r) * (3036.0d0 / t_l - 1.78d0))
  end function thetae_bolton_like_metpy

  pure real(c_double) function dt_dp_moist(p_pa, t_k) result(dt_dp)
    real(c_double), intent(in) :: p_pa, t_k
    real(c_double) :: es, rs, num, den

    es = sat_vapor_pressure_pa_from_tk(t_k)
    rs = mixing_ratio_from_e_pa(p_pa, es)
    num = rd * t_k + lv * rs
    den = cp + (lv * lv * rs * eps) / (rd * t_k * t_k)
    dt_dp = (num / den) / p_pa
  end function dt_dp_moist

  real(c_double) function rk4_integrate_moist_t(p_start_pa, t_start_k, p_target_pa) result(t_out)
    real(c_double), intent(in) :: p_start_pa, t_start_k, p_target_pa
    real(c_double) :: dp_total, h, p, t, k1, k2, k3, k4
    integer :: nstep, istep

    dp_total = p_target_pa - p_start_pa
    if (abs(dp_total) < 1.0d-9) then
      t_out = t_start_k
      return
    end if

    nstep = int(abs(dp_total) / 200.0d0) + 1
    h = dp_total / real(nstep, c_double)
    p = p_start_pa
    t = t_start_k

    do istep = 1, nstep
      k1 = dt_dp_moist(p, t)
      k2 = dt_dp_moist(p + 0.5d0 * h, t + 0.5d0 * h * k1)
      k3 = dt_dp_moist(p + 0.5d0 * h, t + 0.5d0 * h * k2)
      k4 = dt_dp_moist(p + h, t + h * k3)
      t = t + (h / 6.0d0) * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)
      p = p + h
    end do

    t_out = t
  end function rk4_integrate_moist_t

  subroutine lcl_bolton(p_hpa, tc, tdc, plcl_hpa, tlcl_k)
    real(c_double), intent(in) :: p_hpa, tc, tdc
    real(c_double), intent(out) :: plcl_hpa, tlcl_k
    real(c_double) :: t, td

    t = tc + 273.15d0
    td = max(tdc + 273.15d0, 173.15d0)
    tlcl_k = 56.0d0 + 1.0d0 / (1.0d0 / (td - 56.0d0) + log(t / td) / 800.0d0)
    plcl_hpa = p_hpa * (tlcl_k / t) ** (1.0d0 / kappa)
  end subroutine lcl_bolton

  pure real(c_double) function interp_linear_p(p0, v0, p1, v1, p_target) result(v_out)
    real(c_double), intent(in) :: p0, v0, p1, v1, p_target
    real(c_double) :: w

    if (abs(p1 - p0) < 1.0d-12) then
      v_out = v0
      return
    end if

    if (p0 > 0.0d0 .and. p1 > 0.0d0 .and. p_target > 0.0d0) then
      if (abs(log(p1) - log(p0)) < 1.0d-12) then
        v_out = v0
        return
      end if
      w = (log(p_target) - log(p0)) / (log(p1) - log(p0))
    else
      w = (p_target - p0) / (p1 - p0)
    end if

    v_out = v0 + w * (v1 - v0)
  end function interp_linear_p

  subroutine thetae_min_in_700_500(p, t, td, m, found, start_p, start_t, start_td)
    real(c_double), intent(in) :: p(:), t(:), td(:)
    integer, intent(in) :: m
    logical, intent(out) :: found
    real(c_double), intent(out) :: start_p, start_t, start_td
    real(c_double) :: min_te, te, bound, tb, tdb
    real(c_double), dimension(2) :: bounds
    integer :: k, ibound

    min_te = 1.0d99
    start_p = ieee_value(0.0d0, ieee_quiet_nan)
    start_t = ieee_value(0.0d0, ieee_quiet_nan)
    start_td = ieee_value(0.0d0, ieee_quiet_nan)
    found = .false.
    bounds = [700.0d0, 500.0d0]

    do k = 1, m
      if (p(k) >= 500.0d0 .and. p(k) <= 700.0d0) then
        te = thetae_bolton_like_metpy(p(k), t(k), td(k))
        if (te < min_te) then
          min_te = te
          start_p = p(k)
          start_t = t(k)
          start_td = td(k)
          found = .true.
        end if
      end if
    end do

    do ibound = 1, size(bounds)
      bound = bounds(ibound)
      do k = 1, m - 1
        if ((p(k) >= bound .and. p(k + 1) <= bound) .or. &
            (p(k) <= bound .and. p(k + 1) >= bound)) then
          tb = interp_linear_p(p(k), t(k), p(k + 1), t(k + 1), bound)
          tdb = interp_linear_p(p(k), td(k), p(k + 1), td(k + 1), bound)
          te = thetae_bolton_like_metpy(bound, tb, tdb)
          if (te < min_te) then
            min_te = te
            start_p = bound
            start_t = tb
            start_td = tdb
            found = .true.
          end if
          exit
        end if
      end do
    end do
  end subroutine thetae_min_in_700_500

  subroutine dcape_profile_impl(nlev, pressure_hpa, temperature_c, dewpoint_c, dcape)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: pressure_hpa(nlev), temperature_c(nlev), dewpoint_c(nlev)
    real(c_double), intent(out) :: dcape
    real(c_double) :: p(nlev), t(nlev), td(nlev), parcel_tk(nlev)
    real(c_double) :: start_p_hpa, start_t_c, start_td_c
    real(c_double) :: plcl_hpa, tlcl_k, tw_start_k
    real(c_double) :: integral, prev_diff, prev_lnp, diff, lnp
    real(c_double) :: p_pa, e_env, r_env, tv_env, e_par, r_par, tv_par
    logical :: found, first
    integer :: m, k, start_idx, k2

    dcape = ieee_value(0.0d0, ieee_quiet_nan)
    if (nlev < 4) then
      return
    end if

    m = 0
    do k = 1, nlev
      if (ieee_is_finite(pressure_hpa(k)) .and. ieee_is_finite(temperature_c(k)) .and. &
          ieee_is_finite(dewpoint_c(k))) then
        m = m + 1
        p(m) = pressure_hpa(k)
        t(m) = temperature_c(k)
        td(m) = dewpoint_c(k)
      end if
    end do

    if (m < 4) then
      return
    end if

    if (p(1) < p(m)) then
      do k = 1, m / 2
        k2 = m + 1 - k
        call swap_real(p(k), p(k2))
        call swap_real(t(k), t(k2))
        call swap_real(td(k), td(k2))
      end do
    end if

    call thetae_min_in_700_500(p, t, td, m, found, start_p_hpa, start_t_c, start_td_c)
    if (.not. found) then
      return
    end if

    start_idx = -1
    do k = 1, m
      if (p(k) >= start_p_hpa) then
        start_idx = k
      else
        exit
      end if
    end do
    if (start_idx <= 1) then
      return
    end if

    call lcl_bolton(start_p_hpa, start_t_c, start_td_c, plcl_hpa, tlcl_k)
    tw_start_k = rk4_integrate_moist_t(plcl_hpa * 100.0d0, tlcl_k, start_p_hpa * 100.0d0)

    parcel_tk(start_idx) = rk4_integrate_moist_t( &
      start_p_hpa * 100.0d0, tw_start_k, p(start_idx) * 100.0d0 &
    )
    do k = start_idx, 2, -1
      parcel_tk(k - 1) = rk4_integrate_moist_t( &
        p(k) * 100.0d0, parcel_tk(k), p(k - 1) * 100.0d0 &
      )
    end do

    integral = 0.0d0
    first = .true.
    prev_diff = 0.0d0
    prev_lnp = 0.0d0

    do k = 1, start_idx
      p_pa = p(k) * 100.0d0
      e_env = sat_vapor_pressure_hpa(td(k)) * 100.0d0
      r_env = mixing_ratio_from_e_pa(p_pa, e_env)
      tv_env = virtual_temp_k(t(k) + 273.15d0, r_env)

      e_par = sat_vapor_pressure_pa_from_tk(parcel_tk(k))
      r_par = mixing_ratio_from_e_pa(p_pa, e_par)
      tv_par = virtual_temp_k(parcel_tk(k), r_par)

      diff = tv_env - tv_par
      lnp = log(p(k))
      if (first) then
        prev_diff = diff
        prev_lnp = lnp
        first = .false.
      else
        integral = integral + 0.5d0 * (prev_diff + diff) * (lnp - prev_lnp)
        prev_diff = diff
        prev_lnp = lnp
      end if
    end do

    dcape = -rd * integral
    if (.not. ieee_is_finite(dcape)) then
      dcape = ieee_value(0.0d0, ieee_quiet_nan)
    else if (dcape < 0.0d0) then
      dcape = 0.0d0
    end if
  end subroutine dcape_profile_impl

  subroutine dcape_grid_impl(nlev, nlat, nlon, pressure_3d, t_3d, td_3d, out)
    integer, intent(in) :: nlev, nlat, nlon
    real(c_double), intent(in) :: pressure_3d(nlev, nlat, nlon)
    real(c_double), intent(in) :: t_3d(nlev, nlat, nlon)
    real(c_double), intent(in) :: td_3d(nlev, nlat, nlon)
    real(c_double), intent(out) :: out(nlat, nlon)
    integer :: j, i

    do i = 1, nlon
      do j = 1, nlat
        call dcape_profile_impl( &
          nlev, &
          pressure_3d(:, j, i), &
          t_3d(:, j, i), &
          td_3d(:, j, i), &
          out(j, i) &
        )
      end do
    end do
  end subroutine dcape_grid_impl

  subroutine swap_real(a, b)
    real(c_double), intent(inout) :: a, b
    real(c_double) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap_real

end module dcape_mod
