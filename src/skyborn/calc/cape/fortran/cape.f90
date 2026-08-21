module cape_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
  use iso_c_binding, only: c_double
  implicit none
  private
  public :: cape_profile_impl, cape_grid_impl
  public :: parcel_profile_impl, parcel_profile_grid_impl
  public :: most_unstable_parcel_impl, most_unstable_parcel_grid_impl
  public :: mucape_profile_impl, mucape_grid_impl

  ! Physical constants -- same Bolton-based system as the DCAPE backend
  real(c_double), parameter :: eps = 0.622d0
  real(c_double), parameter :: rd = 287.05d0
  real(c_double), parameter :: cp = 1004.0d0
  real(c_double), parameter :: lv = 2.5d6
  real(c_double), parameter :: kappa = rd / cp

contains

  !--------------------------------------------------------------------------
  ! Saturation vapor pressure (hPa) from temperature in Celsius
  pure real(c_double) function sat_vapor_pressure_hpa(tc) result(es)
    real(c_double), intent(in) :: tc
    es = 6.112d0 * exp(17.67d0 * tc / (tc + 243.5d0))
  end function sat_vapor_pressure_hpa

  ! Saturation vapor pressure (Pa) from temperature in Kelvin
  pure real(c_double) function sat_vapor_pressure_pa_from_tk(tk) result(es_pa)
    real(c_double), intent(in) :: tk
    es_pa = sat_vapor_pressure_hpa(tk - 273.15d0) * 100.0d0
  end function sat_vapor_pressure_pa_from_tk

  ! Mixing ratio from vapor pressure and total pressure (both Pa)
  pure real(c_double) function mixing_ratio_from_e_pa(p_pa, e_pa) result(r)
    real(c_double), intent(in) :: p_pa, e_pa
    real(c_double) :: e
    e = max(e_pa, 1.0d-8)
    e = min(e, 0.99d0 * p_pa)
    r = eps * e / (p_pa - e)
  end function mixing_ratio_from_e_pa

  ! Virtual temperature from temperature and mixing ratio
  pure real(c_double) function virtual_temp_k(tk, mixing_ratio) result(tv)
    real(c_double), intent(in) :: tk, mixing_ratio
    tv = tk * (mixing_ratio + eps) / (eps * (1.0d0 + mixing_ratio))
  end function virtual_temp_k

  ! Moist pseudo-adiabat dT/dp (Bakhshaii2013 formulation, same as metpy)
  pure real(c_double) function dt_dp_moist(p_pa, t_k) result(dt_dp)
    real(c_double), intent(in) :: p_pa, t_k
    real(c_double) :: es, rs, num, den
    es = sat_vapor_pressure_pa_from_tk(t_k)
    rs = mixing_ratio_from_e_pa(p_pa, es)
    num = rd * t_k + lv * rs
    den = cp + (lv * lv * rs * eps) / (rd * t_k * t_k)
    dt_dp = (num / den) / p_pa
  end function dt_dp_moist

  ! Fixed-step RK4 integration of the moist adiabat
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

  ! Bolton-approximation LCL
  subroutine lcl_bolton(p_hpa, tc, tdc, plcl_hpa, tlcl_k)
    real(c_double), intent(in) :: p_hpa, tc, tdc
    real(c_double), intent(out) :: plcl_hpa, tlcl_k
    real(c_double) :: t, td
    t = tc + 273.15d0
    td = max(tdc + 273.15d0, 173.15d0)
    tlcl_k = 56.0d0 + 1.0d0 / (1.0d0 / (td - 56.0d0) + log(t / td) / 800.0d0)
    plcl_hpa = p_hpa * (tlcl_k / t) ** (1.0d0 / kappa)
  end subroutine lcl_bolton

  !--------------------------------------------------------------------------
  ! Drop non-finite levels and sort pressure in decreasing order (sfc first)
  subroutine clean_and_sort(nlev, p_in, t_in, td_in, m, p, t, td)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: p_in(nlev), t_in(nlev), td_in(nlev)
    integer, intent(out) :: m
    real(c_double), intent(out) :: p(:), t(:), td(:)
    integer :: k, k2

    m = 0
    do k = 1, nlev
      if (ieee_is_finite(p_in(k)) .and. ieee_is_finite(t_in(k)) .and. &
          ieee_is_finite(td_in(k))) then
        m = m + 1
        p(m) = p_in(k)
        t(m) = t_in(k)
        td(m) = td_in(k)
      end if
    end do

    if (m > 1) then
      if (p(1) < p(m)) then
        do k = 1, m / 2
          k2 = m + 1 - k
          call swap_real(p(k), p(k2))
          call swap_real(t(k), t(k2))
          call swap_real(td(k), td(k2))
        end do
      end if
    end if
  end subroutine clean_and_sort

  !--------------------------------------------------------------------------
  ! Parcel temperature profile (K) at each pressure level.
  ! Dry adiabatic below the LCL, moist pseudo-adiabatic (RK4) above.
  ! The moist segment is integrated sequentially level by level (O(n)
  ! rather than re-integrating from the LCL for every level).
  subroutine parcel_temp_profile(m, p_hpa, t_c, td_c, plcl_hpa, tlcl_k, t_par_k)
    integer, intent(in) :: m
    real(c_double), intent(in) :: p_hpa(m), t_c(m), td_c(m)
    real(c_double), intent(in) :: plcl_hpa, tlcl_k
    real(c_double), intent(out) :: t_par_k(m)
    real(c_double) :: t0_k, p0_hpa
    integer :: k, k2

    t0_k = t_c(1) + 273.15d0
    p0_hpa = p_hpa(1)

    ! Dry adiabat below the LCL (potential temperature conserved)
    k = 1
    do while (k <= m .and. p_hpa(k) >= plcl_hpa)
      t_par_k(k) = t0_k * (p_hpa(k) / p0_hpa) ** kappa
      k = k + 1
    end do

    ! All levels below the LCL? Nothing left to do.
    if (k > m) return

    ! First moist level: integrate once from the LCL
    t_par_k(k) = rk4_integrate_moist_t(plcl_hpa * 100.0d0, tlcl_k, &
                                       p_hpa(k) * 100.0d0)
    ! Remaining moist levels: integrate sequentially from the previous level
    do k2 = k + 1, m
      t_par_k(k2) = rk4_integrate_moist_t(p_hpa(k2 - 1) * 100.0d0, &
                                          t_par_k(k2 - 1), &
                                          p_hpa(k2) * 100.0d0)
    end do
  end subroutine parcel_temp_profile

  !--------------------------------------------------------------------------
  ! Zero crossings of y along log(p). Fills p_ext / y_ext (sorted in
  ! decreasing pressure), n_ext points. Crossing direction:
  !   dir = 1 : y goes from negative to positive
  !   dir = -1: y goes from positive to negative
  !   dir = 0 : all crossings
  subroutine append_zero_crossings(m, p, y, n_ext, p_ext, y_ext)
    integer, intent(in) :: m
    real(c_double), intent(in) :: p(m), y(m)
    integer, intent(out) :: n_ext
    real(c_double), intent(out) :: p_ext(:), y_ext(:)
    real(c_double) :: lnp0, lnp1, w, p_zero
    integer :: k

    n_ext = 0
    ! Start from the surface point (y(1) == 0 for a surface parcel)
    if (m >= 1) then
      n_ext = 1
      p_ext(1) = p(1)
      y_ext(1) = y(1)
    end if

    ! Find crossings between levels k-1 and k (metpy uses x[1:], skipping
    ! the surface point, since the surface parcel matches the environment)
    do k = 2, m
      ! sign convention: treat exact zeros as positive (as metpy does)
      if ((y(k - 1) < 0.0d0 .and. y(k) > 0.0d0) .or. &
          (y(k - 1) < 0.0d0 .and. abs(y(k)) < 1.0d-12) .or. &
          (abs(y(k - 1)) < 1.0d-12 .and. y(k) < 0.0d0)) then
        ! log-pressure linear interpolation to y = 0
        lnp0 = log(p(k - 1))
        lnp1 = log(p(k))
        if (abs(y(k) - y(k - 1)) > 1.0d-12) then
          w = (0.0d0 - y(k - 1)) / (y(k) - y(k - 1))
          p_zero = exp(lnp0 + w * (lnp1 - lnp0))
        else
          p_zero = p(k)
        end if
        ! avoid duplicate pressure values
        if (abs(p_zero - p(k)) > 1.0d-6 .and. abs(p_zero - p(k - 1)) > 1.0d-6) then
          n_ext = n_ext + 1
          p_ext(n_ext) = p_zero
          y_ext(n_ext) = 0.0d0
        end if
      else if ((y(k - 1) > 0.0d0 .and. y(k) < 0.0d0) .or. &
               (abs(y(k - 1)) < 1.0d-12 .and. y(k) > 0.0d0)) then
        lnp0 = log(p(k - 1))
        lnp1 = log(p(k))
        if (abs(y(k) - y(k - 1)) > 1.0d-12) then
          w = (0.0d0 - y(k - 1)) / (y(k) - y(k - 1))
          p_zero = exp(lnp0 + w * (lnp1 - lnp0))
        else
          p_zero = p(k)
        end if
        if (abs(p_zero - p(k)) > 1.0d-6 .and. abs(p_zero - p(k - 1)) > 1.0d-6) then
          n_ext = n_ext + 1
          p_ext(n_ext) = p_zero
          y_ext(n_ext) = 0.0d0
        end if
      end if
      n_ext = n_ext + 1
      p_ext(n_ext) = p(k)
      y_ext(n_ext) = y(k)
    end do
  end subroutine append_zero_crossings

  ! Collect crossings of a given direction. Returns count and pressures.
  subroutine collect_crossings(m, p, y, direction, n_cross, p_cross)
    integer, intent(in) :: m
    real(c_double), intent(in) :: p(m), y(m)
    integer, intent(in) :: direction   ! 1: increasing (neg->pos), -1: decreasing
    integer, intent(out) :: n_cross
    real(c_double), intent(out) :: p_cross(:)
    real(c_double) :: lnp0, lnp1, w, p_zero
    integer :: k

    n_cross = 0
    do k = 2, m
      if (direction > 0) then
        if ((y(k - 1) < 0.0d0 .and. y(k) > 0.0d0) .or. &
            (y(k - 1) < 0.0d0 .and. abs(y(k)) < 1.0d-12)) then
          lnp0 = log(p(k - 1))
          lnp1 = log(p(k))
          if (abs(y(k) - y(k - 1)) > 1.0d-12) then
            w = (0.0d0 - y(k - 1)) / (y(k) - y(k - 1))
            p_zero = exp(lnp0 + w * (lnp1 - lnp0))
          else
            p_zero = p(k)
          end if
          n_cross = n_cross + 1
          p_cross(n_cross) = p_zero
        end if
      else
        if ((y(k - 1) > 0.0d0 .and. y(k) < 0.0d0) .or. &
            (abs(y(k - 1)) < 1.0d-12 .and. y(k) < 0.0d0)) then
          lnp0 = log(p(k - 1))
          lnp1 = log(p(k))
          if (abs(y(k) - y(k - 1)) > 1.0d-12) then
            w = (0.0d0 - y(k - 1)) / (y(k) - y(k - 1))
            p_zero = exp(lnp0 + w * (lnp1 - lnp0))
          else
            p_zero = p(k)
          end if
          n_cross = n_cross + 1
          p_cross(n_cross) = p_zero
        end if
      end if
    end do
  end subroutine collect_crossings

  !--------------------------------------------------------------------------
  ! Single-profile CAPE/CIN (surface-based parcel, metpy-compatible
  ! which_lfc='bottom', which_el='top' semantics).
  subroutine cape_profile_impl(nlev, pressure_hpa, temperature_c, dewpoint_c, &
                               cape, cin, lfc_p_hpa, el_p_hpa)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: pressure_hpa(nlev)
    real(c_double), intent(in) :: temperature_c(nlev)
    real(c_double), intent(in) :: dewpoint_c(nlev)
    real(c_double), intent(out) :: cape, cin, lfc_p_hpa, el_p_hpa

    real(c_double), allocatable :: p(:), t(:), td(:)
    real(c_double), allocatable :: t_par_k(:), r_par(:), r_env(:)
    real(c_double), allocatable :: tv_par(:), tv_env(:), y(:)
    real(c_double), allocatable :: p_ext(:), y_ext(:)
    real(c_double), allocatable :: p_lfc_cross(:), p_el_cross(:)

    real(c_double) :: t0_c, td0_c, p0_hpa
    real(c_double) :: plcl_actual, tlcl_k, plcl_lfc, tv1_c
    real(c_double) :: lfc, el_p, rtol, atol
    real(c_double) :: es0_pa, es_par_pa, es_env_pa
    real(c_double) :: integral, dlnp
    integer :: m, k, n_ext, n_lfc, n_el, i
    logical :: no_lfc

    cape = ieee_value(0.0d0, ieee_quiet_nan)
    cin = ieee_value(0.0d0, ieee_quiet_nan)
    lfc_p_hpa = ieee_value(0.0d0, ieee_quiet_nan)
    el_p_hpa = ieee_value(0.0d0, ieee_quiet_nan)

    allocate(p(nlev), t(nlev), td(nlev))
    call clean_and_sort(nlev, pressure_hpa, temperature_c, dewpoint_c, m, p, t, td)

    if (m < 2) then
      deallocate(p, t, td)
      return
    end if

    allocate(t_par_k(m), r_par(m), r_env(m), tv_par(m), tv_env(m), y(m))

    ! Surface parcel state
    t0_c = t(1)
    td0_c = td(1)
    p0_hpa = p(1)

    ! LCL from actual temperature (used for parcel mixing-ratio split)
    call lcl_bolton(p0_hpa, t0_c, td0_c, plcl_actual, tlcl_k)

    ! Parcel temperature profile
    call parcel_temp_profile(m, p, t, td, plcl_actual, tlcl_k, t_par_k)

    ! Mixing ratios and virtual temperatures
    es0_pa = sat_vapor_pressure_hpa(td0_c) * 100.0d0
    do k = 1, m
      if (p(k) >= plcl_actual) then
        ! Below LCL: parcel mixing ratio fixed at surface saturated value
        r_par(k) = mixing_ratio_from_e_pa(p0_hpa * 100.0d0, es0_pa)
      else
        es_par_pa = sat_vapor_pressure_pa_from_tk(t_par_k(k))
        r_par(k) = mixing_ratio_from_e_pa(p(k) * 100.0d0, es_par_pa)
      end if
      es_env_pa = sat_vapor_pressure_hpa(td(k)) * 100.0d0
      r_env(k) = mixing_ratio_from_e_pa(p(k) * 100.0d0, es_env_pa)
      tv_par(k) = virtual_temp_k(t_par_k(k), r_par(k))
      tv_env(k) = virtual_temp_k(t(k) + 273.15d0, r_env(k))
      y(k) = tv_par(k) - tv_env(k)
    end do

    ! LCL in virtual temperature, used by metpy's lfc()/el() as the
    ! lower bound for physically meaningful intersections
    tv1_c = tv_par(1) - 273.15d0
    call lcl_bolton(p0_hpa, tv1_c, td0_c, plcl_lfc, tlcl_k)

    ! ---- Level of free convection (metpy which='bottom') ----
    no_lfc = .false.
    lfc = ieee_value(0.0d0, ieee_quiet_nan)
    allocate(p_lfc_cross(m), p_el_cross(m))
    call collect_crossings(m, p, y, 1, n_lfc, p_lfc_cross)
    if (n_lfc == 0) then
      ! Any positive area above the LCL?
      do k = 1, m
        if (p(k) < plcl_lfc .and. y(k) > 0.0d0) exit
      end do
      if (k > m) then
        no_lfc = .true.   ! no positive buoyancy above LCL: no LFC
      else
        lfc = plcl_lfc    ! LFC = LCL
      end if
    else
      ! Is any crossing above the LCL (pressure < LCL)?
      do i = 1, n_lfc
        if (p_lfc_cross(i) < plcl_lfc) exit
      end do
      if (i > n_lfc) then
        ! All crossings below the LCL. Check whether an EL also exists
        ! below the LCL: if so, no LFC; otherwise LFC = LCL.
        call collect_crossings(m, p, y, -1, n_el, p_el_cross)
        if (n_el > 0) then
          if (minval(p_el_cross) > plcl_lfc) then
            no_lfc = .true.
          else
            lfc = plcl_lfc
          end if
        else
          lfc = plcl_lfc
        end if
      else
        ! First (highest-pressure) crossing above the LCL
        lfc = p_lfc_cross(i)
      end if
    end if

    if (no_lfc) then
      cape = 0.0d0
      cin = 0.0d0
      lfc_p_hpa = ieee_value(0.0d0, ieee_quiet_nan)
      el_p_hpa = ieee_value(0.0d0, ieee_quiet_nan)
      deallocate(p, t, td, t_par_k, r_par, r_env, tv_par, tv_env, y, &
                 p_lfc_cross, p_el_cross)
      return
    end if

    ! ---- Equilibrium level (metpy which='top') ----
    el_p = ieee_value(0.0d0, ieee_quiet_nan)
    if (y(m) > 0.0d0) then
      ! Parcel warmer than environment at the top of the sounding
      el_p = p(m)
    else
      call collect_crossings(m, p, y, -1, n_el, p_el_cross)
      if (n_el > 0) then
        ! ELs below the LCL; take the lowest pressure one (which='top')
        do i = n_el, 1, -1
          if (p_el_cross(i) < plcl_lfc) exit
        end do
        if (i >= 1) then
          el_p = p_el_cross(i)
        else
          el_p = p(m)
        end if
      else
        el_p = p(m)
      end if
    end if

    ! ---- Build extended series with appended zero crossings ----
    allocate(p_ext(2 * m + 2), y_ext(2 * m + 2))
    call append_zero_crossings(m, p, y, n_ext, p_ext, y_ext)

    ! ---- Integrate (pressure decreasing order) ----
    rtol = 1.0d-5
    atol = 1.0d-8

    ! CIN: from surface to LFC.  cin = -Rd * sum 0.5*(y_i+y_{i+1})*dlnp
    integral = 0.0d0
    do i = 1, n_ext - 1
      ! segment [p_ext(i), p_ext(i+1)], both >= LFC (isclose semantics)
      if ((p_ext(i) >= lfc .or. abs(p_ext(i) - lfc) <= atol + rtol * abs(lfc)) .and. &
          (p_ext(i + 1) >= lfc .or. abs(p_ext(i + 1) - lfc) <= atol + rtol * abs(lfc))) then
        dlnp = log(p_ext(i + 1)) - log(p_ext(i))
        integral = integral + 0.5d0 * (y_ext(i) + y_ext(i + 1)) * dlnp
      end if
    end do
    cin = -rd * integral
    if (cin > 0.0d0) cin = 0.0d0

    ! CAPE: from LFC to EL
    integral = 0.0d0
    do i = 1, n_ext - 1
      if ((p_ext(i) <= lfc .or. abs(p_ext(i) - lfc) <= atol + rtol * abs(lfc)) .and. &
          (p_ext(i) >= el_p .or. abs(p_ext(i) - el_p) <= atol + rtol * abs(el_p)) .and. &
          (p_ext(i + 1) <= lfc .or. abs(p_ext(i + 1) - lfc) <= atol + rtol * abs(lfc)) .and. &
          (p_ext(i + 1) >= el_p .or. abs(p_ext(i + 1) - el_p) <= atol + rtol * abs(el_p))) then
        dlnp = log(p_ext(i + 1)) - log(p_ext(i))
        integral = integral + 0.5d0 * (y_ext(i) + y_ext(i + 1)) * dlnp
      end if
    end do
    cape = -rd * integral
    if (.not. ieee_is_finite(cape)) cape = 0.0d0
    if (.not. ieee_is_finite(cin)) cin = 0.0d0

    lfc_p_hpa = lfc
    el_p_hpa = el_p

    deallocate(p, t, td, t_par_k, r_par, r_env, tv_par, tv_env, y, &
               p_ext, y_ext, p_lfc_cross, p_el_cross)
  end subroutine cape_profile_impl

  !--------------------------------------------------------------------------
  ! Single-profile parcel temperature (Celsius), for reuse by callers.
  subroutine parcel_profile_impl(nlev, pressure_hpa, temperature_c, dewpoint_c, t_par_c)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: pressure_hpa(nlev)
    real(c_double), intent(in) :: temperature_c(nlev)
    real(c_double), intent(in) :: dewpoint_c(nlev)
    real(c_double), intent(out) :: t_par_c(nlev)

    real(c_double), allocatable :: p(:), t(:), td(:), t_par_k(:)
    real(c_double) :: plcl_hpa, tlcl_k
    integer :: m, k, j

    t_par_c = ieee_value(0.0d0, ieee_quiet_nan)
    allocate(p(nlev), t(nlev), td(nlev), t_par_k(nlev))
    call clean_and_sort(nlev, pressure_hpa, temperature_c, dewpoint_c, m, p, t, td)
    if (m < 2) then
      deallocate(p, t, td, t_par_k)
      return
    end if

    call lcl_bolton(p(1), t(1), td(1), plcl_hpa, tlcl_k)
    call parcel_temp_profile(m, p, t, td, plcl_hpa, tlcl_k, t_par_k)

    ! Map back to the original (cleaned) level ordering
    do k = 1, nlev
      do j = 1, m
        if (abs(pressure_hpa(k) - p(j)) < 1.0d-9 .and. &
            abs(temperature_c(k) - t(j)) < 1.0d-9 .and. &
            abs(dewpoint_c(k) - td(j)) < 1.0d-9) then
          t_par_c(k) = t_par_k(j) - 273.15d0
          exit
        end if
      end do
    end do
    deallocate(p, t, td, t_par_k)
  end subroutine parcel_profile_impl

  !--------------------------------------------------------------------------
  ! Grid kernels: loop over horizontal profiles (column-major friendly)
  subroutine parcel_profile_grid_impl(nlev, nlat, nlon, p3, t3, td3, out)
    integer, intent(in) :: nlev, nlat, nlon
    real(c_double), intent(in) :: p3(nlev, nlat, nlon)
    real(c_double), intent(in) :: t3(nlev, nlat, nlon)
    real(c_double), intent(in) :: td3(nlev, nlat, nlon)
    real(c_double), intent(out) :: out(nlev, nlat, nlon)
    integer :: j, i
    do i = 1, nlon
      do j = 1, nlat
        call parcel_profile_impl(nlev, p3(:, j, i), t3(:, j, i), td3(:, j, i), out(:, j, i))
      end do
    end do
  end subroutine parcel_profile_grid_impl

  subroutine cape_grid_impl(nlev, nlat, nlon, p3, t3, td3, cape2, cin2, lfc2, el2)
    integer, intent(in) :: nlev, nlat, nlon
    real(c_double), intent(in) :: p3(nlev, nlat, nlon)
    real(c_double), intent(in) :: t3(nlev, nlat, nlon)
    real(c_double), intent(in) :: td3(nlev, nlat, nlon)
    real(c_double), intent(out) :: cape2(nlat, nlon)
    real(c_double), intent(out) :: cin2(nlat, nlon)
    real(c_double), intent(out) :: lfc2(nlat, nlon)
    real(c_double), intent(out) :: el2(nlat, nlon)
    integer :: j, i
    do i = 1, nlon
      do j = 1, nlat
        call cape_profile_impl(nlev, p3(:, j, i), t3(:, j, i), td3(:, j, i), &
                               cape2(j, i), cin2(j, i), lfc2(j, i), el2(j, i))
      end do
    end do
  end subroutine cape_grid_impl

  subroutine swap_real(a, b)
    real(c_double), intent(inout) :: a, b
    real(c_double) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap_real

  !--------------------------------------------------------------------------
  ! Bolton 1980 equivalent potential temperature (K). Identical to metpy's
  ! equivalent_potential_temperature and to the DCAPE backend's formulation.
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

  !--------------------------------------------------------------------------
  ! Search the layer [p(1)-depth, p(1)] for the level with maximum theta-e.
  ! Operates on a cleaned, pressure-decreasing profile. Returns the level
  ! index k0 (1-based within the cleaned array) or -1 if none qualify.
  subroutine most_unstable_search(m, p, t, td, depth_hpa, k0)
    integer, intent(in) :: m
    real(c_double), intent(in) :: p(m), t(m), td(m)
    real(c_double), intent(in) :: depth_hpa
    integer, intent(out) :: k0
    real(c_double) :: bottom_p, top_p, te_max, te
    integer :: k

    k0 = -1
    if (m < 1) return
    bottom_p = p(1)
    top_p = bottom_p - depth_hpa
    te_max = -huge(0.0d0)
    do k = 1, m
      if (p(k) <= bottom_p .and. p(k) >= top_p - 1.0d-9) then
        te = thetae_bolton_like_metpy(p(k), t(k), td(k))
        if (te > te_max) then
          te_max = te
          k0 = k
        end if
      end if
    end do
  end subroutine most_unstable_search

  !--------------------------------------------------------------------------
  ! Most unstable parcel: pressure / temperature / dewpoint of the level with
  ! the largest theta-e inside the bottom-`depth` layer (metpy default 300 hPa).
  subroutine most_unstable_parcel_impl(nlev, pressure_hpa, temperature_c, dewpoint_c, &
                                       depth_hpa, mup_p, mup_t, mup_td, mup_idx)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: pressure_hpa(nlev)
    real(c_double), intent(in) :: temperature_c(nlev)
    real(c_double), intent(in) :: dewpoint_c(nlev)
    real(c_double), intent(in) :: depth_hpa
    real(c_double), intent(out) :: mup_p, mup_t, mup_td
    integer, intent(out) :: mup_idx

    real(c_double), allocatable :: p(:), t(:), td(:)
    integer :: m, k0

    mup_p = ieee_value(0.0d0, ieee_quiet_nan)
    mup_t = ieee_value(0.0d0, ieee_quiet_nan)
    mup_td = ieee_value(0.0d0, ieee_quiet_nan)
    mup_idx = -1

    allocate(p(nlev), t(nlev), td(nlev))
    call clean_and_sort(nlev, pressure_hpa, temperature_c, dewpoint_c, m, p, t, td)
    if (m < 1) then
      deallocate(p, t, td)
      return
    end if
    call most_unstable_search(m, p, t, td, depth_hpa, k0)
    if (k0 >= 1) then
      mup_p = p(k0)
      mup_t = t(k0)
      mup_td = td(k0)
      mup_idx = k0 - 1   ! 0-based index within the cleaned profile
    end if
    deallocate(p, t, td)
  end subroutine most_unstable_parcel_impl

  !--------------------------------------------------------------------------
  ! Most unstable CAPE/CIN: parcel lifted from the most unstable level.
  subroutine mucape_profile_impl(nlev, pressure_hpa, temperature_c, dewpoint_c, &
                                 depth_hpa, cape, cin, lfc_p, el_p)
    integer, intent(in) :: nlev
    real(c_double), intent(in) :: pressure_hpa(nlev)
    real(c_double), intent(in) :: temperature_c(nlev)
    real(c_double), intent(in) :: dewpoint_c(nlev)
    real(c_double), intent(in) :: depth_hpa
    real(c_double), intent(out) :: cape, cin, lfc_p, el_p

    real(c_double), allocatable :: p(:), t(:), td(:)
    real(c_double), allocatable :: sp(:), st(:), std(:)
    integer :: m, k0

    cape = ieee_value(0.0d0, ieee_quiet_nan)
    cin = ieee_value(0.0d0, ieee_quiet_nan)
    lfc_p = ieee_value(0.0d0, ieee_quiet_nan)
    el_p = ieee_value(0.0d0, ieee_quiet_nan)

    allocate(p(nlev), t(nlev), td(nlev))
    call clean_and_sort(nlev, pressure_hpa, temperature_c, dewpoint_c, m, p, t, td)
    if (m < 2) then
      deallocate(p, t, td)
      return
    end if
    call most_unstable_search(m, p, t, td, depth_hpa, k0)
    if (k0 < 1) then
      deallocate(p, t, td)
      return
    end if

    allocate(sp(m - k0 + 1), st(m - k0 + 1), std(m - k0 + 1))
    sp(1:) = p(k0:m)
    st(1:) = t(k0:m)
    std(1:) = td(k0:m)
    call cape_profile_impl(m - k0 + 1, sp, st, std, cape, cin, lfc_p, el_p)
    deallocate(p, t, td, sp, st, std)
  end subroutine mucape_profile_impl

  !--------------------------------------------------------------------------
  ! Grid versions
  subroutine most_unstable_parcel_grid_impl(nlev, nlat, nlon, p3, t3, td3, depth, &
                                            out_p3, out_t3, out_td3, out_idx3)
    integer, intent(in) :: nlev, nlat, nlon
    real(c_double), intent(in) :: p3(nlev, nlat, nlon)
    real(c_double), intent(in) :: t3(nlev, nlat, nlon)
    real(c_double), intent(in) :: td3(nlev, nlat, nlon)
    real(c_double), intent(in) :: depth
    real(c_double), intent(out) :: out_p3(nlev, nlat, nlon)
    real(c_double), intent(out) :: out_t3(nlev, nlat, nlon)
    real(c_double), intent(out) :: out_td3(nlev, nlat, nlon)
    integer, intent(out) :: out_idx3(nlev, nlat, nlon)
    real(c_double) :: mp, mt, mtd
    integer :: mi, j, i

    out_p3 = ieee_value(0.0d0, ieee_quiet_nan)
    out_t3 = ieee_value(0.0d0, ieee_quiet_nan)
    out_td3 = ieee_value(0.0d0, ieee_quiet_nan)
    out_idx3 = -1
    do i = 1, nlon
      do j = 1, nlat
        call most_unstable_parcel_impl(nlev, p3(:, j, i), t3(:, j, i), td3(:, j, i), &
                                       depth, mp, mt, mtd, mi)
        if (mi >= 0) then
          out_p3(mi + 1, j, i) = mp
          out_t3(mi + 1, j, i) = mt
          out_td3(mi + 1, j, i) = mtd
          out_idx3(mi + 1, j, i) = mi
        end if
      end do
    end do
  end subroutine most_unstable_parcel_grid_impl

  subroutine mucape_grid_impl(nlev, nlat, nlon, p3, t3, td3, depth, cape2, cin2, lfc2, el2)
    integer, intent(in) :: nlev, nlat, nlon
    real(c_double), intent(in) :: p3(nlev, nlat, nlon)
    real(c_double), intent(in) :: t3(nlev, nlat, nlon)
    real(c_double), intent(in) :: td3(nlev, nlat, nlon)
    real(c_double), intent(in) :: depth
    real(c_double), intent(out) :: cape2(nlat, nlon)
    real(c_double), intent(out) :: cin2(nlat, nlon)
    real(c_double), intent(out) :: lfc2(nlat, nlon)
    real(c_double), intent(out) :: el2(nlat, nlon)
    integer :: j, i
    do i = 1, nlon
      do j = 1, nlat
        call mucape_profile_impl(nlev, p3(:, j, i), t3(:, j, i), td3(:, j, i), &
                                 depth, cape2(j, i), cin2(j, i), lfc2(j, i), el2(j, i))
      end do
    end do
  end subroutine mucape_grid_impl

end module cape_mod
