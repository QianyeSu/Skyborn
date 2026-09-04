module causality_core_mod
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    public :: liang_single
    public :: liang_batch
    public :: ar1_filter_batch

contains

subroutine liang_single(y1, y2, nm, npt, t21, tau21, z, dh1_star, dh1_noise, ierr)
    integer, intent(in) :: nm, npt
    real(real64), intent(in) :: y1(nm), y2(nm)
    real(real64), intent(out) :: t21, tau21, z, dh1_star, dh1_noise
    integer, intent(out) :: ierr

    integer :: i, n
    real(real64) :: mean_y1, mean_y2, mean_g1, mean_g2
    real(real64) :: c11, c12, c21, c22, det_c
    real(real64) :: dc11, dc12, dc21, dc22
    real(real64) :: a11, a12, f1, residual, q1, b1
    real(real64) :: inv_npt, inv_n, inv_nm1

    t21 = ieee_value(0.0_real64, ieee_quiet_nan)
    tau21 = t21
    z = t21
    dh1_star = t21
    dh1_noise = t21
    ierr = 0

    if (npt < 1 .or. nm - npt < 2) then
        ierr = 1
        return
    end if

    n = nm - npt
    do i = 1, nm
        if (.not. ieee_is_finite(y1(i)) .or. .not. ieee_is_finite(y2(i))) then
            ierr = 2
            return
        end if
    end do

    inv_n = 1.0_real64 / real(n, real64)
    inv_nm1 = 1.0_real64 / real(n - 1, real64)
    inv_npt = 1.0_real64 / real(npt, real64)

    mean_y1 = 0.0_real64
    mean_y2 = 0.0_real64
    mean_g1 = 0.0_real64
    mean_g2 = 0.0_real64
    do i = 1, n
        mean_y1 = mean_y1 + y1(i)
        mean_y2 = mean_y2 + y2(i)
        mean_g1 = mean_g1 + (y1(i + npt) - y1(i)) * inv_npt
        mean_g2 = mean_g2 + (y2(i + npt) - y2(i)) * inv_npt
    end do
    mean_y1 = mean_y1 * inv_n
    mean_y2 = mean_y2 * inv_n
    mean_g1 = mean_g1 * inv_n
    mean_g2 = mean_g2 * inv_n

    c11 = 0.0_real64
    c12 = 0.0_real64
    c21 = 0.0_real64
    c22 = 0.0_real64
    dc11 = 0.0_real64
    dc12 = 0.0_real64
    dc21 = 0.0_real64
    dc22 = 0.0_real64
    do i = 1, n
        c11 = c11 + (y1(i) - mean_y1) * (y1(i) - mean_y1)
        c12 = c12 + (y1(i) - mean_y1) * (y2(i) - mean_y2)
        c21 = c21 + (y2(i) - mean_y2) * (y1(i) - mean_y1)
        c22 = c22 + (y2(i) - mean_y2) * (y2(i) - mean_y2)
        dc11 = dc11 + (y1(i) - mean_y1) * ( &
            (y1(i + npt) - y1(i)) * inv_npt - mean_g1 &
        )
        dc12 = dc12 + (y1(i) - mean_y1) * ( &
            (y2(i + npt) - y2(i)) * inv_npt - mean_g2 &
        )
        dc21 = dc21 + (y2(i) - mean_y2) * ( &
            (y1(i + npt) - y1(i)) * inv_npt - mean_g1 &
        )
        dc22 = dc22 + (y2(i) - mean_y2) * ( &
            (y2(i + npt) - y2(i)) * inv_npt - mean_g2 &
        )
    end do
    c11 = c11 * inv_nm1
    c12 = c12 * inv_nm1
    c21 = c21 * inv_nm1
    c22 = c22 * inv_nm1
    dc11 = dc11 * inv_nm1
    dc12 = dc12 * inv_nm1
    dc21 = dc21 * inv_nm1
    dc22 = dc22 * inv_nm1

    det_c = c11 * c22 - c12 * c21
    if (c11 == 0.0_real64 .or. det_c == 0.0_real64) then
        ierr = 3
        return
    end if

    a11 = (c22 * dc11 - c12 * dc21) / det_c
    a12 = (-c12 * dc11 + c11 * dc21) / det_c
    t21 = (c12 / c11) * (-c21 * dc11 + c11 * dc21) / det_c

    f1 = mean_g1 - a11 * mean_y1 - a12 * mean_y2
    q1 = 0.0_real64
    do i = 1, n
        residual = (y1(i + npt) - y1(i)) * inv_npt - &
            (f1 + a11 * y1(i) + a12 * y2(i))
        q1 = q1 + residual * residual
    end do

    b1 = sqrt(q1 * inv_n)
    dh1_star = a11
    dh1_noise = b1 * b1 / (2.0_real64 * c11)
    z = abs(t21) + abs(dh1_star) + abs(dh1_noise)
    if (z == 0.0_real64 .or. .not. ieee_is_finite(z)) then
        ierr = 4
        return
    end if

    tau21 = t21 / z
    dh1_star = dh1_star / z
    dh1_noise = dh1_noise / z
end subroutine liang_single


subroutine liang_batch(y1, y2, nm, nsim, npt, t21, tau21, ierr)
    integer, intent(in) :: nm, nsim, npt
    real(real64), intent(in) :: y1(nm, nsim), y2(nm, nsim)
    real(real64), intent(out) :: t21(nsim), tau21(nsim)
    integer, intent(out) :: ierr

    integer :: i, local_ierr
    real(real64) :: z, dh1_star, dh1_noise

    ierr = 0
    !$omp parallel do private(i, local_ierr, z, dh1_star, dh1_noise) &
    !$omp& if(nsim >= 32) schedule(static)
    do i = 1, nsim
        call liang_single( &
            y1(:, i), y2(:, i), nm, npt, t21(i), tau21(i), &
            z, dh1_star, dh1_noise, local_ierr &
        )
        if (local_ierr /= 0) then
            t21(i) = ieee_value(0.0_real64, ieee_quiet_nan)
            tau21(i) = t21(i)
            !$omp critical(causality_error_code)
            if (ierr == 0) ierr = local_ierr
            !$omp end critical(causality_error_code)
        end if
    end do
    !$omp end parallel do
end subroutine liang_batch


subroutine ar1_filter_batch(innovations, output, nnoise, nsim, nout, burnin, g, ierr)
    integer, intent(in) :: nnoise, nsim, nout, burnin
    real(real64), intent(in) :: innovations(nnoise, nsim), g
    real(real64), intent(out) :: output(nout, nsim)
    integer, intent(out) :: ierr

    integer :: i, j
    real(real64) :: state

    ierr = 0
    if (nnoise < 1 .or. nsim < 1 .or. nout < 1 .or. burnin < 0 .or. &
        burnin + nout > nnoise) then
        ierr = 1
        return
    end if

    !$omp parallel do private(i, j, state) if(nsim >= 32) schedule(static)
    do j = 1, nsim
        state = 0.0_real64
        do i = 1, nnoise
            state = g * state + innovations(i, j)
            if (i > burnin .and. i <= burnin + nout) then
                output(i - burnin, j) = state
            end if
        end do
    end do
    !$omp end parallel do
end subroutine ar1_filter_batch

end module causality_core_mod
