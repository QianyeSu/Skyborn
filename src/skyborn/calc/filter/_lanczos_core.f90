!
! Lanczos filter Fortran 90 implementation
! Optimized for SIMD and modern processors
! Based on Duchon (1979) - Lanczos Filtering in One and Two Dimensions
!

! Compute Lanczos filter weights
subroutine compute_lanczos_weights(cutoff_freq, nwt, pass_type, weights, ierr)
    implicit none

    ! Arguments
    real(8), intent(in) :: cutoff_freq
    integer, intent(in) :: nwt
    integer, intent(in) :: pass_type
    real(8), intent(out) :: weights(nwt)
    integer, intent(out) :: ierr

    ! Local variables
    integer :: k, mwt
    real(8) :: pi, fcon, fact, sumw, x, sigma

    ierr = 0
    pi = 4.0d0 * atan(1.0d0)

    ! Validate inputs
    if (nwt < 3) then
        ierr = 1
        return
    end if

    if (mod(nwt, 2) == 0) then
        ierr = 2
        return
    end if

    if (cutoff_freq <= 0.0d0 .or. cutoff_freq >= 0.5d0) then
        ierr = 3
        return
    end if

    mwt = (nwt - 1) / 2
    fcon = 2.0d0 * cutoff_freq

    ! Compute weights using Lanczos window
    do k = -mwt, mwt
        if (k == 0) then
            weights(mwt + 1) = fcon
        else
            x = real(k, 8)
            fact = sin(fcon * pi * x) / (pi * x)
            sigma = sin(pi * x / real(mwt, 8)) / (pi * x / real(mwt, 8))
            weights(k + mwt + 1) = fact * sigma
        end if
    end do

    ! Normalize weights
    if (pass_type == 1) then
        ! Low-pass: sum to 1
        sumw = sum(weights)
        if (abs(sumw) > 1.0d-10) then
            weights = weights / sumw
        else
            ierr = 4
            return
        end if
    else if (pass_type == 2) then
        ! High-pass: subtract from delta
        sumw = sum(weights)
        if (abs(sumw) > 1.0d-10) then
            weights = weights / sumw
        else
            ierr = 4
            return
        end if
        weights = -weights
        weights(mwt + 1) = weights(mwt + 1) + 1.0d0
    else
        ierr = 5
        return
    end if

end subroutine compute_lanczos_weights


! Apply 1D filter with reflect boundary
subroutine apply_filter_1d(data, n, weights, nwt, filtered)
    implicit none

    ! Arguments
    integer, intent(in) :: n, nwt
    real(8), intent(in) :: data(n)
    real(8), intent(in) :: weights(nwt)
    real(8), intent(out) :: filtered(n)

    ! Local variables
    integer :: i, j, k, mwt, idx
    real(8) :: sum_val

    mwt = (nwt - 1) / 2

    ! Apply convolution with reflect boundary
    !$OMP PARALLEL DO PRIVATE(i, j, k, idx, sum_val) SHARED(data, weights, filtered, n, nwt, mwt)
    do i = 1, n
        sum_val = 0.0d0

        do j = -mwt, mwt
            k = i + j

            ! Reflect boundary
            if (k < 1) then
                idx = 2 - k
            else if (k > n) then
                idx = 2 * n - k
            else
                idx = k
            end if

            sum_val = sum_val + data(idx) * weights(j + mwt + 1)
        end do

        filtered(i) = sum_val
    end do
    !$OMP END PARALLEL DO

end subroutine apply_filter_1d


! Apply 1D filter with periodic boundary
subroutine apply_filter_1d_periodic(data, n, weights, nwt, filtered)
    implicit none

    ! Arguments
    integer, intent(in) :: n, nwt
    real(8), intent(in) :: data(n)
    real(8), intent(in) :: weights(nwt)
    real(8), intent(out) :: filtered(n)

    ! Local variables
    integer :: i, j, k, mwt, idx
    real(8) :: sum_val

    mwt = (nwt - 1) / 2

    ! Apply convolution with periodic boundary
    !$OMP PARALLEL DO PRIVATE(i, j, k, idx, sum_val) SHARED(data, weights, filtered, n, nwt, mwt)
    do i = 1, n
        sum_val = 0.0d0

        do j = -mwt, mwt
            k = i + j
            idx = modulo(k - 1, n) + 1
            sum_val = sum_val + data(idx) * weights(j + mwt + 1)
        end do

        filtered(i) = sum_val
    end do
    !$OMP END PARALLEL DO

end subroutine apply_filter_1d_periodic


! Apply 2D separable filter
subroutine apply_filter_2d(data, ny, nx, weights_y, nwt_y, weights_x, nwt_x, filtered)
    implicit none

    ! Arguments
    integer, intent(in) :: ny, nx, nwt_y, nwt_x
    real(8), intent(in) :: data(ny, nx)
    real(8), intent(in) :: weights_y(nwt_y), weights_x(nwt_x)
    real(8), intent(out) :: filtered(ny, nx)

    ! Local variables
    real(8) :: temp(ny, nx)
    integer :: i, j

    ! Filter along y-direction
    !$OMP PARALLEL DO PRIVATE(j) SHARED(data, temp, weights_y, ny, nx, nwt_y)
    do j = 1, nx
        call apply_filter_1d(data(:, j), ny, weights_y, nwt_y, temp(:, j))
    end do
    !$OMP END PARALLEL DO

    ! Filter along x-direction
    !$OMP PARALLEL DO PRIVATE(i) SHARED(temp, filtered, weights_x, ny, nx, nwt_x)
    do i = 1, ny
        call apply_filter_1d(temp(i, :), nx, weights_x, nwt_x, filtered(i, :))
    end do
    !$OMP END PARALLEL DO

end subroutine apply_filter_2d


! Apply 3D separable filter
subroutine apply_filter_3d(data, nz, ny, nx, weights_z, nwt_z, weights_y, nwt_y, weights_x, nwt_x, filtered)
    implicit none

    ! Arguments
    integer, intent(in) :: nz, ny, nx, nwt_z, nwt_y, nwt_x
    real(8), intent(in) :: data(nz, ny, nx)
    real(8), intent(in) :: weights_z(nwt_z), weights_y(nwt_y), weights_x(nwt_x)
    real(8), intent(out) :: filtered(nz, ny, nx)

    ! Local variables
    real(8) :: temp1(nz, ny, nx), temp2(nz, ny, nx)
    integer :: i, j, k

    ! Filter along z-direction
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j) SHARED(data, temp1, weights_z, nz, ny, nx, nwt_z)
    do j = 1, nx
        do i = 1, ny
            call apply_filter_1d(data(:, i, j), nz, weights_z, nwt_z, temp1(:, i, j))
        end do
    end do
    !$OMP END PARALLEL DO

    ! Filter along y-direction
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(k, j) SHARED(temp1, temp2, weights_y, nz, ny, nx, nwt_y)
    do j = 1, nx
        do k = 1, nz
            call apply_filter_1d(temp1(k, :, j), ny, weights_y, nwt_y, temp2(k, :, j))
        end do
    end do
    !$OMP END PARALLEL DO

    ! Filter along x-direction
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(k, i) SHARED(temp2, filtered, weights_x, nz, ny, nx, nwt_x)
    do i = 1, ny
        do k = 1, nz
            call apply_filter_1d(temp2(k, i, :), nx, weights_x, nwt_x, filtered(k, i, :))
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine apply_filter_3d
