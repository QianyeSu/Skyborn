! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-19
! File    : mk_kernels_backend.f90
! Purpose : Provide the private optional Fortran backend for the clean-series
!           Mann-Kendall hot path in `skyborn.calc.mann_kendall`.
! Notes   : Public batch entry points operate on `data(ntime, nseries)` with
!           no NaN values. Missing-value filtering and fallback dispatch stay
!           on the Python side so this backend can focus on dense numeric work.
! =============================================================================
!
module mk_kernels_core
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    public :: real64
    public :: compute_base_variance
    public :: compute_modified_variance
    public :: compute_modified_variance_with_slope
    public :: compute_s_value
    public :: compute_sen_slope

contains

    ! QUICK REFERENCE
    ! PURPOSE
    !    CHOOSE A MEDIAN-OF-THREE PIVOT VALUE FOR THE IN-PLACE PARTITION
    !    ROUTINES USED BY QUICKSORT AND QUICKSELECT.
    !
    ! INPUTS
    !    A, B, C - THREE CANDIDATE VALUES
    !
    ! OUTPUT
    !    VALUE - THE MEDIAN OF THE THREE INPUTS
    pure real(real64) function median_of_three(a, b, c) result(value)
        real(real64), intent(in) :: a, b, c

        if (a < b) then
            if (b < c) then
                value = b
            else if (a < c) then
                value = c
            else
                value = a
            end if
        else
            if (a < c) then
                value = a
            else if (b < c) then
                value = c
            else
                value = b
            end if
        end if
    end function median_of_three


    ! QUICK REFERENCE
    ! PURPOSE
    !    SORT A SMALL SLICE OF A REAL BUFFER IN ASCENDING ORDER USING
    !    INSERTION SORT.
    !
    ! INPUTS
    !    LEFT, RIGHT - 1-BASED INCLUSIVE BOUNDS OF THE ACTIVE WINDOW
    !
    ! INPUT / OUTPUT
    !    VALUES(:) - BUFFER TO SORT IN PLACE OVER VALUES(LEFT:RIGHT)
    subroutine insertion_sort_real(values, left, right)
        real(real64), intent(inout) :: values(:)
        integer, intent(in) :: left, right

        integer :: i, j
        real(real64) :: current

        do i = left + 1, right
            current = values(i)
            j = i - 1
            do while (j >= left .and. values(j) > current)
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = current
        end do
    end subroutine insertion_sort_real


    ! QUICK REFERENCE
    ! PURPOSE
    !    RECURSIVELY PARTITION AND SORT A REAL BUFFER IN ASCENDING ORDER.
    !
    ! INPUTS
    !    LEFT, RIGHT - 1-BASED INCLUSIVE BOUNDS OF THE ACTIVE WINDOW
    !
    ! INPUT / OUTPUT
    !    VALUES(:) - BUFFER TO SORT IN PLACE OVER VALUES(LEFT:RIGHT)
    !
    ! NOTES
    !    SMALL TAILS ARE LEFT FOR `insertion_sort_real(...)` TO FINISH
    !    BECAUSE THAT IS FASTER THAN RECURSING ALL THE WAY DOWN.
    recursive subroutine quicksort_real(values, left, right)
        real(real64), intent(inout) :: values(:)
        integer, intent(in) :: left, right

        integer, parameter :: insertion_threshold = 16
        integer :: i, j
        real(real64) :: pivot, temp

        if (right - left <= insertion_threshold) return

        i = left
        j = right
        pivot = values(left + (right - left) / 2)

        do
            do while (values(i) < pivot)
                i = i + 1
            end do

            do while (values(j) > pivot)
                j = j - 1
            end do

            if (i <= j) then
                temp = values(i)
                values(i) = values(j)
                values(j) = temp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (left < j) call quicksort_real(values, left, j)
        if (i < right) call quicksort_real(values, i, right)
    end subroutine quicksort_real


    ! QUICK REFERENCE
    ! PURPOSE
    !    SORT AN ENTIRE REAL BUFFER IN ASCENDING ORDER.
    !
    ! INPUT / OUTPUT
    !    VALUES(:) - FULL BUFFER TO SORT IN PLACE
    subroutine sort_real_inplace(values)
        real(real64), intent(inout) :: values(:)

        if (size(values) <= 1) return
        call quicksort_real(values, 1, size(values))
        call insertion_sort_real(values, 1, size(values))
    end subroutine sort_real_inplace


    ! QUICK REFERENCE
    ! PURPOSE
    !    PARTITION A REAL BUFFER IN PLACE UNTIL THE REQUESTED ORDER
    !    STATISTIC IS POSITIONED AT `KTH_INDEX`.
    !
    ! INPUTS
    !    KTH_INDEX - 1-BASED ORDER STATISTIC TO SELECT
    !
    ! INPUT / OUTPUT
    !    VALUES(:) - BUFFER TO PARTITION IN PLACE
    !
    ! NOTES
    !    THIS AVOIDS A FULL SORT WHEN ONLY THE MEDIAN OR NEAR-MEDIAN
    !    PAIRWISE SLOPE IS NEEDED.
    subroutine select_kth_real(values, kth_index)
        real(real64), intent(inout) :: values(:)
        integer, intent(in) :: kth_index

        integer, parameter :: insertion_threshold = 16
        integer :: i, j, left, right
        real(real64) :: pivot, temp

        left = 1
        right = size(values)

        do while (right - left > insertion_threshold)
            pivot = median_of_three(values(left), values(left + (right - left) / 2), values(right))
            i = left
            j = right

            do
                do while (values(i) < pivot)
                    i = i + 1
                end do

                do while (values(j) > pivot)
                    j = j - 1
                end do

                if (i <= j) then
                    temp = values(i)
                    values(i) = values(j)
                    values(j) = temp
                    i = i + 1
                    j = j - 1
                end if

                if (i > j) exit
            end do

            if (kth_index <= j) then
                right = j
            else if (kth_index >= i) then
                left = i
            else
                return
            end if
        end do

        call insertion_sort_real(values, left, right)
    end subroutine select_kth_real


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE MANN-KENDALL S STATISTIC FOR ONE CLEAN TIME SERIES.
    !
    ! INPUTS
    !    Y(:) - INPUT SERIES WITH NO NaN VALUES
    !
    ! OUTPUT
    !    S_VALUE - SUM OF PAIRWISE SIGN COMPARISONS
    subroutine compute_s_value(y, s_value)
        real(real64), intent(in) :: y(:)
        real(real64), intent(out) :: s_value

        integer :: i, j, n

        n = size(y)
        s_value = 0.0_real64

        do i = 1, n - 1
            do j = i + 1, n
                if (y(j) > y(i)) then
                    s_value = s_value + 1.0_real64
                else if (y(j) < y(i)) then
                    s_value = s_value - 1.0_real64
                end if
            end do
        end do
    end subroutine compute_s_value


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE TIE-CORRECTED VARIANCE OF THE MANN-KENDALL S
    !    STATISTIC FOR ONE CLEAN TIME SERIES.
    !
    ! INPUTS
    !    Y(:) - INPUT SERIES WITH NO NaN VALUES
    !
    ! INPUT / OUTPUT
    !    SORTED_WORK(:) - CALLER-SUPPLIED TEMPORARY BUFFER OF LENGTH AT
    !                     LEAST SIZE(Y); FILLED WITH A SORTED COPY OF `Y`
    !
    ! OUTPUT
    !    VAR_S - TIE-CORRECTED VARIANCE OF S
    subroutine compute_base_variance(y, sorted_work, var_s)
        real(real64), intent(in) :: y(:)
        real(real64), intent(inout) :: sorted_work(:)
        real(real64), intent(out) :: var_s

        integer :: i, n, run_length
        real(real64) :: n_real, tie_term, run_real

        n = size(y)
        if (n <= 1) then
            var_s = 0.0_real64
            return
        end if

        sorted_work(:n) = y
        call sort_real_inplace(sorted_work(:n))

        tie_term = 0.0_real64
        run_length = 1

        do i = 2, n
            if (sorted_work(i) == sorted_work(i - 1)) then
                run_length = run_length + 1
            else
                if (run_length > 1) then
                    run_real = real(run_length, real64)
                    tie_term = tie_term + run_real * (run_real - 1.0_real64) * &
                        (2.0_real64 * run_real + 5.0_real64)
                end if
                run_length = 1
            end if
        end do

        if (run_length > 1) then
            run_real = real(run_length, real64)
            tie_term = tie_term + run_real * (run_real - 1.0_real64) * &
                (2.0_real64 * run_real + 5.0_real64)
        end if

        n_real = real(n, real64)
        var_s = (n_real * (n_real - 1.0_real64) * (2.0_real64 * n_real + 5.0_real64) - &
            tie_term) / 18.0_real64
    end subroutine compute_base_variance


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE EXACT SEN / THEIL-SEN SLOPE MEDIAN FOR ONE CLEAN
    !    TIME SERIES.
    !
    ! INPUTS
    !    Y(:) - INPUT SERIES WITH NO NaN VALUES AND UNIT TIME SPACING
    !
    ! INPUT / OUTPUT
    !    SLOPE_WORK(:) - CALLER-SUPPLIED BUFFER WITH LENGTH AT LEAST
    !                    N * (N - 1) / 2 SO ALL PAIRWISE SLOPES CAN BE
    !                    MATERIALIZED BEFORE THE MEDIAN IS SELECTED
    !
    ! OUTPUT
    !    SLOPE - EXACT MEDIAN OF ALL PAIRWISE SLOPES; NaN WHEN SIZE(Y) < 2
    subroutine compute_sen_slope(y, slope, slope_work)
        real(real64), intent(in) :: y(:)
        real(real64), intent(out) :: slope
        real(real64), intent(inout) :: slope_work(:)

        integer :: i, j, n, nslopes, mid_index, slope_index

        n = size(y)
        if (n < 2) then
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        nslopes = n * (n - 1) / 2
        slope_index = 0

        do i = 1, n - 1
            do j = i + 1, n
                slope_index = slope_index + 1
                slope_work(slope_index) = (y(j) - y(i)) / real(j - i, real64)
            end do
        end do

        mid_index = nslopes / 2
        if (mod(nslopes, 2) == 0) then
            call select_kth_real(slope_work(:nslopes), mid_index)
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            slope = 0.5_real64 * (slope_work(mid_index) + slope_work(mid_index + 1))
        else
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            slope = slope_work(mid_index + 1)
        end if
    end subroutine compute_sen_slope


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE YUE-WANG MODIFIED VARIANCE AFTER DETRENDING BY A
    !    CALLER-SUPPLIED SEN SLOPE.
    !
    ! INPUTS
    !    Y(:)    - INPUT SERIES WITH NO NaN VALUES
    !    SLOPE   - PRECOMPUTED SEN SLOPE USED FOR DETRENDING
    !
    ! INPUT / OUTPUT
    !    SORTED_WORK(:)    - TEMPORARY BUFFER OF LENGTH AT LEAST SIZE(Y)
    !    CENTERED_WORK(:)  - TEMPORARY BUFFER OF LENGTH AT LEAST SIZE(Y)
    !    DETRENDED_WORK(:) - TEMPORARY BUFFER OF LENGTH AT LEAST SIZE(Y)
    !
    ! OUTPUT
    !    VAR_S - MODIFIED VARIANCE; NaN WHEN THE DETRENDED SERIES HAS
    !            ZERO VARIANCE AND THE AUTOCORRELATION TERM IS UNDEFINED
    subroutine compute_modified_variance_with_slope( &
        y, sorted_work, centered_work, detrended_work, slope, var_s &
    )
        real(real64), intent(in) :: y(:)
        real(real64), intent(inout) :: sorted_work(:)
        real(real64), intent(inout) :: centered_work(:), detrended_work(:)
        real(real64), intent(in) :: slope
        real(real64), intent(out) :: var_s

        integer :: i, lag, n
        real(real64) :: base_var, denom, mean_value, n_star, numer

        n = size(y)
        call compute_base_variance(y, sorted_work, base_var)

        do i = 1, n
            detrended_work(i) = y(i) - real(i, real64) * slope
        end do

        mean_value = sum(detrended_work(:n)) / real(n, real64)
        centered_work(:n) = detrended_work(:n) - mean_value
        denom = sum(centered_work(:n) * centered_work(:n))

        if (denom == 0.0_real64) then
            var_s = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        n_star = 1.0_real64
        do lag = 1, n - 1
            numer = sum(centered_work(1:n - lag) * centered_work(1 + lag:n))
            n_star = n_star + 2.0_real64 * &
                (1.0_real64 - real(lag, real64) / real(n, real64)) * &
                (numer / denom)
        end do

        var_s = base_var * n_star
    end subroutine compute_modified_variance_with_slope


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE YUE-WANG MODIFIED VARIANCE WHEN ONLY THE INPUT
    !    SERIES IS KNOWN, BY FIRST COMPUTING THE SEN SLOPE INTERNALLY.
    !
    ! INPUTS
    !    Y(:) - INPUT SERIES WITH NO NaN VALUES
    !
    ! INPUT / OUTPUT
    !    SORTED_WORK(:)    - TEMPORARY BUFFER OF LENGTH AT LEAST SIZE(Y)
    !    SLOPE_WORK(:)     - TEMPORARY BUFFER OF LENGTH AT LEAST
    !                        N * (N - 1) / 2
    !    CENTERED_WORK(:)  - TEMPORARY BUFFER OF LENGTH AT LEAST SIZE(Y)
    !    DETRENDED_WORK(:) - TEMPORARY BUFFER OF LENGTH AT LEAST SIZE(Y)
    !
    ! OUTPUT
    !    VAR_S - MODIFIED VARIANCE OF THE INPUT SERIES
    subroutine compute_modified_variance(y, sorted_work, slope_work, centered_work, detrended_work, var_s)
        real(real64), intent(in) :: y(:)
        real(real64), intent(inout) :: sorted_work(:), slope_work(:)
        real(real64), intent(inout) :: centered_work(:), detrended_work(:)
        real(real64), intent(out) :: var_s

        real(real64) :: slope

        call compute_sen_slope(y, slope, slope_work)
        call compute_modified_variance_with_slope( &
            y, sorted_work, centered_work, detrended_work, slope, var_s &
        )
    end subroutine compute_modified_variance

end module mk_kernels_core


! QUICK REFERENCE
! PURPOSE
!    BATCH ENTRY POINT FOR THE MANN-KENDALL S STATISTIC AND ITS
!    VARIANCE OVER MULTIPLE CLEAN SERIES.
!
! EXPECTED INPUT SHAPES
!    DATA(NTIME, NSERIES) - EACH COLUMN IS ONE CLEAN TIME SERIES
!
! FLAGS
!    MODIFIED - 0: USE THE ORIGINAL MANN-KENDALL VARIANCE
!               NONZERO: USE THE YUE-WANG MODIFIED VARIANCE
!
! OUTPUT
!    S_VALUES(NSERIES)   - ONE S STATISTIC PER COLUMN
!    VAR_VALUES(NSERIES) - ONE VARIANCE VALUE PER COLUMN
subroutine mk_score_var_batch_clean(data, s_values, var_values, modified, ntime, nseries)
    use mk_kernels_core, only : compute_base_variance, compute_modified_variance, compute_s_value, real64
    implicit none

    integer, intent(in) :: modified, ntime, nseries
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: s_values(nseries), var_values(nseries)

    real(real64) :: centered_work(ntime), detrended_work(ntime), sorted_work(ntime)
    real(real64) :: slope_work(max(1, ntime * (ntime - 1) / 2))
    integer :: col

    do col = 1, nseries
        call compute_s_value(data(:, col), s_values(col))

        if (modified /= 0) then
            call compute_modified_variance( &
                data(:, col), sorted_work, slope_work, centered_work, detrended_work, var_values(col) &
            )
        else
            call compute_base_variance(data(:, col), sorted_work, var_values(col))
        end if
    end do
end subroutine mk_score_var_batch_clean


! QUICK REFERENCE
! PURPOSE
!    BATCH ENTRY POINT FOR THE EXACT SEN / THEIL-SEN SLOPE.
!
! EXPECTED INPUT SHAPES
!    DATA(NTIME, NSERIES) - EACH COLUMN IS ONE CLEAN TIME SERIES WITH
!                           UNIT TIME SPACING
!
! OUTPUT
!    SLOPES(NSERIES) - ONE EXACT SEN SLOPE PER COLUMN
subroutine sen_slope_batch_clean(data, slopes, ntime, nseries)
    use mk_kernels_core, only : compute_sen_slope, real64
    implicit none

    integer, intent(in) :: ntime, nseries
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: slopes(nseries)

    real(real64) :: slope_work(max(1, ntime * (ntime - 1) / 2))
    integer :: col

    do col = 1, nseries
        call compute_sen_slope(data(:, col), slopes(col), slope_work)
    end do
end subroutine sen_slope_batch_clean


! QUICK REFERENCE
! PURPOSE
!    BATCH ENTRY POINT THAT RETURNS S, VARIANCE, AND SEN SLOPE IN ONE
!    PASS SO THE PYTHON LAYER CAN AVOID REDUNDANT SLOPE WORK.
!
! EXPECTED INPUT SHAPES
!    DATA(NTIME, NSERIES) - EACH COLUMN IS ONE CLEAN TIME SERIES
!
! FLAGS
!    MODIFIED - 0: USE THE ORIGINAL MANN-KENDALL VARIANCE
!               NONZERO: USE THE YUE-WANG MODIFIED VARIANCE
!
! OUTPUT
!    S_VALUES(NSERIES)   - ONE S STATISTIC PER COLUMN
!    VAR_VALUES(NSERIES) - ONE VARIANCE VALUE PER COLUMN
!    SLOPES(NSERIES)     - ONE EXACT SEN SLOPE PER COLUMN
subroutine mk_score_var_sen_batch_clean( &
    data, s_values, var_values, slopes, modified, ntime, nseries &
)
    use mk_kernels_core, only : &
        compute_base_variance, compute_modified_variance_with_slope, compute_s_value, &
        compute_sen_slope, real64
    implicit none

    integer, intent(in) :: modified, ntime, nseries
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: s_values(nseries), var_values(nseries), slopes(nseries)

    real(real64) :: centered_work(ntime), detrended_work(ntime), sorted_work(ntime)
    real(real64) :: slope_work(max(1, ntime * (ntime - 1) / 2))
    integer :: col

    do col = 1, nseries
        call compute_s_value(data(:, col), s_values(col))
        call compute_sen_slope(data(:, col), slopes(col), slope_work)

        if (modified /= 0) then
            call compute_modified_variance_with_slope( &
                data(:, col), sorted_work, centered_work, detrended_work, slopes(col), &
                var_values(col) &
            )
        else
            call compute_base_variance(data(:, col), sorted_work, var_values(col))
        end if
    end do
end subroutine mk_score_var_sen_batch_clean
