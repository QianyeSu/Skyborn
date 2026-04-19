! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-19
! File    : mann_kendall_core.f90
! Purpose : Provide the private optional Fortran backend for the clean-series
!           Mann-Kendall hot path in `skyborn.calc.mann_kendall`.
! Notes   : Public batch entry points operate on `data(ntime, nseries)` with
!           no NaN values. Missing-value filtering and fallback dispatch stay
!           on the Python side so this backend can focus on dense numeric work.
! =============================================================================
!
module mann_kendall_core_mod
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    public :: real64
    public :: compute_base_variance
    public :: compute_modified_variance
    public :: compute_modified_variance_with_slope
    public :: compute_s_value
    public :: compute_grouped_slope_count
    public :: compute_correlated_grouped_stats
    public :: compute_s_value_and_sen_slope
    public :: compute_sen_slope
    public :: compute_grouped_sen_slope_with_inv_lag
    public :: compute_partial_stats
    public :: compute_s_value_and_sen_slope_with_inv_lag
    public :: compute_sen_slope_with_inv_lag

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
    !    STATISTIC FROM AN ALREADY SORTED COPY OF THE INPUT SERIES.
    !
    ! INPUTS
    !    N              - LOGICAL LENGTH OF THE ACTIVE SERIES
    !    SORTED_WORK(:) - ASCENDING SORTED COPY OF THE INPUT SERIES IN
    !                     SORTED_WORK(1:N)
    !
    ! OUTPUT
    !    VAR_S - TIE-CORRECTED VARIANCE OF S
    subroutine compute_base_variance_from_sorted(sorted_work, n, var_s)
        real(real64), intent(in) :: sorted_work(:)
        integer, intent(in) :: n
        real(real64), intent(out) :: var_s

        integer :: i, run_length
        real(real64) :: n_real, tie_term, run_real

        if (n <= 1) then
            var_s = 0.0_real64
            return
        end if

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
    end subroutine compute_base_variance_from_sorted


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

        integer :: n

        n = size(y)
        sorted_work(:n) = y
        call sort_real_inplace(sorted_work(:n))
        call compute_base_variance_from_sorted(sorted_work, n, var_s)
    end subroutine compute_base_variance


    ! QUICK REFERENCE
    ! PURPOSE
    !    RETURN THE TOTAL NUMBER OF WITHIN-GROUP PAIRWISE SLOPES WHEN
    !    A CLEAN FLAT SERIES IS SPLIT INTO `PERIOD` INTERLEAVED GROUPS.
    !
    ! INPUTS
    !    N      - TOTAL FLAT SERIES LENGTH
    !    PERIOD - NUMBER OF INTERLEAVED GROUPS
    !
    ! OUTPUT
    !    N_SLOPES - TOTAL COUNT OF WITHIN-GROUP PAIRWISE SLOPES
    pure integer function compute_grouped_slope_count(n, period) result(n_slopes)
        integer, intent(in) :: n, period

        integer :: group_count, group_index

        n_slopes = 0
        if (n < 2 .or. period <= 0) return

        do group_index = 1, min(period, n)
            group_count = 1 + (n - group_index) / period
            if (group_count > 1) then
                n_slopes = n_slopes + group_count * (group_count - 1) / 2
            end if
        end do
    end function compute_grouped_slope_count


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

        real(real64) :: inv_lag(max(1, size(y) - 1))
        integer :: lag, n

        n = size(y)
        if (n < 2) then
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        do lag = 1, n - 1
            inv_lag(lag) = 1.0_real64 / real(lag, real64)
        end do

        call compute_sen_slope_with_inv_lag(y, inv_lag, slope, slope_work)
    end subroutine compute_sen_slope


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE EXACT SEN / THEIL-SEN SLOPE MEDIAN FOR ONE CLEAN
    !    TIME SERIES USING A PRECOMPUTED RECIPROCAL-LAG VECTOR.
    !
    ! INPUTS
    !    Y(:)        - INPUT SERIES WITH NO NaN VALUES AND UNIT TIME SPACING
    !    INV_LAG(:)  - PRECOMPUTED VALUES 1 / LAG FOR LAG = 1..SIZE(Y)-1
    !
    ! INPUT / OUTPUT
    !    SLOPE_WORK(:) - CALLER-SUPPLIED BUFFER WITH LENGTH AT LEAST
    !                    N * (N - 1) / 2
    !
    ! OUTPUT
    !    SLOPE - EXACT MEDIAN OF ALL PAIRWISE SLOPES; NaN WHEN SIZE(Y) < 2
    !
    ! NOTES
    !    FOR EVEN COUNTS WE ONLY PARTITION ONCE AROUND THE UPPER MEDIAN.
    !    AFTER THAT, THE LOWER MEDIAN IS SIMPLY THE MAXIMUM VALUE IN THE
    !    FIRST HALF OF THE PARTITIONED BUFFER, WHICH AVOIDS A SECOND
    !    QUICKSELECT PASS WITHOUT CHANGING THE EXACT RESULT.
    subroutine compute_sen_slope_with_inv_lag(y, inv_lag, slope, slope_work)
        real(real64), intent(in) :: y(:), inv_lag(:)
        real(real64), intent(out) :: slope
        real(real64), intent(inout) :: slope_work(:)

        integer :: i, lag, n, nslopes, mid_index, slope_index
        real(real64) :: factor, lower_median, upper_median

        n = size(y)
        if (n < 2) then
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        nslopes = n * (n - 1) / 2
        slope_index = 0

        do lag = 1, n - 1
            factor = inv_lag(lag)
            do i = 1, n - lag
                slope_index = slope_index + 1
                slope_work(slope_index) = (y(i + lag) - y(i)) * factor
            end do
        end do

        mid_index = nslopes / 2
        if (mod(nslopes, 2) == 0) then
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            upper_median = slope_work(mid_index + 1)
            lower_median = maxval(slope_work(:mid_index))
            slope = 0.5_real64 * (lower_median + upper_median)
        else
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            slope = slope_work(mid_index + 1)
        end if
    end subroutine compute_sen_slope_with_inv_lag


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE EXACT GROUPED SEN / THEIL-SEN SLOPE FOR A CLEAN
    !    FLAT SERIES REPRESENTING `PERIOD` INTERLEAVED GROUPS.
    !
    ! INPUTS
    !    Y(:)        - CLEAN FLAT INPUT SERIES
    !    PERIOD      - NUMBER OF INTERLEAVED GROUPS
    !    INV_LAG(:)  - PRECOMPUTED 1 / LAG VALUES FOR GROUP-INTERNAL
    !                  CYCLE LAGS
    !
    ! INPUT / OUTPUT
    !    SLOPE_WORK(:) - CALLER-SUPPLIED BUFFER WITH LENGTH AT LEAST
    !                    `compute_grouped_slope_count(size(y), period)`
    !
    ! OUTPUT
    !    SLOPE - EXACT GROUPED MEDIAN SLOPE; NaN WHEN FEWER THAN TWO
    !            WITHIN-GROUP PAIRS EXIST
    subroutine compute_grouped_sen_slope_with_inv_lag(y, period, inv_lag, slope, slope_work)
        real(real64), intent(in) :: y(:), inv_lag(:)
        integer, intent(in) :: period
        real(real64), intent(out) :: slope
        real(real64), intent(inout) :: slope_work(:)

        integer :: group_count, group_index, i, lag, mid_index, n, nslopes, slope_index
        real(real64) :: factor, lower_median, upper_median

        n = size(y)
        if (n < 2 .or. period <= 0) then
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        nslopes = compute_grouped_slope_count(n, period)
        if (nslopes <= 0) then
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        slope_index = 0
        do group_index = 1, min(period, n)
            group_count = 1 + (n - group_index) / period
            do lag = 1, group_count - 1
                factor = inv_lag(lag)
                do i = 1, group_count - lag
                    slope_index = slope_index + 1
                    slope_work(slope_index) = ( &
                        y(group_index + (i - 1 + lag) * period) - &
                        y(group_index + (i - 1) * period) &
                    ) * factor
                end do
            end do
        end do

        mid_index = nslopes / 2
        if (mod(nslopes, 2) == 0) then
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            upper_median = slope_work(mid_index + 1)
            lower_median = maxval(slope_work(:mid_index))
            slope = 0.5_real64 * (lower_median + upper_median)
        else
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            slope = slope_work(mid_index + 1)
        end if
    end subroutine compute_grouped_sen_slope_with_inv_lag


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE AVERAGE RANKS FOR ONE GROUP USING THE SAME R DEFINITION
    !    AS `pymannkendall`'s INTERNAL CORRELATED MULTIVARIATE HELPERS.
    !
    ! INPUTS
    !    VALUES(:) - GROUP SERIES OF LENGTH N
    !    N         - ACTIVE LENGTH
    !
    ! OUTPUT
    !    RANKS(:)  - AVERAGE RANKS OVER THE ACTIVE WINDOW 1:N
    subroutine compute_rank_vector(values, n, ranks)
        real(real64), intent(in) :: values(:)
        integer, intent(in) :: n
        real(real64), intent(out) :: ranks(:)

        integer :: i, j, s

        if (n <= 0) return

        do j = 1, n
            s = 0
            do i = 1, n
                if (values(j) > values(i)) then
                    s = s + 1
                else if (values(j) < values(i)) then
                    s = s - 1
                end if
            end do
            ranks(j) = 0.5_real64 * real(n + 1 + s, real64)
        end do
    end subroutine compute_rank_vector


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE K STATISTIC USED BY THE CORRELATED MULTIVARIATE
    !    MANN-KENDALL VARIANCE FORMULA.
    !
    ! INPUTS
    !    LEFT(:), RIGHT(:) - GROUP SERIES OF COMMON LENGTH N
    !    N                 - ACTIVE LENGTH
    !
    ! OUTPUT
    !    K_VALUE - PAIRWISE SIGN-OF-PRODUCT SUM
    subroutine compute_k_stat(left, right, n, k_value)
        real(real64), intent(in) :: left(:), right(:)
        integer, intent(in) :: n
        real(real64), intent(out) :: k_value

        integer :: i, j
        real(real64) :: product

        k_value = 0.0_real64
        if (n <= 1) return

        do i = 1, n - 1
            do j = i + 1, n
                product = (left(j) - left(i)) * (right(j) - right(i))
                if (product > 0.0_real64) then
                    k_value = k_value + 1.0_real64
                else if (product < 0.0_real64) then
                    k_value = k_value - 1.0_real64
                end if
            end do
        end do
    end subroutine compute_k_stat


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE CORRELATED GROUPED MANN-KENDALL S STATISTIC,
    !    VARIANCE, AND TAU DENOMINATOR FOR ONE CLEAN FLAT INTERLEAVED
    !    SERIES WITH A FIXED GROUP COUNT.
    !
    ! INPUTS
    !    Y(:)     - CLEAN FLAT INPUT SERIES
    !    PERIOD   - NUMBER OF INTERLEAVED GROUPS
    !
    ! INPUT / OUTPUT
    !    GROUP_WORK(:,:) - WORK BUFFER OF SHAPE AT LEAST
    !                      (SIZE(Y) / PERIOD, PERIOD)
    !    RANK_WORK(:,:)  - WORK BUFFER OF SAME SHAPE AS GROUP_WORK
    !
    ! OUTPUT
    !    S_VALUE - SUM OF GROUP-WISE MANN-KENDALL SCORES
    !    VAR_S   - CORRELATED MULTIVARIATE VARIANCE
    !    DENOM   - TAU NORMALIZATION DENOMINATOR
    subroutine compute_correlated_grouped_stats(y, period, s_value, var_s, denom, group_work, rank_work)
        real(real64), intent(in) :: y(:)
        integer, intent(in) :: period
        real(real64), intent(out) :: s_value, var_s, denom
        real(real64), intent(inout) :: group_work(:,:), rank_work(:,:)

        integer :: cycle_index, group_i, group_j, n, ncycles
        real(real64) :: gamma_value, k_value, rank_sum
        real(real64) :: n_real

        n = size(y)
        if (period <= 0) then
            s_value = 0.0_real64
            var_s = 0.0_real64
            denom = 0.0_real64
            return
        end if

        ncycles = n / period
        if (ncycles <= 0) then
            s_value = 0.0_real64
            var_s = 0.0_real64
            denom = 0.0_real64
            return
        end if

        do group_i = 1, period
            do cycle_index = 1, ncycles
                group_work(cycle_index, group_i) = y(group_i + (cycle_index - 1) * period)
            end do
        end do

        s_value = 0.0_real64
        denom = 0.0_real64
        if (ncycles > 1) then
            do group_i = 1, period
                call compute_s_value(group_work(:ncycles, group_i), gamma_value)
                s_value = s_value + gamma_value
                denom = denom + 0.5_real64 * real(ncycles * (ncycles - 1), real64)
                call compute_rank_vector(group_work(:ncycles, group_i), ncycles, rank_work(:ncycles, group_i))
            end do
        else
            do group_i = 1, period
                call compute_rank_vector(group_work(:ncycles, group_i), ncycles, rank_work(:ncycles, group_i))
            end do
        end if

        var_s = 0.0_real64
        n_real = real(ncycles, real64)

        do group_i = 1, period
            call compute_k_stat(group_work(:ncycles, group_i), group_work(:ncycles, group_i), ncycles, k_value)
            rank_sum = sum(rank_work(:ncycles, group_i) * rank_work(:ncycles, group_i))
            gamma_value = (k_value + 4.0_real64 * rank_sum - n_real * (n_real + 1.0_real64) ** 2) / 3.0_real64
            var_s = var_s + gamma_value

            do group_j = 1, group_i - 1
                call compute_k_stat(group_work(:ncycles, group_i), group_work(:ncycles, group_j), ncycles, k_value)
                rank_sum = sum(rank_work(:ncycles, group_i) * rank_work(:ncycles, group_j))
                gamma_value = (k_value + 4.0_real64 * rank_sum - n_real * (n_real + 1.0_real64) ** 2) / 3.0_real64
                var_s = var_s + 2.0_real64 * gamma_value
            end do
        end do
    end subroutine compute_correlated_grouped_stats


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE PARTIAL MANN-KENDALL STATISTICS FOR ONE CLEAN
    !    RESPONSE / COVARIATE PAIR.
    !
    ! INPUTS
    !    X(:) - CLEAN RESPONSE SERIES
    !    Y(:) - CLEAN COVARIATE SERIES
    !
    ! INPUT / OUTPUT
    !    RANK_X(:), RANK_Y(:) - WORK BUFFERS OF SIZE AT LEAST SIZE(X)
    !
    ! OUTPUT
    !    S_VALUE - PARTIAL MANN-KENDALL S STATISTIC
    !    VAR_S   - PARTIAL MANN-KENDALL VARIANCE
    !    TAU     - RESPONSE-SERIES KENDALL TAU NORMALIZATION
    subroutine compute_partial_stats(x, y, s_value, var_s, tau, rank_x, rank_y)
        real(real64), intent(in) :: x(:), y(:)
        real(real64), intent(out) :: s_value, var_s, tau
        real(real64), intent(inout) :: rank_x(:), rank_y(:)

        integer :: i, j, n
        real(real64) :: x_score, y_score, k_value, sigma, rho
        real(real64) :: delta_x, delta_y, sign_x, sign_y
        real(real64) :: var_no_ties, n_real

        n = size(x)
        if (size(y) /= n .or. n < 3) then
            s_value = 0.0_real64
            var_s = 0.0_real64
            tau = 0.0_real64
            return
        end if

        n_real = real(n, real64)
        var_no_ties = n_real * real(n - 1, real64) * real(2 * n + 5, real64) / 18.0_real64

        x_score = 0.0_real64
        y_score = 0.0_real64
        k_value = 0.0_real64
        rank_x(:n) = 0.0_real64
        rank_y(:n) = 0.0_real64

        do i = 1, n - 1
            do j = i + 1, n
                delta_x = x(j) - x(i)
                if (delta_x > 0.0_real64) then
                    sign_x = 1.0_real64
                else if (delta_x < 0.0_real64) then
                    sign_x = -1.0_real64
                else
                    sign_x = 0.0_real64
                end if
                x_score = x_score + sign_x
                rank_x(i) = rank_x(i) - sign_x
                rank_x(j) = rank_x(j) + sign_x

                delta_y = y(j) - y(i)
                if (delta_y > 0.0_real64) then
                    sign_y = 1.0_real64
                else if (delta_y < 0.0_real64) then
                    sign_y = -1.0_real64
                else
                    sign_y = 0.0_real64
                end if
                y_score = y_score + sign_y
                rank_y(i) = rank_y(i) - sign_y
                rank_y(j) = rank_y(j) + sign_y

                k_value = k_value + sign_x * sign_y
            end do
        end do

        rank_x(:n) = (n_real + 1.0_real64 + rank_x(:n)) / 2.0_real64
        rank_y(:n) = (n_real + 1.0_real64 + rank_y(:n)) / 2.0_real64

        sigma = ( &
            k_value + 4.0_real64 * sum(rank_x(:n) * rank_y(:n)) - n_real * (n_real + 1.0_real64) ** 2 &
        ) / 3.0_real64
        rho = sigma / var_no_ties

        s_value = x_score - rho * y_score
        var_s = (1.0_real64 - rho * rho) * var_no_ties
        tau = x_score / (0.5_real64 * n_real * real(n - 1, real64))
    end subroutine compute_partial_stats


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE MANN-KENDALL S STATISTIC AND THE EXACT SEN /
    !    THEIL-SEN SLOPE USING ONE SHARED PAIRWISE SCAN.
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
    !    S_VALUE - SUM OF PAIRWISE SIGN COMPARISONS
    !    SLOPE   - EXACT MEDIAN OF ALL PAIRWISE SLOPES; NaN WHEN SIZE(Y) < 2
    subroutine compute_s_value_and_sen_slope(y, s_value, slope, slope_work)
        real(real64), intent(in) :: y(:)
        real(real64), intent(out) :: s_value, slope
        real(real64), intent(inout) :: slope_work(:)

        real(real64) :: inv_lag(max(1, size(y) - 1))
        integer :: lag, n

        n = size(y)
        if (n < 2) then
            s_value = 0.0_real64
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        do lag = 1, n - 1
            inv_lag(lag) = 1.0_real64 / real(lag, real64)
        end do

        call compute_s_value_and_sen_slope_with_inv_lag(y, inv_lag, s_value, slope, slope_work)
    end subroutine compute_s_value_and_sen_slope


    ! QUICK REFERENCE
    ! PURPOSE
    !    COMPUTE THE MANN-KENDALL S STATISTIC AND THE EXACT SEN /
    !    THEIL-SEN SLOPE USING ONE SHARED PAIRWISE SCAN AND A
    !    PRECOMPUTED RECIPROCAL-LAG VECTOR.
    !
    ! INPUTS
    !    Y(:)        - INPUT SERIES WITH NO NaN VALUES AND UNIT TIME SPACING
    !    INV_LAG(:)  - PRECOMPUTED VALUES 1 / LAG FOR LAG = 1..SIZE(Y)-1
    !
    ! INPUT / OUTPUT
    !    SLOPE_WORK(:) - CALLER-SUPPLIED BUFFER WITH LENGTH AT LEAST
    !                    N * (N - 1) / 2
    !
    ! OUTPUT
    !    S_VALUE - SUM OF PAIRWISE SIGN COMPARISONS
    !    SLOPE   - EXACT MEDIAN OF ALL PAIRWISE SLOPES; NaN WHEN SIZE(Y) < 2
    !
    ! NOTES
    !    THE EVEN-SLOPE MEDIAN CASE REUSES THE SAME SINGLE-PARTITION
    !    STRATEGY AS `compute_sen_slope_with_inv_lag(...)` SO THE BATCH
    !    KERNEL DOES NOT PAY FOR TWO FULL QUICKSELECT PASSES.
    subroutine compute_s_value_and_sen_slope_with_inv_lag(y, inv_lag, s_value, slope, slope_work)
        real(real64), intent(in) :: y(:), inv_lag(:)
        real(real64), intent(out) :: s_value, slope
        real(real64), intent(inout) :: slope_work(:)

        integer :: i, lag, n, nslopes, mid_index, slope_index
        real(real64) :: dy, factor, lower_median, upper_median

        n = size(y)
        if (n < 2) then
            s_value = 0.0_real64
            slope = ieee_value(0.0_real64, ieee_quiet_nan)
            return
        end if

        nslopes = n * (n - 1) / 2
        slope_index = 0
        s_value = 0.0_real64

        do lag = 1, n - 1
            factor = inv_lag(lag)
            do i = 1, n - lag
                dy = y(i + lag) - y(i)
                slope_index = slope_index + 1
                slope_work(slope_index) = dy * factor

                if (dy > 0.0_real64) then
                    s_value = s_value + 1.0_real64
                else if (dy < 0.0_real64) then
                    s_value = s_value - 1.0_real64
                end if
            end do
        end do

        mid_index = nslopes / 2
        if (mod(nslopes, 2) == 0) then
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            upper_median = slope_work(mid_index + 1)
            lower_median = maxval(slope_work(:mid_index))
            slope = 0.5_real64 * (lower_median + upper_median)
        else
            call select_kth_real(slope_work(:nslopes), mid_index + 1)
            slope = slope_work(mid_index + 1)
        end if
    end subroutine compute_s_value_and_sen_slope_with_inv_lag


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

        do i = 1, n
            sorted_work(i) = y(i)
            detrended_work(i) = y(i) - real(i, real64) * slope
        end do
        call sort_real_inplace(sorted_work(:n))
        call compute_base_variance_from_sorted(sorted_work, n, base_var)

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

end module mann_kendall_core_mod


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
subroutine mk_score_var_batch(data, s_values, var_values, modified, ntime, nseries)
    use mann_kendall_core_mod, only : compute_base_variance, compute_modified_variance, compute_s_value, real64
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
end subroutine mk_score_var_batch


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
subroutine sen_slope_batch(data, slopes, ntime, nseries)
    use mann_kendall_core_mod, only : compute_sen_slope_with_inv_lag, real64
    implicit none

    integer, intent(in) :: ntime, nseries
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: slopes(nseries)

    real(real64) :: inv_lag(max(1, ntime - 1))
    real(real64) :: slope_work(max(1, ntime * (ntime - 1) / 2))
    integer :: col, lag

    do lag = 1, ntime - 1
        inv_lag(lag) = 1.0_real64 / real(lag, real64)
    end do

    do col = 1, nseries
        call compute_sen_slope_with_inv_lag(data(:, col), inv_lag, slopes(col), slope_work)
    end do
end subroutine sen_slope_batch


! QUICK REFERENCE
! PURPOSE
!    BATCH ENTRY POINT FOR THE EXACT GROUPED SEN / THEIL-SEN SLOPE USED
!    BY SEASONAL AND OTHER INTERLEAVED GROUPED MANN-KENDALL FAMILIES.
!
! EXPECTED INPUT SHAPES
!    DATA(NTIME, NSERIES) - EACH COLUMN IS ONE CLEAN FLAT SERIES
!
! INPUTS
!    PERIOD - NUMBER OF INTERLEAVED GROUPS WITHIN EACH COLUMN
!
! OUTPUT
!    SLOPES(NSERIES) - ONE EXACT GROUPED SEN SLOPE PER COLUMN
subroutine grouped_sen_slope_batch(data, slopes, period, ntime, nseries)
    use mann_kendall_core_mod, only : compute_grouped_sen_slope_with_inv_lag, &
        compute_grouped_slope_count, real64
    implicit none

    integer, intent(in) :: ntime, nseries, period
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: slopes(nseries)

    integer :: col, lag, max_group_size
    real(real64), allocatable :: inv_lag(:), slope_work(:)

    max_group_size = 1 + max(0, ntime - 1) / max(1, period)
    allocate(inv_lag(max(1, max_group_size - 1)))
    if (max_group_size > 1) then
        do lag = 1, max_group_size - 1
            inv_lag(lag) = 1.0_real64 / real(lag, real64)
        end do
    end if

    allocate(slope_work(max(1, compute_grouped_slope_count(ntime, period))))

    do col = 1, nseries
        call compute_grouped_sen_slope_with_inv_lag( &
            data(:, col), period, inv_lag, slopes(col), slope_work &
        )
    end do
end subroutine grouped_sen_slope_batch


! QUICK REFERENCE
! PURPOSE
!    BATCH ENTRY POINT FOR CORRELATED GROUPED MANN-KENDALL STATISTICS.
!
! EXPECTED INPUT SHAPES
!    DATA(NTIME, NSERIES) - EACH COLUMN IS ONE CLEAN FLAT INTERLEAVED SERIES
!
! INPUTS
!    PERIOD - NUMBER OF INTERLEAVED GROUPS WITHIN EACH COLUMN
!
! OUTPUT
!    S_VALUES(NSERIES)   - ONE CORRELATED GROUPED S VALUE PER COLUMN
!    VAR_VALUES(NSERIES) - ONE CORRELATED GROUPED VARIANCE PER COLUMN
!    DENOM               - COMMON TAU DENOMINATOR FOR THE CLEAN BATCH
subroutine grouped_correlated_stats_batch(data, s_values, var_values, denom, period, ntime, nseries)
    use mann_kendall_core_mod, only : compute_correlated_grouped_stats, real64
    implicit none

    integer, intent(in) :: ntime, nseries, period
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: s_values(nseries), var_values(nseries), denom

    integer :: col, ncycles
    real(real64), allocatable :: group_work(:,:), rank_work(:,:)
    real(real64) :: denom_value

    ncycles = ntime / max(1, period)
    allocate(group_work(max(1, ncycles), max(1, period)))
    allocate(rank_work(max(1, ncycles), max(1, period)))

    denom = 0.0_real64
    do col = 1, nseries
        call compute_correlated_grouped_stats( &
            data(:, col), period, s_values(col), var_values(col), denom_value, group_work, rank_work &
        )
        if (col == 1) denom = denom_value
    end do
end subroutine grouped_correlated_stats_batch


! QUICK REFERENCE
! PURPOSE
!    BATCH ENTRY POINT FOR PARTIAL MANN-KENDALL STATISTICS.
!
! EXPECTED INPUT SHAPES
!    RESPONSE(NTIME, NSERIES)  - EACH COLUMN IS ONE CLEAN RESPONSE SERIES
!    COVARIATE(NTIME, NSERIES) - MATCHING CLEAN COVARIATE SERIES
!
! OUTPUT
!    S_VALUES(NSERIES)   - ONE PARTIAL MK S VALUE PER COLUMN
!    VAR_VALUES(NSERIES) - ONE PARTIAL MK VARIANCE PER COLUMN
!    TAU_VALUES(NSERIES) - RESPONSE-SERIES KENDALL TAU PER COLUMN
subroutine partial_stats_batch(response, covariate, s_values, var_values, tau_values, ntime, nseries)
    use mann_kendall_core_mod, only : compute_partial_stats, real64
    implicit none

    integer, intent(in) :: ntime, nseries
    real(real64), intent(in) :: response(ntime, nseries), covariate(ntime, nseries)
    real(real64), intent(out) :: s_values(nseries), var_values(nseries), tau_values(nseries)

    integer :: col
    real(real64), allocatable :: rank_x(:), rank_y(:)

    allocate(rank_x(max(1, ntime)))
    allocate(rank_y(max(1, ntime)))

    do col = 1, nseries
        call compute_partial_stats( &
            response(:, col), covariate(:, col), s_values(col), var_values(col), tau_values(col), rank_x, rank_y &
        )
    end do
end subroutine partial_stats_batch


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
subroutine mk_score_var_sen_batch( &
    data, s_values, var_values, slopes, modified, ntime, nseries &
)
    use mann_kendall_core_mod, only : &
        compute_base_variance, compute_modified_variance_with_slope, &
        compute_s_value_and_sen_slope_with_inv_lag, real64
    implicit none

    integer, intent(in) :: modified, ntime, nseries
    real(real64), intent(in) :: data(ntime, nseries)
    real(real64), intent(out) :: s_values(nseries), var_values(nseries), slopes(nseries)

    real(real64) :: inv_lag(max(1, ntime - 1))
    real(real64) :: centered_work(ntime), detrended_work(ntime), sorted_work(ntime)
    real(real64) :: slope_work(max(1, ntime * (ntime - 1) / 2))
    integer :: col, lag

    do lag = 1, ntime - 1
        inv_lag(lag) = 1.0_real64 / real(lag, real64)
    end do

    do col = 1, nseries
        call compute_s_value_and_sen_slope_with_inv_lag( &
            data(:, col), inv_lag, s_values(col), slopes(col), slope_work &
        )

        if (modified /= 0) then
            call compute_modified_variance_with_slope( &
                data(:, col), sorted_work, centered_work, detrended_work, slopes(col), &
                var_values(col) &
            )
        else
            call compute_base_variance(data(:, col), sorted_work, var_values(col))
        end if
    end do
end subroutine mk_score_var_sen_batch
