module linmsg_core
    implicit none

    integer, parameter :: real64 = selected_real_kind(15, 307)
contains

    ! QUICK REFERENCE
    ! PURPOSE
    !    FILL SHORT MISSING RUNS IN A 1-D SERIES USING LINEAR
    !    INTERPOLATION, PRESERVING THE HISTORICAL NCL/Fortran
    !    DLINMSG EDGE AND GAP-LENGTH SEMANTICS.
    !
    ! INPUTS
    !    X(NPTS)  - SERIES TO UPDATE IN PLACE
    !    XMSG     - MISSING VALUE SENTINEL
    !    MFLAG    - IF < 0, LEADING/TRAILING GAPS USE NEAREST VALID VALUE;
    !               OTHERWISE THEY REMAIN MISSING
    !    MPTCRT   - GAPS LONGER THAN ABS(MPTCRT) ARE LEFT UNFILLED
    !
    ! OUTPUT
    !    X(NPTS) IS UPDATED IN PLACE. A SERIES THAT IS ENTIRELY MISSING IS
    !    LEFT UNCHANGED TO MATCH THE LEGACY ROUTINE.
    subroutine dlinmsg(x, npts, xmsg, mflag, mptcrt)
        integer, intent(in) :: npts, mflag, mptcrt
        real(real64), intent(inout) :: x(npts)
        real(real64), intent(in) :: xmsg

        integer :: n, nbase, nend, nn, nptcrt, nstrt
        real(real64) :: slope
        logical :: seen_valid

        if (npts <= 0) return

        nstrt = 0
        nend = 0
        nptcrt = abs(mptcrt)
        seen_valid = .false.

        do n = 1, npts
            if (x(n) == xmsg) then
                if (nstrt == 0) nstrt = n
                nend = n
            else
                seen_valid = .true.

                if (nstrt /= 0) then
                    if ((nend - nstrt + 1) <= nptcrt) then
                        if (nstrt == 1) then
                            if (mflag < 0) then
                                x(nstrt:nend) = x(n)
                            else
                                x(nstrt:nend) = xmsg
                            end if
                        else
                            nbase = nstrt - 1
                            slope = (x(n) - x(nbase)) / real(n - nbase, real64)
                            do nn = nstrt, nend
                                x(nn) = x(nbase) + slope * real(nn - nbase, real64)
                            end do
                        end if
                    end if

                    nstrt = 0
                    nend = 0
                end if
            end if
        end do

        if (nstrt /= 0 .and. nend == npts) then
            if (.not. seen_valid) return

            if (mflag < 0) then
                x(nstrt:npts) = x(nstrt - 1)
            else
                x(nstrt:npts) = xmsg
            end if
        end if
    end subroutine dlinmsg

    ! QUICK REFERENCE
    ! PURPOSE
    !    FAST PATH FOR THE `rcm2rgrid` ROW CLEANUP STEP. THIS ONLY FILLS
    !    INTERIOR GAPS OF LENGTH 1 OR 2 AND LEAVES LEADING / TRAILING GAPS
    !    MISSING, WHICH MATCHES THE EFFECT OF
    !    DLINMSG(X, NPTS, XMSG, MFLAG=0, MPTCRT=2).
    !
    ! NOTE
    !    THIS IS AN ADDITIVE OPTIMIZATION HELPER. IT DOES NOT REPLACE THE
    !    FULL LEGACY-COMPATIBLE `dlinmsg(...)` ROUTINE ABOVE.
    subroutine dlinmsg_interior_short_gaps(x, npts, xmsg)
        integer, intent(in) :: npts
        real(real64), intent(inout) :: x(npts)
        real(real64), intent(in) :: xmsg

        integer :: gap_len, gap_start, n
        real(real64) :: left, right, step

        if (npts <= 2) return

        n = 2
        do while (n <= npts - 1)
            if (x(n) /= xmsg) then
                n = n + 1
                cycle
            end if

            gap_start = n
            do while (n <= npts .and. x(n) == xmsg)
                n = n + 1
            end do

            gap_len = n - gap_start
            if (n > npts .or. gap_len > 2) cycle
            if (x(gap_start - 1) == xmsg .or. x(n) == xmsg) cycle

            left = x(gap_start - 1)
            right = x(n)

            if (gap_len == 1) then
                x(gap_start) = 0.5_real64 * (left + right)
            else
                step = (right - left) / 3.0_real64
                x(gap_start) = left + step
                x(gap_start + 1) = left + 2.0_real64 * step
            end if
        end do
    end subroutine dlinmsg_interior_short_gaps

end module linmsg_core
