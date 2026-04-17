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

end module linmsg_core
