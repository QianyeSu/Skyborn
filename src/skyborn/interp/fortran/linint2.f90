! QUICK REFERENCE
! PURPOSE
!    MODERN FREE-FORM REPLACEMENT FOR THE HISTORICAL
!    `archive/linint2.f` NCL INTERPOLATION HELPERS.
!    THIS FILE PRESERVES THE LEGACY PUBLIC ENTRY POINT NAMES:
!
!      - `dlinint1`
!      - `dlinint2`
!      - `dlinint2pts`
!
!    PLUS THE SUPPORT ROUTINES THEY DEPEND ON.
!
! NOTES
!    THE GOAL HERE IS SOURCE MODERNIZATION, NOT AN ALGORITHM CHANGE.
!    EXACT VALUE COMPARISONS, MISSING-VALUE SEMANTICS, CYCLIC-X HANDLING,
!    AND LEGACY ERROR CODES ARE KEPT IN PLACE SO EXISTING CALLERS SEE THE
!    SAME BEHAVIOR.

subroutine dlinint1(nxi, xi, fi, icycx, nxo, xo, fo, xiw, fxiw, nxi2, xmsg, iopt, ier)
    implicit none

    integer :: nxi, nxo, nxi2, iopt, icycx, ier
    double precision :: xi(nxi), fi(nxi)
    double precision :: xo(nxo), fo(nxo), xmsg
    double precision :: xiw(nxi2), fxiw(nxi2)

    integer :: nx, npts, iflag, nxstrt, nxlast

    ier = 0
    if (nxo < 1) then
        ier = 1
        return
    end if

    fo = xmsg

    if (nxi <= 1) then
        ier = 1
        return
    end if

    call dmonoid2(nxi, xi, nxo, xo, iflag, ier)
    if (iflag == 0 .or. ier /= 0) return

    if (icycx == 0) then
        if (iopt == 0) then
            call dlin2int1(nxi, xi, fi, nxo, xo, fo, xmsg, iflag)
        else if (iopt == 1) then
            npts = 0
            do nx = 1, nxi
                if (fi(nx) /= xmsg) then
                    npts = npts + 1
                    xiw(npts + 1) = xi(nx)
                    fxiw(npts + 1) = fi(nx)
                end if
            end do

            call dlin2int1(npts, xiw(2), fxiw(2), nxo, xo, fo, xmsg, iflag)
        end if
    else
        if (iopt == 0) then
            do nx = 1, nxi
                xiw(nx + 1) = xi(nx)
                fxiw(nx + 1) = fi(nx)
            end do

            call dlincyc(nxi, xi, fi, 1, nxi, iflag, xiw, fxiw, nxi)
            call dlin2int1(nxi + 2, xiw, fxiw, nxo, xo, fo, xmsg, iflag)
        else if (iopt == 1) then
            npts = 0
            nxstrt = 0
            nxlast = 0

            do nx = 1, nxi
                if (fi(nx) /= xmsg) then
                    npts = npts + 1
                    xiw(npts + 1) = xi(nx)
                    fxiw(npts + 1) = fi(nx)
                    if (npts == 1) nxstrt = nx
                    nxlast = nx
                end if
            end do

            if (npts == 0) then
                ier = 1
                return
            end if

            call dlincyc(nxi, xi, fi, nxstrt, nxlast, iflag, xiw, fxiw, npts)
            call dlin2int1(npts + 2, xiw, fxiw, nxo, xo, fo, xmsg, iflag)
        end if
    end if
end subroutine dlinint1


subroutine dlinint2(nxi, xi, nyi, yi, fi, icycx, nxo, xo, nyo, yo, fo, xiw, fxiw, nxi2, xmsg, iopt, ier)
    implicit none

    integer :: nxi, nyi, nxo, nyo, nxi2, icycx, iopt, ier
    double precision :: xi(nxi), yi(nyi), fi(nxi, nyi)
    double precision :: xo(nxo), yo(nyo), fo(nxo, nyo), xmsg
    double precision :: xiw(nxi2), fxiw(nxi2)

    integer :: iflag, npts, nxstrt, nxlast
    double precision :: yiw(nyi), fyiw(nyi), foyw(nyo)
    double precision :: ftmp(nxo, nyi)
    integer :: nx, ny

    ier = 0
    if (nxo < 1 .or. nyo < 1) then
        ier = 1
        return
    end if

    fo = xmsg

    if (nxi < 2 .or. nyi < 2) then
        ier = 1
        return
    end if

    call dmonoid2(nxi, xi, nxo, xo, iflag, ier)
    if (iflag == 0 .or. ier /= 0) then
        ier = 2
        return
    end if

    call dmonoid2(nyi, yi, nyo, yo, iflag, ier)
    if (iflag == 0 .or. ier /= 0) then
        ier = 3
        return
    end if

    if (icycx == 0) then
        if (iopt == 0) then
            do ny = 1, nyi
                call dlin2int1(nxi, xi, fi(1, ny), nxo, xo, ftmp(1, ny), xmsg, iflag)
            end do

            do nx = 1, nxo
                do ny = 1, nyi
                    fyiw(ny) = ftmp(nx, ny)
                end do

                call dlin2int1(nyi, yi, fyiw, nyo, yo, foyw, xmsg, iflag)

                do ny = 1, nyo
                    fo(nx, ny) = foyw(ny)
                end do
            end do
        else if (iopt == 1) then
            do ny = 1, nyi
                npts = 0
                do nx = 1, nxi
                    if (fi(nx, ny) /= xmsg) then
                        npts = npts + 1
                        xiw(npts + 1) = xi(nx)
                        fxiw(npts + 1) = fi(nx, ny)
                    end if
                end do

                call dlin2int1(npts, xiw(2), fxiw(2), nxo, xo, ftmp(1, ny), xmsg, iflag)
            end do

            do nx = 1, nxo
                npts = 0
                do ny = 1, nyi
                    if (ftmp(nx, ny) /= xmsg) then
                        npts = npts + 1
                        yiw(npts) = yi(ny)
                        fyiw(npts) = ftmp(nx, ny)
                    end if
                end do

                call dlin2int1(npts, yiw, fyiw, nyo, yo, foyw, xmsg, iflag)

                do ny = 1, nyo
                    fo(nx, ny) = foyw(ny)
                end do
            end do
        end if
    else
        if (iopt == 0) then
            do ny = 1, nyi
                do nx = 1, nxi
                    xiw(nx + 1) = xi(nx)
                    fxiw(nx + 1) = fi(nx, ny)
                end do

                npts = nxi
                call dlincyc(nxi, xi, fi(1, ny), 1, nxi, iflag, xiw, fxiw, npts)
                call dlin2int1(nxi + 2, xiw, fxiw, nxo, xo, ftmp(1, ny), xmsg, iflag)
            end do

            do nx = 1, nxo
                do ny = 1, nyi
                    fyiw(ny) = ftmp(nx, ny)
                end do

                call dlin2int1(nyi, yi, fyiw, nyo, yo, foyw, xmsg, iflag)

                do ny = 1, nyo
                    fo(nx, ny) = foyw(ny)
                end do
            end do
        else if (iopt == 1) then
            do ny = 1, nyi
                npts = 0
                nxstrt = 0
                nxlast = 0

                do nx = 1, nxi
                    if (fi(nx, ny) /= xmsg) then
                        npts = npts + 1
                        xiw(npts + 1) = xi(nx)
                        fxiw(npts + 1) = fi(nx, ny)
                        if (npts == 1) nxstrt = nx
                        nxlast = nx
                    end if
                end do

                call dlincyc(nxi, xi, fi(1, ny), nxstrt, nxlast, iflag, xiw, fxiw, npts)
                call dlin2int1(npts + 2, xiw, fxiw, nxo, xo, ftmp(1, ny), xmsg, iflag)
            end do

            do nx = 1, nxo
                npts = 0
                do ny = 1, nyi
                    if (ftmp(nx, ny) /= xmsg) then
                        npts = npts + 1
                        yiw(npts) = yi(ny)
                        fyiw(npts) = ftmp(nx, ny)
                    end if
                end do

                call dlin2int1(npts, yiw, fyiw, nyo, yo, foyw, xmsg, iflag)

                do ny = 1, nyo
                    fo(nx, ny) = foyw(ny)
                end do
            end do
        end if
    end if
end subroutine dlinint2


subroutine dlin2int1(nin, xi, fi, nout, xo, fo, xmsg, iflag)
    implicit none

    integer :: nin, nout, iflag
    double precision :: xi(nin), fi(nin), xo(nout), fo(nout), xmsg

    integer :: ni, no, nistrt, nis
    double precision :: slope

    fo = xmsg

    ! Exact-match pass first, keeping the historical progressive start index.
    nistrt = 1
    nis = nistrt
    do no = 1, nout
        do ni = nistrt, nin
            if (xo(no) == xi(ni)) then
                fo(no) = fi(ni)
                nis = ni + 1
                exit
            end if
        end do
        nistrt = nis
    end do

    if (iflag == 1) then
        do no = 1, nout
            do ni = 1, nin - 1
                if (xo(no) > xi(ni) .and. xo(no) < xi(ni + 1)) then
                    if (fi(ni) /= xmsg .and. fi(ni + 1) /= xmsg) then
                        slope = (fi(ni + 1) - fi(ni)) / (xi(ni + 1) - xi(ni))
                        fo(no) = fi(ni) + slope * (xo(no) - xi(ni))
                    end if
                end if
            end do
        end do
    else if (iflag == -1) then
        do no = 1, nout
            do ni = 1, nin - 1
                if (xo(no) < xi(ni) .and. xo(no) > xi(ni + 1)) then
                    if (fi(ni) /= xmsg .and. fi(ni + 1) /= xmsg) then
                        slope = (fi(ni + 1) - fi(ni)) / (xi(ni + 1) - xi(ni))
                        fo(no) = fi(ni) + slope * (xo(no) - xi(ni))
                    end if
                end if
            end do
        end do
    end if
end subroutine dlin2int1


subroutine dlincyc(nxi, xi, fi, nxstrt, nxlast, iflag, xiw, fxiw, npts)
    implicit none

    integer :: nxi, nxstrt, nxlast, iflag, npts
    double precision :: xi(nxi), fi(nxi), xiw(0:nxi + 1), fxiw(0:nxi + 1)
    double precision :: dx

    if (nxstrt == 1 .and. nxlast == nxi) then
        dx = abs(xi(2) - xi(1))
        xiw(0) = xi(1) - iflag * dx
        fxiw(0) = fi(nxi)

        dx = abs(xi(nxi) - xi(nxi - 1))
        xiw(npts + 1) = xi(nxi) + iflag * dx
        fxiw(npts + 1) = fi(1)
    else
        dx = (nxstrt + (nxi - nxlast)) * abs(xi(2) - xi(1))
        xiw(0) = xiw(1) - iflag * dx
        fxiw(0) = fxiw(npts)

        dx = ((nxi - nxlast) + nxstrt) * abs(xi(nxi) - xi(nxi - 1))
        xiw(npts + 1) = xiw(npts) + iflag * dx
        fxiw(npts + 1) = fxiw(1)
    end if
end subroutine dlincyc


subroutine dmonoinc(x, nx, ner, ier)
    implicit none

    integer :: nx, ner, ier
    double precision :: x(nx)
    integer :: n

    ier = 0
    if (nx <= 1) return

    do n = 1, nx - 1
        if (x(n + 1) <= x(n)) then
            ier = ner
            return
        end if
    end do
end subroutine dmonoinc


subroutine dmonoid1(nin, xi, iflag, ier)
    implicit none

    integer :: nin, iflag, ier
    double precision :: xi(nin)
    integer :: ni

    ier = 0
    iflag = 0

    if (nin < 2) then
        iflag = 1
        return
    end if

    if (xi(2) > xi(1)) then
        do ni = 1, nin - 1
            if (xi(ni + 1) <= xi(ni)) then
                ier = 2
                return
            end if
        end do
        iflag = 1
    else
        do ni = 1, nin - 1
            if (xi(ni + 1) >= xi(ni)) then
                ier = 2
                return
            end if
        end do
        iflag = -1
    end if
end subroutine dmonoid1


subroutine dmonoid2(nin, xi, nout, xo, iflag, ier)
    implicit none

    integer :: nin, nout, iflag, ier
    double precision :: xi(nin), xo(nout)
    integer :: nflag, ner

    ier = 0
    iflag = 0
    nflag = 0

    call dmonoid1(nin, xi, iflag, ier)
    if (iflag == 0) then
        ier = 2
    else
        call dmonoid1(nout, xo, nflag, ner)
        if (nflag == 0) then
            ier = 3
        else if (iflag /= nflag) then
            ier = 4
        end if
    end if
end subroutine dmonoid2


subroutine dlinint2pts(nxi, xi, nyi, yi, fi, icycx, nxyo, xo, yo, fo, xiw, fixw, nxi2, xmsg, ier)
    implicit none

    integer :: nxi, nyi, nxyo, icycx, nxi2, ier
    double precision :: xi(nxi), yi(nyi), fi(nxi, nyi)
    double precision :: xo(nxyo), yo(nxyo), fo(nxyo)
    double precision :: xiw(nxi2), fixw(nxi2, nyi), xmsg

    double precision :: dx
    integer :: nx, ny, nopt

    nopt = -1
    ier = 0

    if (nxi < 2 .or. nyi < 2) then
        ier = 1
        return
    end if

    call dmonoinc(xi, nxi, 2, ier)
    if (ier /= 0) return

    call dmonoinc(yi, nyi, 3, ier)
    if (ier /= 0) return

    if (icycx == 0) then
        call dlint2xy(nxi, xi, nyi, yi, fi, nxyo, xo, yo, fo, xmsg, nopt, ier)
    else
        do nx = 1, nxi
            xiw(nx + 1) = xi(nx)
        end do

        dx = xi(2) - xi(1)
        xiw(1) = xi(1) - dx
        xiw(nxi2) = xi(nxi) + dx

        do ny = 1, nyi
            do nx = 1, nxi
                fixw(nx + 1, ny) = fi(nx, ny)
            end do
            fixw(1, ny) = fi(nxi, ny)
            fixw(nxi2, ny) = fi(1, ny)
        end do

        call dlint2xy(nxi2, xiw, nyi, yi, fixw, nxyo, xo, yo, fo, xmsg, nopt, ier)
    end if
end subroutine dlinint2pts


subroutine dlint2xy(nxi, xi, nyi, yi, fi, nxyo, xo, yo, fo, xmsg, nopt, ier)
    implicit none

    integer :: nxi, nyi, nxyo, nopt, ier
    double precision :: xi(nxi), yi(nyi), fi(nxi, nyi), xmsg
    double precision :: xo(nxyo), yo(nxyo), fo(nxyo)

    integer :: n, m, nxy, nn, mm
    double precision :: tmp1, tmp2, slpx, slpy

    do nxy = 1, nxyo
        fo(nxy) = xmsg

        nn = 0
        do n = 1, nxi - 1
            if (xo(nxy) >= xi(n) .and. xo(nxy) < xi(n + 1)) then
                nn = n
                exit
            end if
        end do

        mm = 0
        do m = 1, nyi - 1
            if (yo(nxy) >= yi(m) .and. yo(nxy) < yi(m + 1)) then
                mm = m
                exit
            end if
        end do

        if (nn /= 0 .and. mm /= 0) then
            if (xo(nxy) == xi(nn) .and. yo(nxy) == yi(mm)) then
                fo(nxy) = fi(nn, mm)
            else if (fi(nn, mm) /= xmsg .and. fi(nn + 1, mm) /= xmsg .and. &
                     fi(nn, mm + 1) /= xmsg .and. fi(nn + 1, mm + 1) /= xmsg) then
                slpx = (xo(nxy) - xi(nn)) / (xi(nn + 1) - xi(nn))
                slpy = (yo(nxy) - yi(mm)) / (yi(mm + 1) - yi(mm))

                tmp1 = fi(nn, mm) + slpx * (fi(nn + 1, mm) - fi(nn, mm))
                tmp2 = fi(nn, mm + 1) + slpx * (fi(nn + 1, mm + 1) - fi(nn, mm + 1))
                fo(nxy) = tmp1 + slpy * (tmp2 - tmp1)
            else if (nopt == -1) then
                call estfow( &
                    fi(nn, mm), fi(nn + 1, mm), fi(nn, mm + 1), fi(nn + 1, mm + 1), &
                    xi(nn), xi(nn + 1), yi(mm), yi(mm + 1), &
                    fo(nxy), xo(nxy), yo(nxy), xmsg &
                )
            end if
        end if
    end do
end subroutine dlint2xy


subroutine estfow(f1, f2, f3, f4, x1, x2, y1, y2, f0, x0, y0, xmsg)
    implicit none

    double precision :: f1, f2, f3, f4, x1, x2, y1, y2, f0, x0, y0, xmsg
    double precision :: fi(2, 2), w(2, 2), sum, swt
    integer :: n, m

    f0 = xmsg
    if (f1 == xmsg .and. f2 == xmsg .and. f3 == xmsg .and. f4 == xmsg) return

    w(1, 1) = 1.0d0 / sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    w(2, 1) = 1.0d0 / sqrt((x2 - x0) ** 2 + (y1 - y0) ** 2)
    w(1, 2) = 1.0d0 / sqrt((x1 - x0) ** 2 + (y2 - y0) ** 2)
    w(2, 2) = 1.0d0 / sqrt((x2 - x0) ** 2 + (y2 - y0) ** 2)

    fi(1, 1) = f1
    fi(2, 1) = f2
    fi(1, 2) = f3
    fi(2, 2) = f4

    sum = 0.0d0
    swt = 0.0d0
    do m = 1, 2
        do n = 1, 2
            if (fi(n, m) /= xmsg) then
                sum = sum + fi(n, m) * w(n, m)
                swt = swt + w(n, m)
            end if
        end do
    end do

    if (swt > 0.0d0) then
        f0 = sum / swt
    end if
end subroutine estfow
