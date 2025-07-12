!
! HRFFT - Real FFT routines for SPHEREPACK
! Converted from FORTRAN 77 to modern Fortran
! Simplified version for compatibility with older gfortran
!
! Original SPHEREPACK copyright notice:
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 1998 by UCAR                 .
!  .                                                             .
!  .       University Corporation for Atmospheric Research       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                         SPHEREPACK                          .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

module hrfft_mod
   implicit none
   private

   ! Public interfaces
   public :: hrffti, hrfftf, hrfftb

contains

subroutine hrffti(n,wsave)
    implicit none
    integer, intent(in) :: n
    real, intent(out) :: wsave(n+15)

    if (n == 1) return
    call hrfti1(n,wsave(1:n),wsave(n+1:n+15))
end subroutine hrffti

subroutine hrfti1(n,wa,ifac)
    implicit none
    integer, intent(in) :: n
    real, intent(out) :: wa(n), ifac(15)

    integer, parameter :: ntryh(4) = [4, 2, 3, 5]
    integer :: nl, nf, j, ntry, nq, nr, ip, ld, l2, ido, ipm, i, ii, is
    integer :: k1, ib, nfm1, l1
    real :: argh, argld, arg, fi

    nl = n
    nf = 0
    j = 0

    ! Factor n
    do
        j = j + 1
        if (j <= 4) then
            ntry = ntryh(j)
        else
            ntry = ntry + 2
        end if

        nq = nl / ntry
        nr = nl - ntry * nq

        if (nr == 0) then
            nf = nf + 1
            ifac(nf + 2) = real(ntry)
            nl = nq

            if (ntry == 2 .and. nf > 1) then
                do i = 2, nf
                    ib = nf - i + 2
                    ifac(ib + 2) = ifac(ib + 1)
                end do
                ifac(3) = 2.0
            end if
        end if

        if (nl == 1) exit
    end do

    ifac(1) = real(n)
    ifac(2) = real(nf)

    argh = 8.0 * atan(1.0) / real(n)
    is = 0
    nfm1 = nf - 1
    l1 = 1

    if (nfm1 == 0) return

    do k1 = 1, nfm1
        ip = int(ifac(k1 + 2))
        ld = 0
        l2 = l1 * ip
        ido = n / l2
        ipm = ip - 1

        do j = 1, ipm
            ld = ld + l1
            i = is
            argld = real(ld) * argh
            fi = 0.0

            do ii = 3, ido, 2
                i = i + 2
                fi = fi + 1.0
                arg = fi * argld
                wa(i - 1) = cos(arg)
                wa(i) = sin(arg)
            end do
            is = is + ido
        end do
        l1 = l2
    end do
end subroutine hrfti1

subroutine hrfftf(m,n,r,mdimr,whrfft,work)
    implicit none
    integer, intent(in) :: m, n, mdimr
    real, intent(inout) :: r(mdimr,n)
    real, intent(in) :: whrfft(n+15)
    real, intent(inout) :: work(*)

    if (n == 1) return
    call hrftf1(m,n,r,mdimr,work,whrfft,whrfft(n+1:))
end subroutine hrfftf

subroutine hrftf1(m,n,c,mdimc,ch,wa,ifac)
    implicit none
    integer, intent(in) :: m, n, mdimc
    real, intent(inout) :: c(mdimc,n)
    real, intent(inout) :: ch(m,n)
    real, intent(in) :: wa(n), ifac(15)

    integer :: nf, na, l2, iw, k1, kh, ip, l1, ido, idl1
    integer :: ix2, ix3, ix4, j, i

    nf = int(ifac(2))
    na = 1
    l2 = n
    iw = n

    do k1 = 1, nf
        kh = nf - k1
        ip = int(ifac(kh + 3))
        l1 = l2 / ip
        ido = n / l2
        idl1 = ido * l1
        iw = iw - (ip - 1) * ido
        na = 1 - na

        select case (ip)
        case (4)
            ix2 = iw + ido
            ix3 = ix2 + ido
            if (na == 0) then
                call hradf4(m,ido,l1,c,mdimc,ch,m,wa(iw:),wa(ix2:),wa(ix3:))
            else
                call hradf4(m,ido,l1,ch,m,c,mdimc,wa(iw:),wa(ix2:),wa(ix3:))
            end if
        case (2)
            if (na == 0) then
                call hradf2(m,ido,l1,c,mdimc,ch,m,wa(iw:))
            else
                call hradf2(m,ido,l1,ch,m,c,mdimc,wa(iw:))
            end if
        case (3)
            ix2 = iw + ido
            if (na == 0) then
                call hradf3(m,ido,l1,c,mdimc,ch,m,wa(iw:),wa(ix2:))
            else
                call hradf3(m,ido,l1,ch,m,c,mdimc,wa(iw:),wa(ix2:))
            end if
        case (5)
            ix2 = iw + ido
            ix3 = ix2 + ido
            ix4 = ix3 + ido
            if (na == 0) then
                call hradf5(m,ido,l1,c,mdimc,ch,m,wa(iw:),wa(ix2:),wa(ix3:),wa(ix4:))
            else
                call hradf5(m,ido,l1,ch,m,c,mdimc,wa(iw:),wa(ix2:),wa(ix3:),wa(ix4:))
            end if
        case default
            if (ido == 1) na = 1 - na
            if (na == 0) then
                call hradfg(m,ido,ip,l1,idl1,c,c,c,mdimc,ch,ch,m,wa(iw:))
                na = 1
            else
                call hradfg(m,ido,ip,l1,idl1,ch,ch,ch,m,c,c,mdimc,wa(iw:))
                na = 0
            end if
        end select

        l2 = l1
    end do

    if (na == 1) return

    do j = 1, n
        do i = 1, m
            c(i,j) = ch(i,j)
        end do
    end do
end subroutine hrftf1

subroutine hrfftb(m,n,r,mdimr,whrfft,work)
    implicit none
    integer, intent(in) :: m, n, mdimr
    real, intent(inout) :: r(mdimr,n)
    real, intent(in) :: whrfft(n+15)
    real, intent(inout) :: work(*)

    if (n == 1) return
    call hrftb1(m,n,r,mdimr,work,whrfft,whrfft(n+1:))
end subroutine hrfftb

subroutine hrftb1(m,n,c,mdimc,ch,wa,ifac)
    implicit none
    integer, intent(in) :: m, n, mdimc
    real, intent(inout) :: c(mdimc,n)
    real, intent(inout) :: ch(m,n)
    real, intent(in) :: wa(n), ifac(15)

    integer :: nf, na, l1, iw, k1, ip, l2, ido, idl1
    integer :: ix2, ix3, ix4, j, i

    nf = int(ifac(2))
    na = 0
    l1 = 1
    iw = 1

    do k1 = 1, nf
        ip = int(ifac(k1 + 2))
        l2 = ip * l1
        ido = n / l2
        idl1 = ido * l1

        select case (ip)
        case (4)
            ix2 = iw + ido
            ix3 = ix2 + ido
            if (na == 0) then
                call hradb4(m,ido,l1,c,mdimc,ch,m,wa(iw:),wa(ix2:),wa(ix3:))
            else
                call hradb4(m,ido,l1,ch,m,c,mdimc,wa(iw:),wa(ix2:),wa(ix3:))
            end if
            na = 1 - na
        case (2)
            if (na == 0) then
                call hradb2(m,ido,l1,c,mdimc,ch,m,wa(iw:))
            else
                call hradb2(m,ido,l1,ch,m,c,mdimc,wa(iw:))
            end if
            na = 1 - na
        case (3)
            ix2 = iw + ido
            if (na == 0) then
                call hradb3(m,ido,l1,c,mdimc,ch,m,wa(iw:),wa(ix2:))
            else
                call hradb3(m,ido,l1,ch,m,c,mdimc,wa(iw:),wa(ix2:))
            end if
            na = 1 - na
        case (5)
            ix2 = iw + ido
            ix3 = ix2 + ido
            ix4 = ix3 + ido
            if (na == 0) then
                call hradb5(m,ido,l1,c,mdimc,ch,m,wa(iw:),wa(ix2:),wa(ix3:),wa(ix4:))
            else
                call hradb5(m,ido,l1,ch,m,c,mdimc,wa(iw:),wa(ix2:),wa(ix3:),wa(ix4:))
            end if
            na = 1 - na
        case default
            if (na == 0) then
                call hradbg(m,ido,ip,l1,idl1,c,c,c,mdimc,ch,ch,m,wa(iw:))
            else
                call hradbg(m,ido,ip,l1,idl1,ch,ch,ch,m,c,c,mdimc,wa(iw:))
            end if
            if (ido == 1) na = 1 - na
        end select

        l1 = l2
        iw = iw + (ip - 1) * ido
    end do

    if (na == 0) return

    do j = 1, n
        do i = 1, m
            c(i,j) = ch(i,j)
        end do
    end do
end subroutine hrftb1

subroutine hradf4(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2,wa3)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,l1,4)
    real, intent(out) :: ch(mdimch,ido,4,l1)
    real, intent(in) :: wa1(ido), wa2(ido), wa3(ido)

    integer :: k, m, i, ic, idp2
    real :: tr1, tr2, tr3, tr4, ti1, ti2, ti3, ti4
    real :: cr2, cr3, cr4, ci2, ci3, ci4

    do k = 1, l1
        do m = 1, mp
            tr1 = cc(m,1,k,2) + cc(m,1,k,4)
            tr2 = cc(m,1,k,1) + cc(m,1,k,3)
            ch(m,1,1,k) = tr1 + tr2
            ch(m,ido,4,k) = tr2 - tr1
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,3)
            ch(m,1,3,k) = cc(m,1,k,4) - cc(m,1,k,2)
        end do
    end do

    if (ido < 2) return
    if (ido == 2) then
        do k = 1, l1
            do m = 1, mp
                ti1 = -cc(m,ido,k,2) - cc(m,ido,k,4)
                tr1 = cc(m,ido,k,2) - cc(m,ido,k,4)
                ch(m,ido,1,k) = tr1 + cc(m,ido,k,1)
                ch(m,ido,3,k) = cc(m,ido,k,1) - tr1
                ch(m,1,2,k) = ti1 - cc(m,ido,k,3)
                ch(m,1,4,k) = ti1 + cc(m,ido,k,3)
            end do
        end do
        return
    end if

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                cr2 = wa1(i-2)*cc(m,i-1,k,2) + wa1(i-1)*cc(m,i,k,2)
                ci2 = wa1(i-2)*cc(m,i,k,2) - wa1(i-1)*cc(m,i-1,k,2)
                cr3 = wa2(i-2)*cc(m,i-1,k,3) + wa2(i-1)*cc(m,i,k,3)
                ci3 = wa2(i-2)*cc(m,i,k,3) - wa2(i-1)*cc(m,i-1,k,3)
                cr4 = wa3(i-2)*cc(m,i-1,k,4) + wa3(i-1)*cc(m,i,k,4)
                ci4 = wa3(i-2)*cc(m,i,k,4) - wa3(i-1)*cc(m,i-1,k,4)
                tr1 = cr2 + cr4
                tr4 = cr4 - cr2
                ti1 = ci2 + ci4
                ti4 = ci2 - ci4
                ti2 = cc(m,i,k,1) + ci3
                ti3 = cc(m,i,k,1) - ci3
                tr2 = cc(m,i-1,k,1) + cr3
                tr3 = cc(m,i-1,k,1) - cr3
                ch(m,i-1,1,k) = tr1 + tr2
                ch(m,ic-1,4,k) = tr2 - tr1
                ch(m,i,1,k) = ti1 + ti2
                ch(m,ic,4,k) = ti1 - ti2
                ch(m,i-1,3,k) = tr4 + tr3
                ch(m,ic-1,2,k) = tr3 - tr4
                ch(m,i,3,k) = ti4 + ti3
                ch(m,ic,2,k) = ti4 - ti3
            end do
        end do
    end do
end subroutine hradf4

subroutine hradf2(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,l1,2)
    real, intent(out) :: ch(mdimch,ido,2,l1)
    real, intent(in) :: wa1(ido)

    integer :: k, m, i, ic, idp2
    real :: tr2, ti2

    do k = 1, l1
        do m = 1, mp
            ch(m,1,1,k) = cc(m,1,k,1) + cc(m,1,k,2)
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,2)
        end do
    end do

    if (ido < 2) return
    if (ido == 2) then
        do k = 1, l1
            do m = 1, mp
                ch(m,1,2,k) = -cc(m,ido,k,2)
                ch(m,ido,1,k) = cc(m,ido,k,1)
            end do
        end do
        return
    end if

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                tr2 = wa1(i-2)*cc(m,i-1,k,2) + wa1(i-1)*cc(m,i,k,2)
                ti2 = wa1(i-2)*cc(m,i,k,2) - wa1(i-1)*cc(m,i-1,k,2)
                ch(m,i,1,k) = cc(m,i,k,1) + ti2
                ch(m,ic,2,k) = ti2 - cc(m,i,k,1)
                ch(m,i-1,1,k) = cc(m,i-1,k,1) + tr2
                ch(m,ic-1,2,k) = cc(m,i-1,k,1) - tr2
            end do
        end do
    end do
end subroutine hradf2

subroutine hradf3(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,l1,3)
    real, intent(out) :: ch(mdimch,ido,3,l1)
    real, intent(in) :: wa1(ido), wa2(ido)

    integer :: k, m, i, ic, idp2
    real, parameter :: taur = -0.5
    real, parameter :: taui = 0.866025403784439
    real :: cr2, ci2, cr3, ci3, dr2, dr3, di2, di3

    do k = 1, l1
        do m = 1, mp
            cr2 = cc(m,1,k,2) + cc(m,1,k,3)
            ch(m,1,1,k) = cc(m,1,k,1) + cr2
            ch(m,1,3,k) = taui*(cc(m,1,k,3) - cc(m,1,k,2))
            ch(m,ido,2,k) = cc(m,1,k,1) + taur*cr2
        end do
    end do

    if (ido == 1) return

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                cr2 = wa1(i-2)*cc(m,i-1,k,2) + wa1(i-1)*cc(m,i,k,2)
                ci2 = wa1(i-2)*cc(m,i,k,2) - wa1(i-1)*cc(m,i-1,k,2)
                cr3 = wa2(i-2)*cc(m,i-1,k,3) + wa2(i-1)*cc(m,i,k,3)
                ci3 = wa2(i-2)*cc(m,i,k,3) - wa2(i-1)*cc(m,i-1,k,3)
                dr2 = cr2 + cr3
                dr3 = cr2 - cr3
                di2 = ci2 + ci3
                di3 = ci2 - ci3
                ch(m,i-1,1,k) = cc(m,i-1,k,1) + dr2
                ch(m,i,1,k) = cc(m,i,k,1) + di2
                ch(m,i-1,3,k) = (cc(m,i-1,k,1) + taur*dr2) + (taui*di3)
                ch(m,ic-1,2,k) = (cc(m,i-1,k,1) + taur*dr2) - (taui*di3)
                ch(m,i,3,k) = (cc(m,i,k,1) + taur*di2) + (taui*dr3)
                ch(m,ic,2,k) = (taui*dr3) - (cc(m,i,k,1) + taur*di2)
            end do
        end do
    end do
end subroutine hradf3

subroutine hradf5(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2,wa3,wa4)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,l1,5)
    real, intent(out) :: ch(mdimch,ido,5,l1)
    real, intent(in) :: wa1(ido), wa2(ido), wa3(ido), wa4(ido)

    integer :: k, m, i, ic, idp2
    real, parameter :: tr11 = 0.309016994374947
    real, parameter :: ti11 = 0.951056516295154
    real, parameter :: tr12 = -0.809016994374947
    real, parameter :: ti12 = 0.587785252292473
    real :: cr2, ci2, cr3, ci3, cr4, ci4, cr5, ci5
    real :: dr2, dr3, dr4, dr5, di2, di3, di4, di5

    do k = 1, l1
        do m = 1, mp
            cr2 = cc(m,1,k,5) + cc(m,1,k,2)
            ci5 = cc(m,1,k,5) - cc(m,1,k,2)
            cr3 = cc(m,1,k,4) + cc(m,1,k,3)
            ci4 = cc(m,1,k,4) - cc(m,1,k,3)
            ch(m,1,1,k) = cc(m,1,k,1) + cr2 + cr3
            ch(m,ido,2,k) = cc(m,1,k,1) + tr11*cr2 + tr12*cr3
            ch(m,1,3,k) = ti11*ci5 + ti12*ci4
            ch(m,ido,4,k) = cc(m,1,k,1) + tr12*cr2 + tr11*cr3
            ch(m,1,5,k) = ti12*ci5 - ti11*ci4
        end do
    end do

    if (ido == 1) return

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                cr2 = wa1(i-2)*cc(m,i-1,k,2) + wa1(i-1)*cc(m,i,k,2)
                ci2 = wa1(i-2)*cc(m,i,k,2) - wa1(i-1)*cc(m,i-1,k,2)
                cr3 = wa2(i-2)*cc(m,i-1,k,3) + wa2(i-1)*cc(m,i,k,3)
                ci3 = wa2(i-2)*cc(m,i,k,3) - wa2(i-1)*cc(m,i-1,k,3)
                cr4 = wa3(i-2)*cc(m,i-1,k,4) + wa3(i-1)*cc(m,i,k,4)
                ci4 = wa3(i-2)*cc(m,i,k,4) - wa3(i-1)*cc(m,i-1,k,4)
                cr5 = wa4(i-2)*cc(m,i-1,k,5) + wa4(i-1)*cc(m,i,k,5)
                ci5 = wa4(i-2)*cc(m,i,k,5) - wa4(i-1)*cc(m,i-1,k,5)
                dr2 = cr2 + cr5
                dr5 = cr5 - cr2
                di2 = ci2 + ci5
                di5 = ci2 - ci5
                dr3 = cr3 + cr4
                dr4 = cr4 - cr3
                di3 = ci3 + ci4
                di4 = ci3 - ci4
                ch(m,i-1,1,k) = cc(m,i-1,k,1) + dr2 + dr3
                ch(m,i,1,k) = cc(m,i,k,1) + di2 + di3
                ch(m,i-1,3,k) = (cc(m,i-1,k,1) + tr11*dr2 + tr12*dr3) + (ti11*di5 + ti12*di4)
                ch(m,ic-1,2,k) = (cc(m,i-1,k,1) + tr11*dr2 + tr12*dr3) - (ti11*di5 + ti12*di4)
                ch(m,i,3,k) = (cc(m,i,k,1) + tr11*di2 + tr12*di3) + (ti11*dr5 + ti12*dr4)
                ch(m,ic,2,k) = (ti11*dr5 + ti12*dr4) - (cc(m,i,k,1) + tr11*di2 + tr12*di3)
                ch(m,i-1,5,k) = (cc(m,i-1,k,1) + tr12*dr2 + tr11*dr3) + (ti12*di5 - ti11*di4)
                ch(m,ic-1,4,k) = (cc(m,i-1,k,1) + tr12*dr2 + tr11*dr3) - (ti12*di5 - ti11*di4)
                ch(m,i,5,k) = (cc(m,i,k,1) + tr12*di2 + tr11*di3) + (ti12*dr5 - ti11*dr4)
                ch(m,ic,4,k) = (ti12*dr5 - ti11*dr4) - (cc(m,i,k,1) + tr12*di2 + tr11*di3)
            end do
        end do
    end do
end subroutine hradf5

subroutine hradfg(mp,ido,ip,l1,idl1,cc,c1,c2,mdimcc,ch,ch2,mdimch,wa)
    implicit none
    integer, intent(in) :: mp, ido, ip, l1, idl1, mdimcc, mdimch
    real, intent(inout) :: cc(mdimcc,ido,ip,l1)
    real, intent(inout) :: c1(mdimcc,ido,l1,ip)
    real, intent(inout) :: c2(mdimcc,idl1,ip)
    real, intent(inout) :: ch(mdimch,ido,l1,ip)
    real, intent(inout) :: ch2(mdimch,idl1,ip)
    real, intent(in) :: wa(ido)

    integer :: ipph, ipp2, idp2, nbd, is, j, k, i, ik, m, idij, jc, l, lc, j2, ic
    real :: dcp, dsp, arg, ar1, ai1, ar1h, ar2, ai2, ar2h, dc2, ds2

    arg = 8.0*atan(1.0) / real(ip)
    dcp = cos(arg)
    dsp = sin(arg)
    ipph = (ip + 1) / 2
    ipp2 = ip + 2
    idp2 = ido + 2
    nbd = (ido - 1) / 2

    if (ido == 1) then
        do ik = 1, idl1
            do m = 1, mp
                c2(m,ik,1) = ch2(m,ik,1)
            end do
        end do
    else
        do ik = 1, idl1
            do m = 1, mp
                ch2(m,ik,1) = c2(m,ik,1)
            end do
        end do

        do j = 2, ip
            do k = 1, l1
                do m = 1, mp
                    ch(m,1,k,j) = c1(m,1,k,j)
                end do
            end do
        end do

        if (nbd <= l1) then
            is = -ido
            do j = 2, ip
                is = is + ido
                idij = is
                do i = 3, ido, 2
                    idij = idij + 2
                    do k = 1, l1
                        do m = 1, mp
                            ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) + wa(idij)*c1(m,i,k,j)
                            ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) - wa(idij)*c1(m,i-1,k,j)
                        end do
                    end do
                end do
            end do
        else
            is = -ido
            do j = 2, ip
                is = is + ido
                do k = 1, l1
                    idij = is
                    do i = 3, ido, 2
                        idij = idij + 2
                        do m = 1, mp
                            ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) + wa(idij)*c1(m,i,k,j)
                            ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) - wa(idij)*c1(m,i-1,k,j)
                        end do
                    end do
                end do
            end do
        end if

        if (nbd >= l1) then
            do j = 2, ipph
                jc = ipp2 - j
                do k = 1, l1
                    do i = 3, ido, 2
                        do m = 1, mp
                            c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                            c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                            c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                            c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
                        end do
                    end do
                end do
            end do
        else
            do j = 2, ipph
                jc = ipp2 - j
                do i = 3, ido, 2
                    do k = 1, l1
                        do m = 1, mp
                            c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                            c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                            c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                            c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
                        end do
                    end do
                end do
            end do
        end if
    end if

    do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
            do m = 1, mp
                c1(m,1,k,j) = ch(m,1,k,j) + ch(m,1,k,jc)
                c1(m,1,k,jc) = ch(m,1,k,jc) - ch(m,1,k,j)
            end do
        end do
    end do

    ar1 = 1.0
    ai1 = 0.0

    do l = 2, ipph
        lc = ipp2 - l
        ar1h = dcp*ar1 - dsp*ai1
        ai1 = dcp*ai1 + dsp*ar1
        ar1 = ar1h

        do ik = 1, idl1
            do m = 1, mp
                ch2(m,ik,l) = c2(m,ik,1) + ar1*c2(m,ik,2)
                ch2(m,ik,lc) = ai1*c2(m,ik,ip)
            end do
        end do

        dc2 = ar1
        ds2 = ai1
        ar2 = ar1
        ai2 = ai1

        do j = 3, ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h

            do ik = 1, idl1
                do m = 1, mp
                    ch2(m,ik,l) = ch2(m,ik,l) + ar2*c2(m,ik,j)
                    ch2(m,ik,lc) = ch2(m,ik,lc) + ai2*c2(m,ik,jc)
                end do
            end do
        end do
    end do

    do j = 2, ipph
        do ik = 1, idl1
            do m = 1, mp
                ch2(m,ik,1) = ch2(m,ik,1) + c2(m,ik,j)
            end do
        end do
    end do

    if (ido >= l1) then
        do k = 1, l1
            do i = 1, ido
                do m = 1, mp
                    cc(m,i,1,k) = ch(m,i,k,1)
                end do
            end do
        end do
    else
        do i = 1, ido
            do k = 1, l1
                do m = 1, mp
                    cc(m,i,1,k) = ch(m,i,k,1)
                end do
            end do
        end do
    end if

    do j = 2, ipph
        jc = ipp2 - j
        j2 = j + j
        do k = 1, l1
            do m = 1, mp
                cc(m,ido,j2-2,k) = ch(m,1,k,j)
                cc(m,1,j2-1,k) = ch(m,1,k,jc)
            end do
        end do
    end do

    if (ido == 1) return

    if (nbd >= l1) then
        do j = 2, ipph
            jc = ipp2 - j
            j2 = j + j
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    do m = 1, mp
                        cc(m,i-1,j2-1,k) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                        cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j) - ch(m,i-1,k,jc)
                        cc(m,i,j2-1,k) = ch(m,i,k,j) + ch(m,i,k,jc)
                        cc(m,ic,j2-2,k) = ch(m,i,k,jc) - ch(m,i,k,j)
                    end do
                end do
            end do
        end do
    else
        do j = 2, ipph
            jc = ipp2 - j
            j2 = j + j
            do i = 3, ido, 2
                ic = idp2 - i
                do k = 1, l1
                    do m = 1, mp
                        cc(m,i-1,j2-1,k) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                        cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j) - ch(m,i-1,k,jc)
                        cc(m,i,j2-1,k) = ch(m,i,k,j) + ch(m,i,k,jc)
                        cc(m,ic,j2-2,k) = ch(m,i,k,jc) - ch(m,i,k,j)
                    end do
                end do
            end do
        end do
    end if
end subroutine hradfg

subroutine hradb4(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2,wa3)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,4,l1)
    real, intent(out) :: ch(mdimch,ido,l1,4)
    real, intent(in) :: wa1(ido), wa2(ido), wa3(ido)

    integer :: k, m, i, ic, idp2
    real :: tr1, tr2, tr3, tr4, ti1, ti2, ti3, ti4
    real :: cr2, cr3, cr4, ci2, ci3, ci4

    do k = 1, l1
        do m = 1, mp
            tr1 = cc(m,1,1,k) - cc(m,ido,4,k)
            tr2 = cc(m,1,1,k) + cc(m,ido,4,k)
            tr4 = cc(m,ido,2,k) + cc(m,ido,2,k)
            tr3 = cc(m,1,3,k) + cc(m,1,3,k)
            ch(m,1,k,1) = tr2 + tr3
            ch(m,1,k,3) = tr2 - tr3
            ch(m,1,k,2) = tr1 + tr4
            ch(m,1,k,4) = tr1 - tr4
        end do
    end do

    if (ido < 2) return
    if (ido == 2) then
        do k = 1, l1
            do m = 1, mp
                ti1 = cc(m,1,2,k) + cc(m,1,4,k)
                ti2 = cc(m,1,2,k) - cc(m,1,4,k)
                tr1 = cc(m,ido,1,k) - cc(m,ido,3,k)
                tr2 = cc(m,ido,1,k) + cc(m,ido,3,k)
                ch(m,ido,k,1) = tr2 + tr2
                ch(m,ido,k,2) = sqrt(2.0)*(tr1 - ti1)
                ch(m,ido,k,3) = ti2 + ti2
                ch(m,ido,k,4) = -sqrt(2.0)*(tr1 + ti1)
            end do
        end do
        return
    end if

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                tr1 = cc(m,i-1,1,k) + cc(m,ic-1,4,k)
                tr2 = cc(m,i-1,1,k) - cc(m,ic-1,4,k)
                ti1 = cc(m,i,1,k) - cc(m,ic,4,k)
                ti2 = cc(m,i,1,k) + cc(m,ic,4,k)
                ti3 = cc(m,i,3,k) + cc(m,ic,2,k)
                tr3 = cc(m,i-1,3,k) + cc(m,ic-1,2,k)
                tr4 = cc(m,i-1,3,k) - cc(m,ic-1,2,k)
                ti4 = cc(m,i,3,k) - cc(m,ic,2,k)
                ch(m,i-1,k,1) = tr1 + tr3
                ch(m,i,k,1) = ti1 + ti3
                cr2 = tr2 - tr4
                ci2 = ti2 - ti4
                ch(m,i-1,k,2) = wa1(i-2)*cr2 - wa1(i-1)*ci2
                ch(m,i,k,2) = wa1(i-2)*ci2 + wa1(i-1)*cr2
                cr3 = tr1 - tr3
                ci3 = ti1 - ti3
                ch(m,i-1,k,3) = wa2(i-2)*cr3 - wa2(i-1)*ci3
                ch(m,i,k,3) = wa2(i-2)*ci3 + wa2(i-1)*cr3
                cr4 = tr2 + tr4
                ci4 = ti2 + ti4
                ch(m,i-1,k,4) = wa3(i-2)*cr4 - wa3(i-1)*ci4
                ch(m,i,k,4) = wa3(i-2)*ci4 + wa3(i-1)*cr4
            end do
        end do
    end do
end subroutine hradb4

subroutine hradb2(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,2,l1)
    real, intent(out) :: ch(mdimch,ido,l1,2)
    real, intent(in) :: wa1(ido)

    integer :: k, m, i, ic, idp2
    real :: tr2, ti2

    do k = 1, l1
        do m = 1, mp
            ch(m,1,k,1) = cc(m,1,1,k) + cc(m,ido,2,k)
            ch(m,1,k,2) = cc(m,1,1,k) - cc(m,ido,2,k)
        end do
    end do

    if (ido < 2) return
    if (ido == 2) then
        do k = 1, l1
            do m = 1, mp
                ch(m,ido,k,1) = cc(m,ido,1,k) + cc(m,ido,1,k)
                ch(m,ido,k,2) = -(cc(m,1,2,k) + cc(m,1,2,k))
            end do
        end do
        return
    end if

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                ch(m,i-1,k,1) = cc(m,i-1,1,k) + cc(m,ic-1,2,k)
                ch(m,i,k,1) = cc(m,i,1,k) - cc(m,ic,2,k)
                tr2 = cc(m,i-1,1,k) - cc(m,ic-1,2,k)
                ti2 = cc(m,i,1,k) + cc(m,ic,2,k)
                ch(m,i-1,k,2) = wa1(i-2)*tr2 - wa1(i-1)*ti2
                ch(m,i,k,2) = wa1(i-2)*ti2 + wa1(i-1)*tr2
            end do
        end do
    end do
end subroutine hradb2

subroutine hradb3(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,3,l1)
    real, intent(out) :: ch(mdimch,ido,l1,3)
    real, intent(in) :: wa1(ido), wa2(ido)

    integer :: k, m, i, ic, idp2
    real, parameter :: taur = -0.5
    real, parameter :: taui = 0.866025403784439
    real :: tr2, ti2, tr3, ti3, cr2, ci2, cr3, ci3

    do k = 1, l1
        do m = 1, mp
            tr2 = cc(m,ido,2,k) + cc(m,ido,2,k)
            ch(m,1,k,1) = cc(m,1,1,k) + tr2
            ch(m,1,k,2) = cc(m,1,1,k) + taur*tr2 - taui*(cc(m,1,3,k) + cc(m,1,3,k))
            ch(m,1,k,3) = cc(m,1,1,k) + taur*tr2 + taui*(cc(m,1,3,k) + cc(m,1,3,k))
        end do
    end do

    if (ido == 1) return

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                tr2 = cc(m,i-1,3,k) + cc(m,ic-1,2,k)
                ti2 = cc(m,i,3,k) - cc(m,ic,2,k)
                ch(m,i-1,k,1) = cc(m,i-1,1,k) + tr2
                ch(m,i,k,1) = cc(m,i,1,k) + ti2
                cr2 = cc(m,i-1,1,k) + taur*tr2
                ci2 = cc(m,i,1,k) + taur*ti2
                cr3 = taui*(cc(m,i,3,k) + cc(m,ic,2,k))
                ci3 = taui*(cc(m,i-1,3,k) - cc(m,ic-1,2,k))
                ch(m,i-1,k,2) = wa1(i-2)*(cr2-ci3) - wa1(i-1)*(ci2+cr3)
                ch(m,i,k,2) = wa1(i-2)*(ci2+cr3) + wa1(i-1)*(cr2-ci3)
                ch(m,i-1,k,3) = wa2(i-2)*(cr2+ci3) - wa2(i-1)*(ci2-cr3)
                ch(m,i,k,3) = wa2(i-2)*(ci2-cr3) + wa2(i-1)*(cr2+ci3)
            end do
        end do
    end do
end subroutine hradb3

subroutine hradb5(mp,ido,l1,cc,mdimcc,ch,mdimch,wa1,wa2,wa3,wa4)
    implicit none
    integer, intent(in) :: mp, ido, l1, mdimcc, mdimch
    real, intent(in) :: cc(mdimcc,ido,5,l1)
    real, intent(out) :: ch(mdimch,ido,l1,5)
    real, intent(in) :: wa1(ido), wa2(ido), wa3(ido), wa4(ido)

    integer :: k, m, i, ic, idp2
    real, parameter :: tr11 = 0.309016994374947
    real, parameter :: ti11 = 0.951056516295154
    real, parameter :: tr12 = -0.809016994374947
    real, parameter :: ti12 = 0.587785252292473
    real :: tr2, ti2, tr3, ti3, tr4, ti4, tr5, ti5
    real :: cr2, ci2, cr3, ci3, cr4, ci4, cr5, ci5

    do k = 1, l1
        do m = 1, mp
            ti5 = cc(m,1,3,k) + cc(m,1,3,k)
            ti4 = cc(m,1,5,k) + cc(m,1,5,k)
            tr2 = cc(m,ido,2,k) + cc(m,ido,2,k)
            tr3 = cc(m,ido,4,k) + cc(m,ido,4,k)
            ch(m,1,k,1) = cc(m,1,1,k) + tr2 + tr3
            ch(m,1,k,2) = (cc(m,1,1,k) + tr11*tr2 + tr12*tr3) - (ti11*ti5 + ti12*ti4)
            ch(m,1,k,3) = (cc(m,1,1,k) + tr12*tr2 + tr11*tr3) - (ti12*ti5 - ti11*ti4)
            ch(m,1,k,4) = (cc(m,1,1,k) + tr12*tr2 + tr11*tr3) + (ti12*ti5 - ti11*ti4)
            ch(m,1,k,5) = (cc(m,1,1,k) + tr11*tr2 + tr12*tr3) + (ti11*ti5 + ti12*ti4)
        end do
    end do

    if (ido == 1) return

    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            do m = 1, mp
                tr2 = cc(m,i-1,3,k) + cc(m,ic-1,2,k)
                ti2 = cc(m,i,3,k) - cc(m,ic,2,k)
                tr3 = cc(m,i-1,5,k) + cc(m,ic-1,4,k)
                ti3 = cc(m,i,5,k) - cc(m,ic,4,k)
                ch(m,i-1,k,1) = cc(m,i-1,1,k) + tr2 + tr3
                ch(m,i,k,1) = cc(m,i,1,k) + ti2 + ti3
                cr2 = cc(m,i-1,1,k) + tr11*tr2 + tr12*tr3
                ci2 = cc(m,i,1,k) + tr11*ti2 + tr12*ti3
                cr3 = cc(m,i-1,1,k) + tr12*tr2 + tr11*tr3
                ci3 = cc(m,i,1,k) + tr12*ti2 + tr11*ti3
                tr4 = ti11*(cc(m,i,3,k) + cc(m,ic,2,k)) + ti12*(cc(m,i,5,k) + cc(m,ic,4,k))
                ti4 = ti11*(cc(m,i-1,3,k) - cc(m,ic-1,2,k)) + ti12*(cc(m,i-1,5,k) - cc(m,ic-1,4,k))
                tr5 = ti12*(cc(m,i,3,k) + cc(m,ic,2,k)) - ti11*(cc(m,i,5,k) + cc(m,ic,4,k))
                ti5 = ti12*(cc(m,i-1,3,k) - cc(m,ic-1,2,k)) - ti11*(cc(m,i-1,5,k) - cc(m,ic-1,4,k))
                ch(m,i-1,k,2) = wa1(i-2)*(cr2-tr4) - wa1(i-1)*(ci2+ti4)
                ch(m,i,k,2) = wa1(i-2)*(ci2+ti4) + wa1(i-1)*(cr2-tr4)
                ch(m,i-1,k,3) = wa2(i-2)*(cr3-tr5) - wa2(i-1)*(ci3+ti5)
                ch(m,i,k,3) = wa2(i-2)*(ci3+ti5) + wa2(i-1)*(cr3-tr5)
                ch(m,i-1,k,4) = wa3(i-2)*(cr3+tr5) - wa3(i-1)*(ci3-ti5)
                ch(m,i,k,4) = wa3(i-2)*(ci3-ti5) + wa3(i-1)*(cr3+tr5)
                ch(m,i-1,k,5) = wa4(i-2)*(cr2+tr4) - wa4(i-1)*(ci2-ti4)
                ch(m,i,k,5) = wa4(i-2)*(ci2-ti4) + wa4(i-1)*(cr2+tr4)
            end do
        end do
    end do
end subroutine hradb5

subroutine hradbg(mp,ido,ip,l1,idl1,cc,c1,c2,mdimcc,ch,ch2,mdimch,wa)
    implicit none
    integer, intent(in) :: mp, ido, ip, l1, idl1, mdimcc, mdimch
    real, intent(inout) :: cc(mdimcc,ido,ip,l1)
    real, intent(inout) :: c1(mdimcc,ido,l1,ip)
    real, intent(inout) :: c2(mdimcc,idl1,ip)
    real, intent(inout) :: ch(mdimch,ido,l1,ip)
    real, intent(inout) :: ch2(mdimch,idl1,ip)
    real, intent(in) :: wa(ido)

    integer :: idp2, nbd, ipp2, ipph, is, j, k, i, ik, m, idij, jc, l, lc, j2, ic
    real :: dcp, dsp, arg, ar1, ai1, ar1h, ar2, ai2, ar2h, dc2, ds2

    arg = 8.0*atan(1.0) / real(ip)
    dcp = cos(arg)
    dsp = sin(arg)
    idp2 = ido + 2
    nbd = (ido - 1) / 2
    ipp2 = ip + 2
    ipph = (ip + 1) / 2

    if (ido >= l1) then
        do k = 1, l1
            do i = 1, ido
                do m = 1, mp
                    ch(m,i,k,1) = cc(m,i,1,k)
                end do
            end do
        end do
    else
        do i = 1, ido
            do k = 1, l1
                do m = 1, mp
                    ch(m,i,k,1) = cc(m,i,1,k)
                end do
            end do
        end do
    end if

    do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
            do m = 1, mp
                ch(m,1,k,j) = cc(m,ido,2*j-2,k) + cc(m,ido,2*j-2,k)
                ch(m,1,k,jc) = cc(m,1,2*j-1,k) + cc(m,1,2*j-1,k)
            end do
        end do
    end do

    if (ido /= 1) then
        if (nbd >= l1) then
            do j = 2, ipph
                jc = ipp2 - j
                do k = 1, l1
                    do i = 3, ido, 2
                        ic = idp2 - i
                        do m = 1, mp
                            ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                            ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) - cc(m,ic-1,2*j-2,k)
                            ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                            ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
                        end do
                    end do
                end do
            end do
        else
            do j = 2, ipph
                jc = ipp2 - j
                do i = 3, ido, 2
                    ic = idp2 - i
                    do k = 1, l1
                        do m = 1, mp
                            ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                            ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) - cc(m,ic-1,2*j-2,k)
                            ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                            ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
                        end do
                    end do
                end do
            end do
        end if
    end if

    ar1 = 1.0
    ai1 = 0.0

    do l = 2, ipph
        lc = ipp2 - l
        ar1h = dcp*ar1 - dsp*ai1
        ai1 = dcp*ai1 + dsp*ar1
        ar1 = ar1h

        do ik = 1, idl1
            do m = 1, mp
                c2(m,ik,l) = ch2(m,ik,1) + ar1*ch2(m,ik,2)
                c2(m,ik,lc) = ai1*ch2(m,ik,ip)
            end do
        end do

        dc2 = ar1
        ds2 = ai1
        ar2 = ar1
        ai2 = ai1

        do j = 3, ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h

            do ik = 1, idl1
                do m = 1, mp
                    c2(m,ik,l) = c2(m,ik,l) + ar2*ch2(m,ik,j)
                    c2(m,ik,lc) = c2(m,ik,lc) + ai2*ch2(m,ik,jc)
                end do
            end do
        end do
    end do

    do j = 2, ipph
        do ik = 1, idl1
            do m = 1, mp
                ch2(m,ik,1) = ch2(m,ik,1) + ch2(m,ik,j)
            end do
        end do
    end do

    do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
            do m = 1, mp
                ch(m,1,k,j) = c1(m,1,k,j) - c1(m,1,k,jc)
                ch(m,1,k,jc) = c1(m,1,k,j) + c1(m,1,k,jc)
            end do
        end do
    end do

    if (ido /= 1) then
        if (nbd >= l1) then
            do j = 2, ipph
                jc = ipp2 - j
                do k = 1, l1
                    do i = 3, ido, 2
                        do m = 1, mp
                            ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                            ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                            ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                            ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
                        end do
                    end do
                end do
            end do
        else
            do j = 2, ipph
                jc = ipp2 - j
                do i = 3, ido, 2
                    do k = 1, l1
                        do m = 1, mp
                            ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                            ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                            ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                            ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
                        end do
                    end do
                end do
            end do
        end if
    end if

    if (ido == 1) return

    do ik = 1, idl1
        do m = 1, mp
            c2(m,ik,1) = ch2(m,ik,1)
        end do
    end do

    do j = 2, ip
        do k = 1, l1
            do m = 1, mp
                c1(m,1,k,j) = ch(m,1,k,j)
            end do
        end do
    end do

    if (nbd <= l1) then
        is = -ido
        do j = 2, ip
            is = is + ido
            idij = is
            do i = 3, ido, 2
                idij = idij + 2
                do k = 1, l1
                    do m = 1, mp
                        c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) - wa(idij)*ch(m,i,k,j)
                        c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) + wa(idij)*ch(m,i-1,k,j)
                    end do
                end do
            end do
        end do
    else
        is = -ido
        do j = 2, ip
            is = is + ido
            do k = 1, l1
                idij = is
                do i = 3, ido, 2
                    idij = idij + 2
                    do m = 1, mp
                        c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) - wa(idij)*ch(m,i,k,j)
                        c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) + wa(idij)*ch(m,i-1,k,j)
                    end do
                end do
            end do
        end do
    end if
end subroutine hradbg

end module hrfft_mod
