subroutine qzval(nm, n, a, b, alfr, alfi, beta)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    !
    integer i,j,n,en,na,nm,nn,isw
    real(real64) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n)
    real(real64) c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1, &
    & u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr, &
    & ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb
    !
    !     this subroutine is the third step of the qz algorithm
    !     for solving generalized matrix eigenvalue problems,
    !     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
    !
    !     this subroutine accepts a pair of real matrices, one of them
    !     in quasi-triangular form and the other in upper triangular form.
    !     it reduces the quasi-triangular matrix further, so that any
    !     remaining 2-by-2 blocks correspond to pairs of complex
    !     eigenvalues, and returns quantities whose ratios give the
    !     generalized eigenvalues.  it is usually preceded by  qzhes
    !     and  qzit  and may be followed by  qzvec.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrices.
    !
    !        a contains a real upper quasi-triangular matrix.
    !
    !        b contains a real upper triangular matrix.  in addition,
    !          location b(n,1) contains the tolerance quantity (epsb)
    !          computed and saved in  qzit.
    !
    !     on output
    !
    !        a has been reduced further to a quasi-triangular matrix
    !          in which all nonzero subdiagonal elements correspond to
    !          pairs of complex eigenvalues.
    !
    !        b is still in upper triangular form, although its elements
    !          have been altered.  b(n,1) is unaltered.
    !
    !        alfr and alfi contain the real and imaginary parts of the
    !          diagonal elements of the triangular matrix that would be
    !          obtained if a were reduced completely to triangular form
    !          by unitary transformations.  non-zero values of alfi occur
    !          in pairs, the first member positive and the second negative.
    !
    !        beta contains the diagonal elements of the corresponding b,
    !          normalized to be real and non-negative.  the generalized
    !          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    epsb = b(n,1)
    isw = 1
    !     .......... find eigenvalues of quasi-triangular matrices.
    !                for en=n step -1 until 1 do -- ..........
    do nn = 1, n
        en = n + 1 - nn
        na = en - 1
        if (isw .eq. 2) go to 505
        if (en .eq. 1) go to 410
        if (a(en,na) .ne. 0.0d0) go to 420
        !     .......... 1-by-1 block, one real root ..........
        410 continue
        alfr(en) = a(en,en)
        if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
        beta(en) = abs(b(en,en))
        alfi(en) = 0.0d0
        go to 510
        !     .......... 2-by-2 block ..........
        420 continue
        if (abs(b(na,na)) .le. epsb) go to 455
        if (abs(b(en,en)) .gt. epsb) go to 430
        a1 = a(en,en)
        a2 = a(en,na)
        bn = 0.0d0
        go to 435
        430 continue
        an = abs(a(na,na)) + abs(a(na,en)) + abs(a(en,na)) &
        & + abs(a(en,en))
        bn = abs(b(na,na)) + abs(b(na,en)) + abs(b(en,en))
        a11 = a(na,na) / an
        a12 = a(na,en) / an
        a21 = a(en,na) / an
        a22 = a(en,en) / an
        b11 = b(na,na) / bn
        b12 = b(na,en) / bn
        b22 = b(en,en) / bn
        e = a11 / b11
        ei = a22 / b22
        s = a21 / (b11 * b22)
        t = (a22 - e * b22) / b22
        if (abs(e) .le. abs(ei)) go to 431
        e = ei
        t = (a11 - e * b11) / b11
        431 continue
        c = 0.5d0 * (t - s * b12)
        d = c * c + s * (a12 - e * b12)
        if (d .lt. 0.0d0) go to 480
        !     .......... two real roots.
        !                zero both a(en,na) and b(en,na) ..........
        e = e + (c + sign(sqrt(d),c))
        a11 = a11 - e * b11
        a12 = a12 - e * b12
        a22 = a22 - e * b22
        if (abs(a11) + abs(a12) .lt. &
        & abs(a21) + abs(a22)) go to 432
        a1 = a12
        a2 = a11
        go to 435
        432 continue
        a1 = a22
        a2 = a21
        !     .......... choose and apply real z ..........
        435 continue
        s = abs(a1) + abs(a2)
        u1 = a1 / s
        u2 = a2 / s
        r = sign(sqrt(u1*u1+u2*u2),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1
        !
        do i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
        end do
        !
        450 continue
        if (bn .eq. 0.0d0) go to 475
        if (an .lt. abs(e) * bn) go to 455
        a1 = b(na,na)
        a2 = b(en,na)
        go to 460
        455 continue
        a1 = a(na,na)
        a2 = a(en,na)
        !     .......... choose and apply real q ..........
        460 continue
        s = abs(a1) + abs(a2)
        if (s .eq. 0.0d0) go to 475
        u1 = a1 / s
        u2 = a2 / s
        r = sign(sqrt(u1*u1+u2*u2),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1
        !
        do j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
        end do
        !
        475 continue
        a(en,na) = 0.0d0
        b(en,na) = 0.0d0
        alfr(na) = a(na,na)
        alfr(en) = a(en,en)
        if (b(na,na) .lt. 0.0d0) alfr(na) = -alfr(na)
        if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
        beta(na) = abs(b(na,na))
        beta(en) = abs(b(en,en))
        alfi(en) = 0.0d0
        alfi(na) = 0.0d0
        go to 505
        !     .......... two complex roots ..........
        480 continue
        e = e + c
        ei = sqrt(-d)
        a11r = a11 - e * b11
        a11i = ei * b11
        a12r = a12 - e * b12
        a12i = ei * b12
        a22r = a22 - e * b22
        a22i = ei * b22
        if (abs(a11r) + abs(a11i) + abs(a12r) + abs(a12i) .lt. &
        & abs(a21) + abs(a22r) + abs(a22i)) go to 482
        a1 = a12r
        a1i = a12i
        a2 = -a11r
        a2i = -a11i
        go to 485
        482 continue
        a1 = a22r
        a1i = a22i
        a2 = -a21
        a2i = 0.0d0
        !     .......... choose complex z ..........
        485 continue
        cz = sqrt(a1*a1+a1i*a1i)
        if (cz .eq. 0.0d0) go to 487
        szr = (a1 * a2 + a1i * a2i) / cz
        szi = (a1 * a2i - a1i * a2) / cz
        r = sqrt(cz*cz+szr*szr+szi*szi)
        cz = cz / r
        szr = szr / r
        szi = szi / r
        go to 490
        487 continue
        szr = 1.0d0
        szi = 0.0d0
        490 continue
        if (an .lt. (abs(e) + ei) * bn) go to 492
        a1 = cz * b11 + szr * b12
        a1i = szi * b12
        a2 = szr * b22
        a2i = szi * b22
        go to 495
        492 continue
        a1 = cz * a11 + szr * a12
        a1i = szi * a12
        a2 = cz * a21 + szr * a22
        a2i = szi * a22
        !     .......... choose complex q ..........
        495 continue
        cq = sqrt(a1*a1+a1i*a1i)
        if (cq .eq. 0.0d0) go to 497
        sqr = (a1 * a2 + a1i * a2i) / cq
        sqi = (a1 * a2i - a1i * a2) / cq
        r = sqrt(cq*cq+sqr*sqr+sqi*sqi)
        cq = cq / r
        sqr = sqr / r
        sqi = sqi / r
        go to 500
        497 continue
        sqr = 1.0d0
        sqi = 0.0d0
        !     .......... compute diagonal elements that would result
        !                if transformations were applied ..........
        500 continue
        ssr = sqr * szr + sqi * szi
        ssi = sqr * szi - sqi * szr
        i = 1
        tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 &
        & + ssr * a22
        ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
        dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
        di = cq * szi * b12 + ssi * b22
        go to 503
        502 continue
        i = 2
        tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 &
        & + cq * cz * a22
        ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
        dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
        di = -ssi * b11 - sqi * cz * b12
        503 continue
        t = ti * dr - tr * di
        j = na
        if (t .lt. 0.0d0) j = en
        r = sqrt(dr*dr+di*di)
        beta(j) = bn * r
        alfr(j) = an * (tr * dr + ti * di) / r
        alfi(j) = an * t / r
        if (i .eq. 1) go to 502
        505 continue
        isw = 3 - isw
        510 continue
    end do
    b(n,1) = epsb
    !
    return
end subroutine qzval
