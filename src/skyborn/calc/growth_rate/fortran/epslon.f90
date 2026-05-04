real(real64) function epslon(x)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none

    real(real64), intent(in) :: x
    real(real64) :: a, b, c, eps

    !     estimate unit roundoff in quantities of size x.
    !
    !     this program should function properly on all systems
    !     satisfying the following two assumptions,
    !        1.  the base used in representing floating point
    !            numbers is not a power of three.
    !        2.  the quantity  a  in statement 10 is represented to
    !            the accuracy used in floating point variables
    !            that are stored in memory.
    !     the statement number 10 and the go to 10 are intended to
    !     force optimizing compilers to generate code satisfying
    !     assumption 2.
    !     under these assumptions, it should be true that,
    !            a  is not exactly equal to four-thirds,
    !            b  has a zero for its last bit or digit,
    !            c  is not exactly equal to one,
    !            eps  measures the separation of 1.0 from
    !                 the next larger floating point number.
    !     the developers of eispack would appreciate being informed
    !     about any systems where these assumptions do not hold.
    !
    !     this version dated 4/6/83.

    a = 4.0_real64 / 3.0_real64
    do
        b = a - 1.0_real64
        c = b + b + b
        eps = abs(c - 1.0_real64)
        if (eps /= 0.0_real64) exit
    end do

    epslon = eps * abs(x)
end function epslon
