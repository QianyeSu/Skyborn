! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-17
! File    : grid2triple.f90
! Purpose : Pack a rectilinear grid into the legacy triple-array layout used
!           by the interp Fortran entry points.
! Notes   : This free-form implementation preserves the historical exact
!           missing-value comparison semantics from the archived fixed-form
!           routine.
! =============================================================================
!
subroutine grid2triple(x, y, z, mx, ny, d, ldmax, ld, zmsg, ier)
    implicit none

    ! QUICK REFERENCE
    ! PURPOSE
    !    PACK A 2-D RECTILINEAR GRID INTO A 3 X LD "TRIPLE" ARRAY
    !    CONTAINING X, Y, AND Z VALUES FOR EACH NON-MISSING GRID CELL.
    !
    ! INPUTS
    !    X(MX)      - X COORDINATE VALUES FOR THE RIGHTMOST GRID DIMENSION
    !    Y(NY)      - Y COORDINATE VALUES FOR THE LEFTMOST GRID DIMENSION
    !    Z(MX, NY)  - INPUT GRID VALUES
    !    LDMAX      - MAXIMUM AVAILABLE ROW COUNT IN D; THE F2PY WRAPPER
    !                 PASSES MX * NY SO THE OUTPUT BUFFER CAN HOLD ALL CELLS
    !    ZMSG       - MISSING VALUE SENTINEL. CELLS EQUAL TO ZMSG ARE SKIPPED
    !
    ! OUTPUTS
    !    D(LDMAX, 3)
    !               COLUMN 1 STORES X, COLUMN 2 STORES Y, COLUMN 3 STORES Z
    !               FOR EACH NON-MISSING INPUT CELL
    !    LD         - NUMBER OF VALID TRIPLES WRITTEN TO D
    !    IER        - 0 ON SUCCESS, -10 WHEN ALL INPUT CELLS ARE MISSING

    integer, intent(in) :: mx, ny, ldmax
    integer, intent(out) :: ld, ier
    integer :: m, n
    real(8), intent(in) :: x(mx), y(ny), z(mx, ny), zmsg
    real(8), intent(out) :: d(ldmax, 3)
    real(8) :: y_n

    ier = 0
    ld = 0

    do n = 1, ny
        y_n = y(n)
        do m = 1, mx
            ! Preserve the historical exact sentinel comparison so missing-cell
            ! detection stays bit-for-bit compatible with the legacy routine.
            if (z(m, n) /= zmsg) then
                ld = ld + 1
                d(ld, 1) = x(m)
                d(ld, 2) = y_n
                d(ld, 3) = z(m, n)
            end if
        end do
    end do

    if (ld == 0) then
        ier = -10
    end if
end subroutine grid2triple
