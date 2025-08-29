! SIMD optimized Fortran version of geostrophic wind calculation
! Using OpenMP SIMD directives for vectorization
!
! NCL: uv = z2geouv(z,lat,lon,iopt)           ; uv(2,...,nlat,mlon)
!
! Given geopotential height (gpm): calculate geostrophic wind components
!
! INPUT
! .   z     - geopotential height (south to north order)
! .   mlon  - number of longitudes
! .   nlat  - number of latitudes
! .   zmsg  - missing value code
! .   glat  - latitudes  of grid points
! .   glon  - longitudes of grid points
! .   iopt  - option: iopt=0 z is not cyclic in longitude
! .                   iopt=1 z is cyclic in longitude
! OUTPUT
! .   ug,vg - geostrophic wind components

subroutine z2geouv(z, mlon, nlat, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! Input parameters
    integer, intent(in) :: mlon, nlat, iopt
    real(kind=8), intent(in) :: z(mlon,nlat), zmsg
    real(kind=8), intent(in) :: glat(nlat), glon(mlon)

    ! Output parameters
    real(kind=8), intent(out) :: ug(mlon,nlat), vg(mlon,nlat)

    if (glat(2) > glat(1)) then
        ! South to north ordering
        call z2guv(z, mlon, nlat, zmsg, glon, glat, ug, vg, iopt)
    else
        ! North to south ordering - need to reorder
        call zuvnew(z, mlon, nlat, zmsg, glon, glat, ug, vg, iopt)
    end if

end subroutine z2geouv

! Subroutine to handle north-to-south data by reordering
!
! Data array needs to be rearranged to go from south to north
! before calling the main calculation routine
subroutine zuvnew(z, mlon, nlat, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! Input parameters
    integer, intent(in) :: mlon, nlat, iopt
    real(kind=8), intent(in) :: z(mlon,nlat), zmsg
    real(kind=8), intent(in) :: glat(nlat), glon(mlon)

    ! Output parameters
    real(kind=8), intent(out) :: ug(mlon,nlat), vg(mlon,nlat)

    ! Local temporary arrays
    real(kind=8) :: ztmp(mlon,nlat), glattmp(nlat)
    real(kind=8) :: utmp(nlat), vtmp(nlat)
    integer :: nl, ml, nlat1

    nlat1 = nlat + 1

    ! Reorder from north-south to south-north [S ==> N] with SIMD
    do nl = 1, nlat
        glattmp(nlat1-nl) = glat(nl)
        !$omp simd
        do ml = 1, mlon
            ztmp(ml, nlat1-nl) = z(ml, nl)
        end do
        !$omp end simd
    end do

    ! Calculate geostrophic winds using reordered arrays
    call z2guv(ztmp, mlon, nlat, zmsg, glon, glattmp, ug, vg, iopt)

    ! Put results back in original order with SIMD
    do ml = 1, mlon
        !$omp simd
        do nl = 1, nlat
            utmp(nlat1-nl) = ug(ml, nl)
            vtmp(nlat1-nl) = vg(ml, nl)
        end do
        !$omp end simd
        !$omp simd
        do nl = 1, nlat
            ug(ml, nl) = utmp(nl)
            vg(ml, nl) = vtmp(nl)
        end do
        !$omp end simd
    end do

end subroutine zuvnew

! Main calculation subroutine for geostrophic wind - SIMD OPTIMIZED
!
! Given geopotential height (gpm): calculate geostrophic components
!
! INPUT
! .   z     - geopotential height (south to north order)
! .   mlon  - number of longitudes
! .   nlat  - number of latitudes
! .   zmsg  - missing value code
! .   glat  - latitudes  of grid points
! .   glon  - longitudes of grid points
! .   iopt  - option: iopt=0 z is not cyclic in longitude
! .                   iopt=1 z is cyclic in longitude
! OUTPUT
! .   ug,vg - geostrophic wind components
subroutine z2guv(z, mlon, nlat, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! Input parameters
    integer, intent(in) :: mlon, nlat, iopt
    real(kind=8), intent(in) :: z(mlon,nlat), zmsg
    real(kind=8), intent(in) :: glat(nlat), glon(mlon)

    ! Output parameters
    real(kind=8), intent(out) :: ug(mlon,nlat), vg(mlon,nlat)

    ! Local variables
    real(kind=8) :: zz(0:mlon+1, 0:nlat+1)
    real(kind=8) :: gf(nlat), dx(nlat), dy(nlat)
    real(kind=8) :: inv_gf(nlat), inv_dx(nlat), inv_dy(nlat)
    real(kind=8), parameter :: gr = 9.80616_8
    real(kind=8), parameter :: omega = 7.292e-5_8
    real(kind=8), parameter :: re = 6371220.0_8
    real(kind=8), parameter :: pi = 4.0_8 * atan(1.0_8)
    real(kind=8), parameter :: rad = pi / 180.0_8
    real(kind=8), parameter :: two_omega = 2.0_8 * omega
    real(kind=8), parameter :: gr_over_two_omega = gr / two_omega
    real(kind=8) :: rr, dlon, dlon2
    real(kind=8) :: sin_lat, cos_lat, abs_lat
    integer :: nl, ml
    logical :: valid_lat(nlat)

    ! Calculate 'constant' quantities
    rr = re * rad
    dlon = rr * (glon(2) - glon(1))
    dlon2 = 2.0_8 * dlon

    ! Calculate the quantity [gravity/coriolis_force] and precompute inverses - SIMD optimized
    ! dx and dy represent the distances over which finite differences apply
    ! Note: dy(1) and dy(nlat) are half width for one sided differences

    !$omp simd
    do nl = 1, nlat
        abs_lat = abs(glat(nl))
        valid_lat(nl) = (glat(nl) /= 0.0_8) .and. (abs_lat <= 89.999_8)

        if (valid_lat(nl)) then
            sin_lat = sin(glat(nl) * rad)
            cos_lat = cos(glat(nl) * rad)
            gf(nl) = gr_over_two_omega / sin_lat
            dx(nl) = dlon2 * cos_lat

            ! Precompute inverses to avoid division in inner loop
            inv_gf(nl) = 1.0_8 / gf(nl)
            if (dx(nl) /= 0.0_8) then
                inv_dx(nl) = 1.0_8 / dx(nl)
            else
                inv_dx(nl) = 0.0_8
            end if
        else
            gf(nl) = zmsg
            dx(nl) = 0.0_8
            inv_gf(nl) = 0.0_8
            inv_dx(nl) = 0.0_8
        end if
    end do
    !$omp end simd

    ! Get around 'rounding' problem at poles
    if (glat(1) < -89.999_8) then
        dx(1) = 0.0_8
        inv_dx(1) = 0.0_8
    end if
    if (glat(nlat) > 89.999_8) then
        dx(nlat) = 0.0_8
        inv_dx(nlat) = 0.0_8
    end if

    ! Latitude can have variable spacing - SIMD optimized
    ! dy(1) and dy(nlat) are half size for use in one-sided differences
    dy(1) = (glat(2) - glat(1)) * rr
    inv_dy(1) = 1.0_8 / dy(1)

    !$omp simd
    do nl = 2, nlat - 1
        dy(nl) = (glat(nl+1) - glat(nl-1)) * rr
        inv_dy(nl) = 1.0_8 / dy(nl)
    end do
    !$omp end simd

    dy(nlat) = (glat(nlat) - glat(nlat-1)) * rr
    inv_dy(nlat) = 1.0_8 / dy(nlat)

    ! Set ug and vg to default values - vectorized
    ug = zmsg
    vg = zmsg

    ! Use zz array to eliminate the need for special
    ! loops/statements to handle the east and west boundaries

    ! Fill zz 'middle' - SIMD optimized
    do nl = 1, nlat
        !$omp simd
        do ml = 1, mlon
            zz(ml, nl) = z(ml, nl)
        end do
        !$omp end simd
    end do

    ! Bottom and top boundaries - SIMD optimized
    !$omp simd
    do ml = 1, mlon
        zz(ml, 0) = zz(ml, 1)
        zz(ml, nlat+1) = zz(ml, nlat)
    end do
    !$omp end simd

    ! Left and right boundaries - SIMD optimized
    if (iopt == 1) then
        ! Specify cyclic values
        !$omp simd
        do nl = 0, nlat+1
            zz(0, nl) = zz(mlon, nl)
            zz(mlon+1, nl) = zz(1, nl)
        end do
        !$omp end simd
    else
        ! Duplicate adjacent points
        !$omp simd
        do nl = 0, nlat+1
            zz(0, nl) = zz(1, nl)
            zz(mlon+1, nl) = zz(mlon, nl)
        end do
        !$omp end simd
    end if

    ! Calculate the winds at all grid points of array zz - HIGHLY SIMD OPTIMIZED
    ! ug calculated via one-sided differences [crude approx]
    ! Note: a 'trick' is used
    ! (1) the top and bottom rows of zz are duplicates of adjacent rows
    ! (2) dy(1) and dy(nlat) are already set for one-sided differences
    ! (3) the combination will result in numerically correct results

    ! Main computation loop - SIMD optimized with precomputed values
    do nl = 1, nlat
        if (valid_lat(nl)) then
            !$omp simd
            do ml = 1, mlon
                ! Check for missing values and compute both components
                if (zz(ml, nl+1) /= zmsg .and. zz(ml, nl-1) /= zmsg .and. &
                    zz(ml+1, nl) /= zmsg .and. zz(ml-1, nl) /= zmsg) then

                    ! Use precomputed inverse values to avoid division
                    ug(ml, nl) = -gf(nl) * (zz(ml, nl+1) - zz(ml, nl-1)) * inv_dy(nl)
                    vg(ml, nl) = gf(nl) * (zz(ml+1, nl) - zz(ml-1, nl)) * inv_dx(nl)
                else
                    ! Handle partial missing data
                    if (zz(ml, nl+1) /= zmsg .and. zz(ml, nl-1) /= zmsg) then
                        ug(ml, nl) = -gf(nl) * (zz(ml, nl+1) - zz(ml, nl-1)) * inv_dy(nl)
                    end if
                    if (zz(ml+1, nl) /= zmsg .and. zz(ml-1, nl) /= zmsg) then
                        vg(ml, nl) = gf(nl) * (zz(ml+1, nl) - zz(ml-1, nl)) * inv_dx(nl)
                    end if
                end if
            end do
            !$omp end simd
        end if
    end do

    ! If not cyclic in longitude then "adjust" the vg for the left/right columns - SIMD optimized
    ! Basically, they should be one sided differences BUT they were
    ! divided by 2*dx. The multiplication by 2 below undoes this.
    ! The reason for the "if" is to handle the case at the EQ

    if (iopt /= 1) then
        !$omp simd
        do nl = 1, nlat
            if (vg(1, nl) /= zmsg) vg(1, nl) = vg(1, nl) * 2.0_8
            if (vg(mlon, nl) /= zmsg) vg(mlon, nl) = vg(mlon, nl) * 2.0_8
        end do
        !$omp end simd
    end if

end subroutine z2guv
