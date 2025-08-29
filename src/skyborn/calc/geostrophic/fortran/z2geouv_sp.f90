! Multi-dimensional SIMD optimized geostrophic wind calculation - SINGLE PRECISION
! Supporting 2D, 3D, and 4D input data with shape (nlat, mlon, ...)
! Direct processing without unnecessary transposes
!
! NCL equivalent: uv = z2geouv(z,lat,lon,iopt)
!
! Given geopotential height (gpm): calculate geostrophic wind components
! Uses finite differences and geostrophic approximation:
!   ug = -(g/f) * dZ/dy   (zonal component)
!   vg =  (g/f) * dZ/dx   (meridional component)
! where g = gravity, f = Coriolis parameter, Z = geopotential height

!==============================================================================
! 2D version - shape (nlat, mlon) - SINGLE PRECISION
!==============================================================================
subroutine z2geouv(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! INPUT - single precision
    integer, intent(in) :: nlat        ! number of latitudes (south to north)
    integer, intent(in) :: mlon        ! number of longitudes
    integer, intent(in) :: iopt        ! longitude boundary option:
                                       ! iopt=0: z is not cyclic in longitude
                                       ! iopt=1: z is cyclic in longitude
    real(kind=4), intent(in) :: z(nlat, mlon)  ! geopotential height [gpm] (nlat x mlon)
    real(kind=4), intent(in) :: zmsg           ! missing value code
    real(kind=4), intent(in) :: glat(nlat)     ! latitudes of grid points [degrees] (south to north)
    real(kind=4), intent(in) :: glon(mlon)     ! longitudes of grid points [degrees]

    ! OUTPUT - single precision
    real(kind=4), intent(out) :: ug(nlat, mlon) ! zonal geostrophic wind [m/s] (nlat x mlon)
    real(kind=4), intent(out) :: vg(nlat, mlon) ! meridional geostrophic wind [m/s] (nlat x mlon)

    if (glat(2) > glat(1)) then
        ! South to north ordering - direct calculation
        call z2guv(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    else
        ! North to south ordering - need to reorder first
        call zuvnew(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    end if

end subroutine z2geouv

!==============================================================================
! 3D version - shape (nlat, mlon, n3rd) - SINGLE PRECISION
! Processes multiple 2D slices (e.g., pressure levels, time steps)
!==============================================================================
subroutine z2geouv_3d(z, nlat, mlon, n3rd, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! INPUT - single precision
    integer, intent(in) :: nlat               ! number of latitudes
    integer, intent(in) :: mlon               ! number of longitudes
    integer, intent(in) :: n3rd               ! size of 3rd dimension (e.g., pressure levels)
    integer, intent(in) :: iopt               ! longitude boundary option (0=non-cyclic, 1=cyclic)
    real(kind=4), intent(in) :: z(nlat, mlon, n3rd)    ! geopotential height [gpm]
    real(kind=4), intent(in) :: zmsg                   ! missing value code
    real(kind=4), intent(in) :: glat(nlat)             ! latitudes [degrees] (south to north)
    real(kind=4), intent(in) :: glon(mlon)             ! longitudes [degrees]

    ! OUTPUT - single precision
    real(kind=4), intent(out) :: ug(nlat, mlon, n3rd)  ! zonal geostrophic wind [m/s]
    real(kind=4), intent(out) :: vg(nlat, mlon, n3rd)  ! meridional geostrophic wind [m/s]

    ! Local variables
    integer :: k

    ! Process each 2D slice along the 3rd dimension
    ! OpenMP parallelization for multiple levels/times
    !$omp parallel do private(k) if(n3rd > 4)
    do k = 1, n3rd
        call z2geouv(z(:, :, k), nlat, mlon, zmsg, glon, glat, &
                     ug(:, :, k), vg(:, :, k), iopt)
    end do
    !$omp end parallel do

end subroutine z2geouv_3d

!==============================================================================
! 4D version - shape (nlat, mlon, n3rd, n4th) - SINGLE PRECISION
! Processes multiple 2D slices (e.g., pressure levels × time steps)
!==============================================================================
subroutine z2geouv_4d(z, nlat, mlon, n3rd, n4th, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! INPUT - single precision
    integer, intent(in) :: nlat                    ! number of latitudes
    integer, intent(in) :: mlon                    ! number of longitudes
    integer, intent(in) :: n3rd                    ! size of 3rd dimension (e.g., pressure levels)
    integer, intent(in) :: n4th                    ! size of 4th dimension (e.g., time steps)
    integer, intent(in) :: iopt                    ! longitude boundary option (0=non-cyclic, 1=cyclic)
    real(kind=4), intent(in) :: z(nlat, mlon, n3rd, n4th)    ! geopotential height [gpm]
    real(kind=4), intent(in) :: zmsg                         ! missing value code
    real(kind=4), intent(in) :: glat(nlat)                   ! latitudes [degrees] (south to north)
    real(kind=4), intent(in) :: glon(mlon)                   ! longitudes [degrees]

    ! OUTPUT - single precision
    real(kind=4), intent(out) :: ug(nlat, mlon, n3rd, n4th)  ! zonal geostrophic wind [m/s]
    real(kind=4), intent(out) :: vg(nlat, mlon, n3rd, n4th)  ! meridional geostrophic wind [m/s]

    ! Local variables
    integer :: k, t

    ! Process each 2D slice along the 3rd and 4th dimensions
    ! OpenMP parallelization with loop collapsing for better load balancing
    !$omp parallel do private(k, t) collapse(2) if(n3rd * n4th > 16)
    do t = 1, n4th
        do k = 1, n3rd
            call z2geouv(z(:, :, k, t), nlat, mlon, zmsg, glon, glat, &
                         ug(:, :, k, t), vg(:, :, k, t), iopt)
        end do
    end do
    !$omp end parallel do

end subroutine z2geouv_4d

!==============================================================================
! North-to-south data reordering routine - SINGLE PRECISION
! Handles input data that is ordered from north to south by reordering
! to south-to-north before calculation, then reordering results back
!==============================================================================
subroutine zuvnew(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! INPUT - single precision
    integer, intent(in) :: nlat           ! number of latitudes
    integer, intent(in) :: mlon           ! number of longitudes
    integer, intent(in) :: iopt           ! longitude boundary option (0=non-cyclic, 1=cyclic)
    real(kind=4), intent(in) :: z(nlat, mlon)    ! geopotential height [gpm] (north to south)
    real(kind=4), intent(in) :: zmsg             ! missing value code
    real(kind=4), intent(in) :: glat(nlat)       ! latitudes [degrees] (north to south)
    real(kind=4), intent(in) :: glon(mlon)       ! longitudes [degrees]

    ! OUTPUT - single precision
    real(kind=4), intent(out) :: ug(nlat, mlon)  ! zonal geostrophic wind [m/s]
    real(kind=4), intent(out) :: vg(nlat, mlon)  ! meridional geostrophic wind [m/s]

    ! Local temporary arrays - single precision
    real(kind=4) :: ztmp(nlat, mlon), glattmp(nlat)
    real(kind=4) :: utmp(nlat), vtmp(nlat)
    integer :: nl, ml, nlat1

    nlat1 = nlat + 1

    ! Reorder from north-south to south-north [S ==> N] with SIMD
    do nl = 1, nlat
        glattmp(nlat1-nl) = glat(nl)
        !$omp simd
        do ml = 1, mlon
            ztmp(nlat1-nl, ml) = z(nl, ml)
        end do
        !$omp end simd
    end do

    ! Calculate geostrophic winds using reordered arrays (now south to north)
    call z2guv(ztmp, nlat, mlon, zmsg, glon, glattmp, ug, vg, iopt)

    ! Put results back in original order with SIMD
    ! Need to match original .f file logic: loop by longitude first
    do ml = 1, mlon
        !$omp simd
        do nl = 1, nlat
            utmp(nlat1-nl) = ug(nl, ml)
            vtmp(nlat1-nl) = vg(nl, ml)
        end do
        !$omp end simd
        !$omp simd
        do nl = 1, nlat
            ug(nl, ml) = utmp(nl)
            vg(nl, ml) = vtmp(nl)
        end do
        !$omp end simd
    end do

end subroutine zuvnew

!==============================================================================
! Core geostrophic wind calculation routine - SINGLE PRECISION
! Main calculation subroutine for geostrophic wind - SIMD OPTIMIZED
!
! Calculates geostrophic wind components using finite differences:
!   ug = -(g/f) * dZ/dy   where dZ/dy = (Z(j+1) - Z(j-1)) / (2*dy)
!   vg =  (g/f) * dZ/dx   where dZ/dx = (Z(i+1) - Z(i-1)) / (2*dx)
!
! Physical constants used:
!   g = 9.80616 m/s²     (gravity)
!   Ω = 7.292×10⁻⁵ s⁻¹   (Earth rotation rate)
!   Re = 6371220 m       (Earth radius)
!   f = 2Ω sin(lat)      (Coriolis parameter)
!==============================================================================
subroutine z2guv(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    ! INPUT - single precision
    integer, intent(in) :: nlat           ! number of latitudes
    integer, intent(in) :: mlon           ! number of longitudes
    integer, intent(in) :: iopt           ! longitude boundary option (0=non-cyclic, 1=cyclic)
    real(kind=4), intent(in) :: z(nlat, mlon)    ! geopotential height [gpm] (south to north)
    real(kind=4), intent(in) :: zmsg             ! missing value code
    real(kind=4), intent(in) :: glat(nlat)       ! latitudes [degrees] (south to north)
    real(kind=4), intent(in) :: glon(mlon)       ! longitudes [degrees]

    ! OUTPUT - single precision
    real(kind=4), intent(out) :: ug(nlat, mlon)  ! zonal geostrophic wind [m/s]
    real(kind=4), intent(out) :: vg(nlat, mlon)  ! meridional geostrophic wind [m/s]

    ! Local variables - all single precision
    real(kind=4) :: zz(0:nlat+1, 0:mlon+1)
    real(kind=4) :: gf(nlat), dx(nlat), dy(nlat)
    real(kind=4) :: inv_gf(nlat), inv_dx(nlat), inv_dy(nlat)
    real(kind=4), parameter :: gr = 9.80616e0
    real(kind=4), parameter :: omega = 7.292e-5
    real(kind=4), parameter :: re = 6371220.0e0
    real(kind=4), parameter :: pi = 4.0e0 * atan(1.0e0)
    real(kind=4), parameter :: rad = pi / 180.0e0
    real(kind=4), parameter :: two_omega = 2.0e0 * omega
    real(kind=4), parameter :: gr_over_two_omega = gr / two_omega
    real(kind=4) :: rr, dlon, dlon2
    real(kind=4) :: sin_lat, cos_lat, abs_lat
    integer :: nl, ml
    logical :: valid_lat(nlat)

    ! Calculate 'constant' quantities
    rr = re * rad
    dlon = rr * (glon(2) - glon(1))
    dlon2 = 2.0e0 * dlon

    ! Calculate the quantity [gravity/coriolis_force] and precompute inverses - SIMD optimized
    !$omp simd
    do nl = 1, nlat
        abs_lat = abs(glat(nl))
        valid_lat(nl) = (glat(nl) /= 0.0e0) .and. (abs_lat <= 89.999e0)

        if (valid_lat(nl)) then
            sin_lat = sin(glat(nl) * rad)
            cos_lat = cos(glat(nl) * rad)
            gf(nl) = gr_over_two_omega / sin_lat
            dx(nl) = dlon2 * cos_lat

            ! Precompute inverses to avoid division in inner loop
            inv_gf(nl) = 1.0e0 / gf(nl)
            if (dx(nl) /= 0.0e0) then
                inv_dx(nl) = 1.0e0 / dx(nl)
            else
                inv_dx(nl) = 0.0e0
            end if
        else
            gf(nl) = zmsg
            dx(nl) = 0.0e0
            inv_gf(nl) = 0.0e0
            inv_dx(nl) = 0.0e0
        end if
    end do
    !$omp end simd

    ! Get around 'rounding' problem at poles
    if (glat(1) < -89.999e0) then
        dx(1) = 0.0e0
        inv_dx(1) = 0.0e0
    end if
    if (glat(nlat) > 89.999e0) then
        dx(nlat) = 0.0e0
        inv_dx(nlat) = 0.0e0
    end if

    ! Latitude can have variable spacing - SIMD optimized
    dy(1) = (glat(2) - glat(1)) * rr
    inv_dy(1) = 1.0e0 / dy(1)

    !$omp simd
    do nl = 2, nlat - 1
        dy(nl) = (glat(nl+1) - glat(nl-1)) * rr
        inv_dy(nl) = 1.0e0 / dy(nl)
    end do
    !$omp end simd

    dy(nlat) = (glat(nlat) - glat(nlat-1)) * rr
    inv_dy(nlat) = 1.0e0 / dy(nlat)

    ! Set ug and vg to default values - vectorized
    ug = zmsg
    vg = zmsg

    ! Fill zz 'middle' - SIMD optimized
    do nl = 1, nlat
        !$omp simd
        do ml = 1, mlon
            zz(nl, ml) = z(nl, ml)
        end do
        !$omp end simd
    end do

    ! North and south boundaries (latitude boundaries)
    !$omp simd
    do ml = 1, mlon
        zz(0, ml) = zz(1, ml)            ! South boundary
        zz(nlat+1, ml) = zz(nlat, ml)    ! North boundary
    end do
    !$omp end simd

    ! East and west boundaries (longitude boundaries) - SIMD optimized
    if (iopt == 1) then
        ! Specify cyclic values
        !$omp simd
        do nl = 0, nlat+1
            zz(nl, 0) = zz(nl, mlon)     ! West boundary = East end
            zz(nl, mlon+1) = zz(nl, 1)   ! East boundary = West start
        end do
        !$omp end simd
    else
        ! Duplicate adjacent points
        !$omp simd
        do nl = 0, nlat+1
            zz(nl, 0) = zz(nl, 1)        ! West boundary
            zz(nl, mlon+1) = zz(nl, mlon)! East boundary
        end do
        !$omp end simd
    end if

    ! Main computation loop - SIMD optimized with precomputed values
    do nl = 1, nlat
        if (valid_lat(nl)) then
            !$omp simd
            do ml = 1, mlon
                ! Check for missing values and compute both components
                if (zz(nl+1, ml) /= zmsg .and. zz(nl-1, ml) /= zmsg .and. &
                    zz(nl, ml+1) /= zmsg .and. zz(nl, ml-1) /= zmsg) then

                    ! Use precomputed inverse values to avoid division
                    ! U component: -1/f * dz/dy (derivative in latitude direction)
                    ug(nl, ml) = -gf(nl) * (zz(nl+1, ml) - zz(nl-1, ml)) * inv_dy(nl)
                    ! V component: 1/f * dz/dx (derivative in longitude direction)
                    vg(nl, ml) = gf(nl) * (zz(nl, ml+1) - zz(nl, ml-1)) * inv_dx(nl)
                else
                    ! Handle partial missing data
                    if (zz(nl+1, ml) /= zmsg .and. zz(nl-1, ml) /= zmsg) then
                        ug(nl, ml) = -gf(nl) * (zz(nl+1, ml) - zz(nl-1, ml)) * inv_dy(nl)
                    end if
                    if (zz(nl, ml+1) /= zmsg .and. zz(nl, ml-1) /= zmsg) then
                        vg(nl, ml) = gf(nl) * (zz(nl, ml+1) - zz(nl, ml-1)) * inv_dx(nl)
                    end if
                end if
            end do
            !$omp end simd
        end if
    end do

    ! Adjust V component for non-cyclic longitude boundaries - SIMD optimized
    if (iopt /= 1) then
        !$omp simd
        do nl = 1, nlat
            if (vg(nl, 1) /= zmsg) vg(nl, 1) = vg(nl, 1) * 2.0e0
            if (vg(nl, mlon) /= zmsg) vg(nl, mlon) = vg(nl, mlon) * 2.0e0
        end do
        !$omp end simd
    end if

end subroutine z2guv
