! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-17
! File    : rcm_geodesy.f90
! Purpose : Provide great-circle distance helpers shared by the interp
!           curvilinear remapping kernels.
! Notes   : Both the internal core routine and the public wrapper keep the
!           same spherical-Earth arithmetic as the historical fixed-form code.
! =============================================================================
!
module rcm_geodesy_core
    implicit none

    integer, parameter :: real64 = selected_real_kind(15, 307)
contains

    ! Compute the great-circle separation between two latitude/longitude
    ! points using the same spherical-Earth arithmetic as the historical
    ! fixed-form kernels. The unit selector matches the original interface:
    ! 1=radians, 2=degrees, 3=meters, 4=kilometers.
    real(real64) function dgcdist_core(rlat1, rlon1, rlat2, rlon2, iu) result(distance)
        real(real64), intent(in) :: rlat1, rlon1, rlat2, rlon2
        integer, intent(in) :: iu
        real(real64), parameter :: units(5) = [ &
            1.0_real64, 57.29577995691645_real64, 6371220.0_real64, &
            6371.2200_real64, 0.0_real64 &
        ]
        real(real64), parameter :: rad = 0.01745329238474369_real64
        real(real64) :: dlonr, rlat1r, rlat2r

        if (rlat1 == rlat2 .and. rlon1 == rlon2) then
            distance = 0.0_real64
            return
        end if

        rlat1r = rlat1 * rad
        rlat2r = rlat2 * rad
        dlonr = min(abs(rlon1 - rlon2), abs(360.0_real64 - rlon1 + rlon2), &
                    abs(360.0_real64 - rlon2 + rlon1)) * rad

        distance = atan2(sqrt((cos(rlat2r) * sin(dlonr)) ** 2 + &
                              (cos(rlat1r) * sin(rlat2r) - &
                               sin(rlat1r) * cos(rlat2r) * cos(dlonr)) ** 2), &
                         sin(rlat1r) * sin(rlat2r) + &
                         cos(rlat1r) * cos(rlat2r) * cos(dlonr)) * units(iu)
    end function dgcdist_core
end module rcm_geodesy_core


! QUICK REFERENCE
! PURPOSE
!    PUBLIC F2PY-FRIENDLY ENTRY POINT FOR GREAT-CIRCLE DISTANCE ON THE
!    SAME SPHERICAL EARTH USED BY THE LEGACY RCM INTERPOLATION KERNELS.
!
! INPUTS
!    RLAT1, RLON1 - FIRST POINT IN DEGREES
!    RLAT2, RLON2 - SECOND POINT IN DEGREES
!    IU           - OUTPUT UNIT SELECTOR
!                   1: RADIANS
!                   2: DEGREES
!                   3: METERS
!                   4: KILOMETERS
!
! OUTPUT
!    DGCDIST RETURNS THE ARC DISTANCE BETWEEN THE TWO INPUT POINTS.
real(real64) function dgcdist(rlat1, rlon1, rlat2, rlon2, iu)
    use rcm_geodesy_core, only : real64, dgcdist_core
    implicit none

    real(real64), intent(in) :: rlat1, rlon1, rlat2, rlon2
    integer, intent(in) :: iu

    dgcdist = dgcdist_core(rlat1, rlon1, rlat2, rlon2, iu)
end function dgcdist
