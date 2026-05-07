! (C) Copyright 1991- ECMWF.
! (C) Copyright 2013- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LDFOU2_MOD
CONTAINS
SUBROUTINE LDFOU2(KM,KF_UV,PAIA,PSIA)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FIELDS      ,ONLY : F
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM,KF_UV

REAL(KIND=JPRB) ,INTENT(INOUT) :: PSIA(:,:),   PAIA(:,:)

INTEGER(KIND=JPIM) :: J, JGL ,IFLD ,ISL

ISL  = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IFLD = 4*KF_UV

DO JGL=ISL,R%NDGNH
  DO J=1,IFLD
    PAIA(J,JGL) = PAIA(J,JGL)*F%RACTHE(JGL)
    PSIA(J,JGL) = PSIA(J,JGL)*F%RACTHE(JGL)
  ENDDO
ENDDO

END SUBROUTINE LDFOU2
END MODULE LDFOU2_MOD
