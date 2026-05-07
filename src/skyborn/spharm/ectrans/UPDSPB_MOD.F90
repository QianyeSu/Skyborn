! (C) Copyright 1988- ECMWF.
! (C) Copyright 2013- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSPB_MOD
CONTAINS
SUBROUTINE UPDSPB(KM,KFIELD,POA,PSPEC,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KM,KFIELD
REAL(KIND=JPRB)   ,INTENT(IN)  :: POA(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPEC(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN, ISMAX, ITMAX, IASM0,IFLD

ISMAX = R%NSMAX
ITMAX = R%NTMAX
IASM0 = D%NASM0(KM)

IF(KM == 0) THEN
  IF(PRESENT(KFLDPTR)) THEN
    DO JFLD=1,KFIELD
      IR = 2*JFLD-1
      IFLD = KFLDPTR(JFLD)
      DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
        INM = IASM0+(ITMAX+2-JN)*2
        PSPEC(IFLD,INM)   = POA(JN,IR)
        PSPEC(IFLD,INM+1) = 0.0_JPRB
      ENDDO
    ENDDO
  ELSE
    DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
      INM = IASM0+(ITMAX+2-JN)*2
!DIR$ IVDEP
!OCL NOVREC
      DO JFLD=1,KFIELD
        IR = 2*JFLD-1
        PSPEC(JFLD,INM)   = POA(JN,IR)
        PSPEC(JFLD,INM+1) = 0.0_JPRB
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF(PRESENT(KFLDPTR)) THEN
    DO JFLD=1,KFIELD
      IR = 2*JFLD-1
      II = IR+1
      IFLD = KFLDPTR(JFLD)
      DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
        INM = IASM0+((ITMAX+2-JN)-KM)*2
        PSPEC(IFLD,INM)   = POA(JN,IR)
        PSPEC(IFLD,INM+1) = POA(JN,II)
      ENDDO
    ENDDO
  ELSE
    DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
      INM = IASM0+((ITMAX+2-JN)-KM)*2
!DIR$ IVDEP
!OCL NOVREC
      DO JFLD=1,KFIELD
        IR = 2*JFLD-1
        II = IR+1
        PSPEC(JFLD,INM)   = POA(JN,IR)
        PSPEC(JFLD,INM+1) = POA(JN,II)
      ENDDO
    ENDDO
  ENDIF
ENDIF

END SUBROUTINE UPDSPB
END MODULE UPDSPB_MOD
