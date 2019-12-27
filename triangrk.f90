  LOGICAL FUNCTION TRIANGRK (LA, K, LB) 
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!   M o d u l e s 
!-----------------------------------------------
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER(kind=1), INTENT(IN) :: LA, LB
    INTEGER, INTENT(IN) :: K 
!
!   Perform the triangularity check
!
    IF (MOD(K + LA + LB,2) /= 0) THEN 
       TRIANGRK = .FALSE. 
    ELSE 
       IF (ABS(LA - LB) > K) THEN 
          TRIANGRK = .FALSE. 
       ELSE IF (LA + LB < K) THEN 
          TRIANGRK = .FALSE. 
       ELSE 
          TRIANGRK = .TRUE. 
       ENDIF
    ENDIF
    RETURN  
  END FUNCTION TRIANGRK
