!***********************************************************************
!                                                                      *
        SUBROUTINE INTRPQ(PA, QA, MA, RA, J, DNORM) 
!                                                                      *
!   This  subprogram  interpolates  the  arrays  PA(1:MA), QA(1:MA),   *
!   tabulated on grid RA(1:MA) into the COMMON arrays PF(1:MF(J),J),   *
!   QF(1:MF(J),J). (Aitken's  algorithm is used. See F B Hildebrand,   *
!   Introduction  to  Numerical  Analysis, 2nd ed., McGraw-Hill, New   *
!   York, NY, 1974.) The orbital is renormalized.                      *
!                                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:37:49   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17

          IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
          INTEGER, INTENT(IN) :: MA 
          INTEGER  :: J 
          REAL(kind=8), INTENT(OUT) :: DNORM 
          REAL(kind=8), DIMENSION(*), INTENT(IN) :: PA, QA, RA
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
          INTEGER, PARAMETER :: MXORD = 13 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
          INTEGER :: I, MFJ, NRSTLO, KOUNT, IROW, K, NRSTHI, LLO, LHI, LOCNXT, &
               ILIROK, ILDIAG, ILOTHR, MFJP1, n, nnnp
          REAL(kind=8) :: RAMA, accy
          INTEGER, DIMENSION(NNNW) :: MF       
          REAL(kind=8), DIMENSION(MXORD) :: X, DX 
          REAL(kind=8), DIMENSION((MXORD*(MXORD + 1))/2) :: POLYP, POLYQ 
          REAL(kind=8) :: RN, XBAR, PESTL, QESTL, DIFF, DIFFT, DXKMN1, DXIROW, &
               FACTOR, PESTT, QESTT, DPBP, DQBQ, DNFAC
          REAL(kind = 8), DIMENSION(npx), ALLOCATABLE :: ta
          LOGICAL :: SET 
          LOGICAL, DIMENSION(npx) :: USED 
!-----------------------------------------------
!
!
!   This function dispenses with the need for a two-dimensional
!   array
!
!
!   Initialization
!
          nnnp = npx
          n    = nnnp
          ACCY = h**6 
          RAMA = RA(MA) 
          RN = R(N) 
!
!   This is always true in GRASP
!
          P(1,J) = 0.0D00 
          Q(1,J) = 0.0D00 
!
!   Checks
!
          IF (RAMA > RN) THEN 
             WRITE (*, 300) RN, RAMA 
             PRINT*, "N =",N,"MA =",MA
             STOP  
          ENDIF
          !
!   Determine end of grid
!
          I = N 
1         CONTINUE 
          I = I - 1 
          IF (R(I) <= RAMA) THEN 
             MFJ = I 
          ELSE 
             GO TO 1 
          ENDIF
          MF(J) = MFJ 
!
!   Overall initialization for interpolation
!
          NRSTLO = 0 
          KOUNT = 0 
!
!   Perform interpolation
!
          DO I = 2, MFJ 
!
!   Initialization for interpolation
!
             XBAR = R(I) 
             IROW = 0 
             PESTL = 0.0D00 
             QESTL = 0.0D00 
!
!   Determine the nearest two grid points bounding the present
!   grid point
!
2            CONTINUE 
             K = NRSTLO + 1 
             IF (RA(K) < XBAR) THEN 
                NRSTLO = K 
                GO TO 2 
             ELSE 
                NRSTHI = K 
             ENDIF
!
!   Clear relevant piece of use-indicator array
!
             LLO = MAX(NRSTLO - MXORD,1) 
             LHI = MIN(NRSTHI + MXORD,MA) 
             USED(LLO:LHI) = .FALSE. 
!
!   Determine next-nearest grid point
!
4            CONTINUE 
             IROW = IROW + 1 
             LLO = MAX(NRSTLO - IROW + 1,1) 
             LHI = MIN(NRSTHI + IROW - 1,MA) 
             SET = .FALSE. 
             DO K = LLO, LHI 
                IF (USED(K)) CYCLE  
                IF (.NOT.SET) THEN 
                   DIFF = RA(K) - XBAR 
                   LOCNXT = K 
                   SET = .TRUE. 
                ELSE 
                   DIFFT = RA(K) - XBAR 
                   IF (ABS(DIFFT) < ABS(DIFF)) THEN 
                      DIFF = DIFFT 
                      LOCNXT = K 
                   ENDIF
                ENDIF
             END DO
             USED(LOCNXT) = .TRUE. 
             X(IROW) = RA(LOCNXT) 
             DX(IROW) = DIFF 
!
!   Fill table for this row
!
             DO K = 1, IROW 
                ILIROK = (irow*(irow - 1))/2 + k 
 
                IF (K == 1) THEN 
                   POLYP(ILIROK) = PA(LOCNXT) 
                   POLYQ(ILIROK) = QA(LOCNXT) 
                ELSE 
                   ILDIAG = ((k-1)*((k-1) - 1))/2 + (k-1) 
                   ILOTHR = (irow*(irow - 1))/2 + (k - 1) 
                   DXKMN1 = DX(K-1) 
                   DXIROW = DX(IROW) 
                   FACTOR = 1.0D00/(X(IROW)-X(K-1)) 
                   POLYP(ILIROK) = (POLYP(ILDIAG)*DXIROW-POLYP(ILOTHR)*DXKMN1)*&
                        FACTOR 
                   POLYQ(ILIROK) = (POLYQ(ILDIAG)*DXIROW-POLYQ(ILOTHR)*DXKMN1)*&
                        FACTOR 
                ENDIF
             END DO
!
!   Check for convergence
!
             ILDIAG = (irow*(irow - 1))/2 + irow 
             PESTT = POLYP(ILDIAG) 
             QESTT = POLYQ(ILDIAG) 
             IF (PESTT==0.0D00 .OR. QESTT==0.0D00) THEN 
                IF (IROW < MXORD) THEN 
                   GO TO 4 
                ELSE 
                   P(I,J) = PESTT 
                   Q(I,J) = QESTT 
                ENDIF
             ELSE 
                DPBP = ABS((PESTT - PESTL)/PESTT) 
                DQBQ = ABS((QESTT - QESTL)/QESTT) 
                IF (DQBQ<ACCY .AND. DPBP<ACCY) THEN 
                   P(I,J) = PESTT 
                   Q(I,J) = QESTT 
                ELSE 
                   PESTL = PESTT 
                   QESTL = QESTT 
                   IF (IROW < MXORD) THEN 
                      GO TO 4 
                   ELSE 
                      P(I,J) = PESTT 
                      Q(I,J) = QESTT 
                      KOUNT = KOUNT + 1 
                   ENDIF
                ENDIF
             ENDIF
!
          END DO
!
!   Ensure that all points of the array are defined by setting the
!   tail to zero
!
          MFJP1 = MFJ + 1 
          P(MFJP1:N,J) = 0.0D00 
          Q(MFJP1:N,J) = 0.0D00 
!
!          IF (LDBPR(3) .AND. KOUNT>0) WRITE (99, 301) ACCY, KOUNT, MFJ 
!
!   Normalization
!

          ta(:) = (p(:,j)*p(:,j) + q(:,j)*q(:,j))*s(:)           
          DNORM = quad(ta) 
          DNFAC = 1.0D00/SQRT(DNORM) 
          P(:MFJ,J) = P(:MFJ,J)*DNFAC 
          Q(:MFJ,J) = Q(:MFJ,J)*DNFAC 
!
          RETURN  
!
300       FORMAT(/,'INTRPQ: Grid of insufficient extent:'/,&
               ' Present grid has R(N) = ',1P,1D19.12,' Bohr radii'/,&
               '          Require R(N) = ',1D19.12,' Bohr radii') 
301       FORMAT(/,'INTRPQ: Interpolation procedure not converged to',1P,1D19.12,&
               ' for ',1I3,' of ',1I3,' tabulation points') 
          RETURN  

        END SUBROUTINE INTRPQ
