SUBROUTINE  getrlist(nlist, list, str, nrlist, rlist)


!>   Determine the list of orbitals the user selected.
!>   List   :: orbitls to be determined
!>   nlist  :: number of orbitals
!>   str    :: user's reponse to orbitals for an option
!!     The basic input format is a list of orbitals, separated by 
!!     space or comma, with a  possible wild card.  * is equivalent to 
!!     all orbitals, n* to all orbitals with a given n, or orbitals 
!!     in non-relativistic notation
!>   nrlist :: number of orbitals selected
!>   rlist  :: list of orbitals selected

  USE  load_mod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlist
  INTEGER, DIMENSION(nw), INTENT(in) :: list
  CHARACTER(len=4*nnnw), INTENT(inout)      :: str
  INTEGER, INTENT(out) ::  nrlist
  INTEGER, DIMENSION(nw), INTENT(out) :: rlist
  
  !      local variables
  INTEGER :: i, ipos, iorb,  len, lens, n
  
  
  nrlist = 0
  rlist = 0
  DO
     str = ADJUSTL(str)   
     len = len_TRIM(str)
     PRINT *, 'len =',  len
     IF (len .EQ. 0) EXIT
     ! scan for blank or comma
     ipos = SCAN(str(1:len), ' ,')   
     PRINT *, 'S1: ipos =', str(1:len), ipos
     IF (ipos .GT. 0) THEN
        lens = ipos-1
     ELSE
        str(len+1:len+1) = ' '
        lens = len + 1
     END IF
     ! when all orbitals are slected
     IF (str(1:lens) .EQ. '*') THEN
        nrlist = nlist
        rlist = list
        ! all orbitals selected
        EXIT   
     ELSE
        ipos = INDEX(str(1:lens),'*')       
        !  we have n* format
        IF (ipos .GT. 0) THEN     
           PRINT *, ipos
           READ(str(1:ipos-1), *) n 
           ! find all orbits with np(i) = n 
           DO i = 1,nw
              IF (np(i) .EQ. n) THEN
                 nrlist = nrlist +1
                 rlist(nrlist) = i
              END IF
           END DO
        ELSE
           ! we have a separate orbital  (say 2p which includes ( 2p-, 2p )
           iorb = INDEX(orblist, str(1:lens)) 
           IF (iorb .GT. 0) THEN
              nrlist = nrlist+1
              rlist(nrlist) = iorb/4 +1
              IF (lens .LT. 3) THEN
                 iorb = INDEX(orblist(ipos+1:nw), str(1:lens))
                 IF (iorb .GT. 0) THEN
                    nrlist = nrlist+1
                    rlist(nrlist) = iorb/4 + 1
                 END IF
              END IF
           ELSE 
              PRINT *, 'Orbital ', str(1:lens), ' not found in list'
              STOP
           END IF
        END IF
     END IF
     ! blank out the first lens characters of str
     
     DO i = 1, MIN(lens+1,len) 
        str(i:i) = ' '
     END DO
  END DO
  
END SUBROUTINE getrlist
