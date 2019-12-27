!***********************************************************************
!>   Funftion that converts a binary integer to a character string
      Function Convert (intnum) 
!                                                                      *
!>   Written by  C. Froese Fischer                          June, 2019  
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!>    binary integer to be converted
      INTEGER,   INTENT(IN)    :: intnum
!>    character string for the decimal number
      Character(Len=12)        :: convert
!-----------------------------------------------

      write(convert,*) intnum
      convert = adjustl(convert)

      Return
      End Function Convert
