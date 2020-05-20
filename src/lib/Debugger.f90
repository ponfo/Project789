module DebuggerM
  use UtilitiesM
  implicit none
  public
  
  interface debugLog
     procedure :: printStr
     procedure :: printStrInt
     procedure :: printStrR
     procedure :: printStrIntStr
     procedure :: printStrStr
  end interface debugLog
  
  integer(ikind)             , parameter    :: logUnit = 95
  integer(ikind)             , dimension(8) :: dateAndTime
  procedure(printStrON)      , pointer      :: printStr
  procedure(printStrIntON)   , pointer      :: printStrInt
  procedure(printStrRON)     , pointer      :: printStrR
  procedure(printStrIntStrON), pointer      :: printStrIntStr
  procedure(printStrStrON)   , pointer      :: printStrStr
  
contains
  
  subroutine initLog(isWorking, logFile)
    implicit none
    logical, intent(in) :: isWorking
    character(*), optional, intent(in) :: logFile
    if(isWorking) then
       printStr => printStrON
       printStrInt => printStrIntON
       printStrR => printStrRON
       printStrIntStr => printStrIntStrON
       printStrStr => printStrStrON
       open(logUnit, file = trim(logFile))
    else
       printStr => printStrOFF
       printStrInt => printStrIntOFF
       printStrR => printStrROFF
       printStrIntStr => printStrIntStrOFF
       printStrStr => printStrStrOFF
    end if
    call date_and_time(VALUES = dateAndTime)
    write(logUnit,'(I2,A,I2,A,I2,2X,A)') dateAndTime(5),":",dateAndTime(6),":",dateAndTime(7) &
         , 'Initiating Program'
  end subroutine initLog
  
  subroutine printStrON(str)
    implicit none
    character(*), intent(in) :: str
    call date_and_time(VALUES = dateAndTime)
    write(logUnit,'(I2,A,I2,A,I2,2X,A)') dateAndTime(5),":",dateAndTime(6),":",dateAndTime(7), str
  end subroutine printStrON
  
  subroutine printStrOFF(str)
    implicit none
    character(*), intent(in) :: str
  end subroutine printStrOFF
  
  subroutine printStrIntON(str, int)
    implicit none
    character(*)  , intent(in) :: str
    integer(ikind), intent(in) :: int
    call date_and_time(VALUES = dateAndTime)
    write(logUnit,'(I2,A,I2,A,I2,2X,A,I0)') dateAndTime(5),":",dateAndTime(6),":",dateAndTime(7), str, int
  end subroutine printStrIntON
  
  subroutine printStrIntOFF(str, int)
    implicit none
    character(*)  , intent(in) :: str
    integer(ikind), intent(in) :: int
  end subroutine printStrIntOFF

  subroutine printStrRON(str, r)
    implicit none
    character(*), intent(in) :: str
    real(rkind) , intent(in) :: r
    call date_and_time(VALUES = dateAndTime)
    write(logUnit,'(I2,A,I2,A,I2,2X,A,E16.8)') dateAndTime(5),":",dateAndTime(6),":",dateAndTime(7), str, r
  end subroutine printStrRON

  subroutine printStrROFF(str, r)
    implicit none
    character(*), intent(in) :: str
    real(rkind) , intent(in) :: r
  end subroutine printStrROFF
  
  subroutine printStrIntStrON(str1, int, str2)
    implicit none
    character(*)  , intent(in) :: str1
    integer(ikind), intent(in) :: int
    character(*)  , intent(in) :: str2
    call date_and_time(VALUES = dateAndTime)
    write(logUnit,'(I2,A,I2,A,I2,2X,A,I0,A)') dateAndTime(5),":",dateAndTime(6),":",dateAndTime(7), str1, int, str2
  end subroutine printStrIntStrON
  
  subroutine printStrIntStrOFF(str1, int, str2)
    implicit none
    character(*)  , intent(in) :: str1
    integer(ikind), intent(in) :: int
    character(*)  , intent(in) :: str2
  end subroutine printStrIntStrOFF
  
  subroutine printStrStrON(str1, str2)
    implicit none
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    call date_and_time(VALUES = dateAndTime)
    write(logUnit,'(I2,A,I2,A,I2,2X,A,A)') dateAndTime(5),":",dateAndTime(6),":",dateAndTime(7), str1, str2
  end subroutine printStrStrON
  
  subroutine printStrStrOFF(str1, str2)
    implicit none
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
  end subroutine printStrStrOFF
  
end module DebuggerM

