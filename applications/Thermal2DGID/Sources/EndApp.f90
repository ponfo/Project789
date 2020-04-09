module finishM
  use tools
  implicit none
  private
  public : finishProgram
  interface finishProgram
     procedure :: finishProgram
  end interface finishProgram
  integer(ikind), dimension(8) :: date_time
contains
  subroutine finishProgram()
    implicit none
    call date_and_time(VALUES=date_time)
    print'(A)','::::::::::::::: Finish Thermal2D :::::::::::::::'
    print'(A,I0,A,I0,A,I0)', 'Date: ', date_time(3), "/", date_time(2), "/", date_time(1)
    print'(A,I0,A,I0,A,I0)', 'Hour: ', date_time(5), ":", date_time(6), ":", date_time(7)
  end subroutine finishProgram
end module finishM
