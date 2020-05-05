module GIDDataOutputM
  use UtilitiesM
  implicit none
  private
  public :: printResults, printForCoupling, finishProgram
  interface printResults
     procedure :: printResultsVec1
     procedure :: printResults1DVec2
     procedure :: printResults2DVec2
     procedure :: printResults2DVec1
  end interface printResults
  interface printForCoupling
     procedure :: printForCoupling
  end interface printForCoupling
  interface finishProgram
     procedure :: finishProgram
  end interface finishProgram
  integer(ikind), parameter    :: projectData = 80
  integer(ikind), parameter    :: results = 3
  integer(ikind), dimension(8) :: date_time
  character(100)               :: projectName, path
contains
  subroutine init()
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(A)') projectName
    read(projectData, '(A)') path
    close(projectData)
  end subroutine init
  subroutine printResultsVec1(resultName, step, graphType, locationName, resultNumber&
       , component1)
    implicit none
    integer(ikind)                 :: iPoint
    integer(ikind), intent(in)     :: step, resultNumber
    real(rkind), intent(in), dimension(resultNumber*2) :: component1
    character(*), intent(in) :: resultName, graphType, locationName
    call init()
    open(results, file = trim(projectName)//'.flavia.res')
    write(results,'(A)')      'GiD Post Result File 1.0'
    write(results,*)   'Result "',trim(resultName),'" "',trim(projectName)&
         ,'" ',step,' ',trim(graphType),' ',trim(locationName)
    write(results,*)   'Values'
    if(trim(graphType) == 'Scalar') then
       Do iPoint = 1, resultNumber
          Write(results,*) iPoint, component1(iPoint)
       End Do
    else if (trim(graphType) == 'Vector') then
       do iPoint = 1, resultNumber
          write(results,*) iPoint, component1(2*iPoint-1), component1(2*iPoint)
       end do
    end if
    write(results,'(A)')   'End Values'
  end subroutine printResultsVec1
  subroutine printResults1DVec2(resultName, type, step, graphType, locationName, gaussPoints &
       , resultNumber, elemID, component1, component2)
    implicit none
    character(*), intent(in) :: resultName
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: step
    character(*), intent(in) :: graphType
    character(*), intent(in) :: locationName
    real(rkind), dimension(:), intent(in) :: gaussPoints
    integer(ikind), intent(in) :: resultNumber
    integer(ikind), dimension(resultNumber), intent(in) :: elemID
    real(rkind), dimension(:), intent(in) :: component1
    real(rkind), dimension(:), intent(in) :: component2
    real(rkind) :: prom(2)
    integer(ikind) :: i, j, k, count, numberGP
    if(resultNumber == 0) return
    call init()
    write(results,'(/,3A)') 'GaussPoints "Points'//trim(resultName), '" ElemType ', trim(type)
    write(results,'(A,I0)') 'Number of GaussPoints: ', size(gaussPoints)
    write(results,'(A)') 'Nodes not included'
    write(results,'(A)') 'Natural Coordinates: Internal'
!!$    do i = 1, size(gaussPoints,1)
!!$       write(results,'(F26.16,2X,F26.16)') gaussPoints(i), 0.0
!!$    end do
    write(results,'(A)') 'End gausspoints'
    write(results,'(5A,I0,6A)') 'Result "', trim(resultName), '" "', trim(projectName), '" ', step &
         , ' ', trim(graphType), ' ', trim(locationName), ' "Points'//trim(resultName) , '"'
    write(results,'(A)') 'Values'
    count = 0
    do i = 1, resultNumber
       count = count + 1
       write(results,'(I0,2X,F26.16,2X,F26.16)') elemID(i), component1(count), 0.0
       do j = 2, size(gaussPoints)
          count = count + 1
          write(results,'(6X,F26.16,2X,F26.16)') component1(count), 0.0
       end do
    end do
    write(results,'(A)') 'End Values'
  end subroutine printResults1DVec2
  subroutine printResults2DVec2(resultName, type, step, graphType, locationName, gaussPoints &
       , resultNumber, elemID, component1, component2)
    implicit none
    character(*), intent(in) :: resultName
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: step
    character(*), intent(in) :: graphType
    character(*), intent(in) :: locationName
    real(rkind), dimension(:,:), intent(in) :: gaussPoints
    integer(ikind), intent(in) :: resultNumber
    integer(ikind), dimension(resultNumber), intent(in) :: elemID
    real(rkind), dimension(:), intent(in) :: component1
    real(rkind), dimension(:), intent(in) :: component2
    real(rkind) :: prom(2)
    integer(ikind) :: i, j, k, count, numberGP
    if(resultNumber == 0) return
    call init()
    write(results,'(/,3A)') 'GaussPoints "Points'//trim(resultName), '" ElemType ', trim(type)
    write(results,'(A,I0)') 'Number of GaussPoints: ', size(gaussPoints,1)
    write(results,'(A)') 'Natural Coordinates: Given'
    do i = 1, size(gaussPoints,1)
       write(results,'(F26.16,2X,F26.16)') gaussPoints(i,1), gaussPoints(i,2)
    end do
    write(results,'(A)') 'End gausspoints'
    write(results,'(5A,I0,6A)') 'Result "', trim(resultName), '" "', trim(projectName), '" ', step &
         , ' ', trim(graphType), ' ', trim(locationName), ' "Points'//trim(resultName) , '"'
    write(results,'(A)') 'Values'
    count = 0
    do i = 1, resultNumber
       count = count + 1
       write(results,'(I0,2X,F26.16,2X,F26.16)') elemID(i), component1(count), component2(count)
       do j = 2, size(gaussPoints,1)
          count = count + 1
          write(results,'(6X,F26.16,2X,F26.16)') component1(count), component2(count)
       end do
    end do
    write(results,'(A)') 'End Values'
  end subroutine printResults2DVec2
  subroutine printResults2DVec1(resultName, type, step, graphType, locationName, gaussPoints &
       , resultNumber, elemID, component1)
    implicit none
    character(*), intent(in) :: resultName
    character(*), intent(in) :: type
    integer(ikind), intent(in) :: step
    character(*), intent(in) :: graphType
    character(*), intent(in) :: locationName
    real(rkind), dimension(:,:), intent(in) :: gaussPoints
    integer(ikind), intent(in) :: resultNumber
    integer(ikind), dimension(resultNumber), intent(in) :: elemID
    real(rkind), dimension(:), intent(in) :: component1
    real(rkind) :: prom(2)
    integer(ikind) :: i, j, k, count, numberGP
    if(resultNumber == 0) return
    call init()
    write(results,'(/,3A)') 'GaussPoints "Points'//trim(resultName), '" ElemType ', trim(type)
    write(results,'(A,I0)') 'Number of GaussPoints: ', size(gaussPoints,1)
    write(results,'(A)') 'Natural Coordinates: Given'
    do i = 1, size(gaussPoints,1)
       write(results,'(F26.16,2X,F26.16)') gaussPoints(i,1), gaussPoints(i,2)
    end do
    write(results,'(A)') 'End gausspoints'
    write(results,'(5A,I0,6A)') 'Result "', trim(resultName), '" "', trim(projectName), '" ', step &
         , ' ', trim(graphType), ' ', trim(locationName), ' "Points'//trim(resultName) , '"'
    write(results,'(A)') 'Values'
    count = 0
    do i = 1, resultNumber
       count = count + 1
       write(results,'(I0,2X,E26.16)') elemID(i), component1(count)
       do j = 2, size(gaussPoints,1)
          count = count + 1
          write(results,'(6X,E26.16)') component1(count)
       end do
    end do
    write(results,'(A)') 'End Values'
  end subroutine printResults2DVec1

  subroutine printForCoupling(nPoint, temperature)
    implicit none
    integer(ikind), intent(in) :: nPoint
    real(rkind), dimension(:), intent(in) :: temperature
    integer(ikind) :: i
    open(8, file = 'Temperatures.dat')
    write(8,'(A)') 'Temperatures for thermal-structural coupling'
    write(8,'(A)') 'Number of nodes'
    write(8,'(I0)') nPoint
    write(8,'(A)')'Node    Temperature'
    do i = 1, nPoint
       write(8,'(I7,E16.8)') i, temperature(i)
    end do
  end subroutine printForCoupling

  subroutine finishProgram()
    implicit none
    close(results)
    !call free()
    call date_and_time(VALUES=date_time)
    print'(A)','::::::::::::::: Finish FEM2D :::::::::::::::'
    print'(A,I0,A,I0,A,I0)', 'Date: ', date_time(3), "/", date_time(2), "/", date_time(1)
    print'(A,I0,A,I0,A,I0)', 'Hour: ', date_time(5), ":", date_time(6), ":", date_time(7)
  end subroutine finishProgram
end module GIDDataOutputM
