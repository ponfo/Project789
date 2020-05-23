module GIDDataOutputM
  use UtilitiesM
  implicit none
  private
  public :: printResults, finishProgram, initDataOutput
  interface printResults
     procedure :: printResultsVec1
     procedure :: printResultsVec2
     procedure :: printResults1DVec2
     procedure :: printResults2DVec1
     procedure :: printResults2DVec2
     procedure :: printResults3DVec3
     procedure :: printResultsVecOnNodesNDof
  end interface printResults
  interface finishProgram
     procedure :: finishProgram
  end interface finishProgram
  integer(ikind), parameter    :: projectData = 80
  integer(ikind), parameter    :: results = 3
  integer(ikind), dimension(8) :: date_time
  character(100)               :: projectName
contains
  subroutine initDataOutput()
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(A)') projectName
    close(projectData)
    open(results, file = trim(projectName)//'.flavia.res')
    write(results,'(A)')      'GiD Post Result File 1.0'
  end subroutine initDataOutput
  
  subroutine printResultsVec1(resultName, step, graphType, locationName, resultNumber&
       , component1)
    implicit none
    integer(ikind)                 :: iPoint
    integer(ikind), intent(in)     :: step, resultNumber
    real(rkind), intent(in), dimension(resultNumber) :: component1
    character(*), intent(in) :: resultName, graphType, locationName
    write(results,*)   'Result "',trim(resultName),'" "',trim(projectName)&
         ,'" ',step,' ',trim(graphType),' ',trim(locationName)
    write(results,*)   'Values'
    Do iPoint = 1, resultNumber
       Write(results,*) iPoint, component1(iPoint)
    End Do
    write(results,*)   'End Values'
  end subroutine printResultsVec1

  subroutine printResultsVec2(resultName, step, graphType, locationName, resultNumber&
       , component1, component2)
    implicit none
    integer(ikind)                 :: iPoint
    integer(ikind), intent(in)     :: step, resultNumber
    real(rkind), intent(in), dimension(resultNumber) :: component1
    real(rkind), intent(in), dimension(resultNumber) :: component2
    character(*), intent(in) :: resultName, graphType, locationName
    write(results,*)   'Result "',trim(resultName),'" "',trim(projectName)&
         ,'" ',step,' ',trim(graphType),' ',trim(locationName)
    write(results,*)   'Values'
    Do iPoint = 1, resultNumber
       Write(results,*) iPoint, component1(iPoint), component2(iPoint)
    End Do
    write(results,*)   'End Values'
  end subroutine printResultsVec2
  
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
    write(results,'(/,3A)') 'GaussPoints "Points'//trim(resultName), '" ElemType ', trim(type)
    write(results,'(A,I0)') 'Number of GaussPoints: ', size(gaussPoints,1)
    write(results,'(A)') 'Natural Coordinates: Given'
    do i = 1, size(gaussPoints,1)
       write(results,'(E26.16,2X,E26.16)') gaussPoints(i,1), gaussPoints(i,2)
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
  
  subroutine printResults3DVec3(resultName, type, step, graphType, locationName, gaussPoints &
       , resultNumber, elemID, component1, component2, component3)
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
    real(rkind), dimension(:), intent(in) :: component3
    real(rkind) :: prom(2)
    integer(ikind) :: i, j, k, count, numberGP
    if(resultNumber == 0) return
    write(results,'(/,3A)') 'GaussPoints "Points'//trim(resultName), '" ElemType ', trim(type)
    write(results,'(A,I0)') 'Number of GaussPoints: ', size(gaussPoints,1)
    write(results,'(A)') 'Natural Coordinates: Given'
    do i = 1, size(gaussPoints,1)
       write(results,'(E26.16,2(2X,E26.16))') gaussPoints(i,1), gaussPoints(i,2), gaussPoints(i,3)
    end do
    write(results,'(A)') 'End gausspoints'
    write(results,'(5A,I0,6A)') 'Result "', trim(resultName), '" "', trim(projectName), '" ', step &
         , ' ', trim(graphType), ' ', trim(locationName), ' "Points'//trim(resultName) , '"'
    write(results,'(A)') 'Values'
    count = 0
    do i = 1, resultNumber
       count = count + 1
       write(results,'(I0,3(2X,E26.16))') elemID(i), component1(count), component2(count), component3(count)
       do j = 2, size(gaussPoints,1)
          count = count + 1
          write(results,'(6X,E26.16,2X,E26.16,2X,E26.16)') component1(count), component2(count), component3(count)
       end do
    end do
    write(results,'(A)') 'End Values'
  end subroutine printResults3DVec3

  subroutine finishProgram()
    implicit none
    close(results)
    call date_and_time(VALUES=date_time)
    print'(A)','::::::::::::::: Finish FEM2D :::::::::::::::'
    print'(A,I0,A,I0,A,I0)', 'Date: ', date_time(3), "/", date_time(2), "/", date_time(1)
    print'(A,I0,A,I0,A,I0)', 'Hour: ', date_time(5), ":", date_time(6), ":", date_time(7)
  end subroutine finishProgram

  subroutine printResultsVecOnNodesNDof(resultName, step, graphType, locationName, resultNumber&
       , nDof, component1)
    implicit none
    integer(ikind)                 :: iPoint
    integer(ikind), intent(in)     :: step, resultNumber, nDof
    real(rkind), intent(in), dimension(:) :: component1
    character(*), intent(in) :: resultName, graphType, locationName
    integer(ikind) :: i, count
    write(results,*)   'Result "',trim(resultName),'" "',trim(projectName)&
         ,'" ',step,' ',trim(graphType),' ',trim(locationName)
    write(results,*)   'Values'
    count = 0
    Do iPoint = 1, resultNumber
       count = count + 1
       write(results,'(I0,*(E16.8))') iPoint, (component1(count*nDof-(nDof-i)),i=1,nDof)
    End Do
    write(results,*)   'End Values'
  end subroutine printResultsVecOnNodesNDof
  
end module GIDDataOutputM
