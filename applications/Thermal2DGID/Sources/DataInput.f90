module DataInputMOD
  use DebuggerMOD
  use IODataMOD
  implicit none
  private
  public :: initFEM2D
  logical                      :: isMaterialAsigned = .true.
  interface initFEM2D
     procedure :: initFEM2D
  end interface initFEM2D
contains
  subroutine initFEM2D(io)
    implicit none
    type(IODataTYPE), intent(inout) :: io
    call io%setUp(io)
      print'(A)', 'Initiating Thermal2D'
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(io)
    if(isMaterialAsigned) then
       call debugLog('  Reading materials properties')
       call initMaterials(io)
       call debugLog('  Reading elements')
       call initElements(io)
    else
       call debugLog('  Auto asigning properties')
       call autoAsignMaterial(io)
       call debugLog('  Reading elements')
       call initElementsDefaultMat(io)
    end if
    call debugLog('  Reading Sources')
    call readPointLineSurfaceSources(io)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(io)
    call debugLog('End loading data')
  end subroutine initFEM2D
end module DataInputMOD


