module DataInputM
  use UtilitiesM
  use DebuggerM

  use MeshM
  use ModelM
  
  implicit none
  
  private
  public :: initFEM2D
  
  integer(ikind), parameter    :: projectData = 1
  integer(ikind), parameter    :: project = 2
  integer(ikind), parameter    :: functions = 4
  integer(ikind), dimension(8) :: date_time
  integer(ikind)               :: nElem
  integer(ikind)               :: nTriangElem
  integer(ikind)               :: nRectElem
  integer(ikind)               :: nPoint
  integer(ikind)               :: iPoint
  integer(ikind)               :: nNormalFlux
  integer(ikind)               :: nDirichlet
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: nConvection
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nSourceOP
  integer(ikind)               :: nSourceOL
  integer(ikind)               :: nSourceOS
  integer(ikind)               :: nPointSource
  integer(ikind)               :: nLineSource
  integer(ikind)               :: nSurfaceSource
  character(100)               :: projectName
  character(100)               :: path
  character(100)               :: aux
  logical       , parameter    :: verbose = .false.
  logical                      :: isMaterialAsigned = .true.
  
  interface initFEM2D
     procedure :: initFEM2D
  end interface initFEM2D
contains
  subroutine initFEM2D(thermalAppl)
    implicit none
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(thermalAppl)
    if(isMaterialAsigned) then
       call debugLog('  Reading materials properties')
       call initMaterials(thermalAppl)
       call debugLog('  Reading elements')
       call initElements(thermalAppl)
    else
       call debugLog('  Auto asigning properties')
       call autoAsignMaterial(thermalAppl)
       call debugLog('  Reading elements')
       call initElementsDefaultMat(thermalAppl)
    end if
    call debugLog('  Reading point and line Sources')
    call readPointLineSurfaceSources(thermalAppl)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(thermalAppl)
    call debugLog('End loading data')
  end subroutine initFEM2D
  
  subroutine readProjectData
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(*(A))') projectName
    read(projectData, '(*(A))') path
    close(projectData)
    call debugLog('    Project name: ', trim(projectName))
    call debugLog('    Path: ', trim(path))
  end subroutine readProjectData
  
  subroutine initMesh(mesh)
    implicit none
    type(MeshDT), intent(inout) :: mesh
    integer(ikind) :: i
    real(rkind)    :: x, y, z
    open(project, file = trim(projectName)//'.dat')
    do i = 1, 9
       read(project,*)
    end do
    read(project,*)  aux, nElem
    read(project,*)  aux, nPoint
    read(project,*)  aux, isQuadratic
    read(project,*)  aux, nTriangElem
    read(project,*)  aux, nRectElem
    read(project,*)  aux, nMaterial
    call checknMaterial(nMaterial)
    read(project,*)  aux, nGauss    
    read(project,*)  aux, nDirichlet
    read(project,*)  aux, nNormalFlux
    read(project,*)  aux, nConvection
    read(project,*)  aux, nPointSource
    read(project,*)  aux, nSourceOP
    read(project,*)  aux, nLineSource
    read(project,*)  aux, nSourceOL
    read(project,*)  aux, nSurfaceSource
    read(project,*)  aux, nSourceOS
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of Triangular elements..................: ', nTriangElem)
    call debugLog('    Number of Rectangular elements.................: ', nRectElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Dirichlet conditions.................: ', nDirichlet)
    call debugLog('    Number of NormalFluxOnLines conditions.........: ', nNormalFlux)   
    call debugLog('    Number of ConvectionOnLines conditions.........: ', nConvection)    
    call debugLog('    Number of Sources on points....................: ', nSourceOP)
    call debugLog('    Number of points with pointSource..............: ', nPointSource)
    call debugLog('    Number of Sources on lines.....................: ', nSourceOL)
    call debugLog('    Number of points with lineSource...............: ', nLineSource)
    call debugLog('    Number of Sources on surfaces..................: ', nSourceOS)
    call debugLog('    Number of Surfaces with surfaceSource..........: ', nSurfaceSource)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)
    
    thermalAppl = thermal2DApplication(                        &
           nNode = nPoint                                      &
         , nElement = nTriangElem + nRectElem                  &
         , nCondition = nDirichlet + nNormalFlux + nConvection &
         , nSource = nSourceOP + nSourceOL + nSourceOS         &
         , nMaterial = nMaterial                               &
         , nGauss = nGauss                                     )
    
    thermalAppl%mesh = mesh(                                   &
           id = 1                                              &
         , nNode = nPoint                                      &
         , nElement = nTriangElem + nRectElem                  &
         , nCondition = nDirichlet + nNormalFlux + nConvection )
    
    call thermalAppl%initThermalProblem(nPoint, isQuadratic, nTriangElem, nRectElem, nGauss &
         , nMaterial, nPointSource, nLineSource, nSurfaceSource, nDirichlet        &
         , nNormalFlux, nConvection                                                )
    do i = 1, 6
       read(project,*)
    end do
    if(verbose) print'(A)', 'Coordinates:'
    if(verbose) print'(1X,A)', 'NODES      X          Y'
    do i = 1, nPoint
       read(project,*) iPoint, x, y
       if(verbose) print'(3X,I0,3X,E10.3,3X,E10.3,3X,E10.3)', iPoint, x, y
       thermalAppl%node(iPoint) = node(iPoint, 1, x, y)
       call mesh%addNode(iPoint, thermalAppl%node(iPoint))
    end do
  end subroutine initMesh
  
  subroutine initMaterials(thermalAppl)
    implicit none
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind) :: i, iMat
    real(rkind) :: kx, ky
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         Kx            Ky    '
    do i = 1, nMaterial
       read(project,*) iMat, kx, ky
       thermalAppl%material(iMat) = thermalMaterial(kx, ky)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X))', iMat, kx, ky 
    end do
  end subroutine initMaterials
  
  subroutine initElements(thermalAppl)
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    type(NodePtrDT), dimension(:), allocatable :: auxNode
    integer(ikind) :: i, j, iElem, iMat, nNode, Conectivities(8)
    character(len=13) :: type
    Conectivities = 0
    do i = 1, 28+nSourceOP+nSourceOL+nSourceOS
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index  |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, nNode, (Conectivities(j),j=1,nNode)
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, nNode, (Conectivities(j),j=1,nNode)
       allocate(auxNode(nNode))
       do j = 1, nNode
          auxNode(j)%ptr => thermalAppl%node(conectivities(j))
       end do
       thermalAppl%element(iElem) = thermalElement(iElem, auxNode, thermalAppl%material(iMat))
       deallocate(auxNode)
    end do
  end subroutine initElements
  
  subroutine readPointLineSurfaceSources(thermalAppl)
    implicit none
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind)                     :: i
    integer(ikind), dimension(:), allocatable :: iNode, iElem, iSource
    allocate(iNode(max(nPointSource, nLineSource)))
    allocate(iElem(nSurfaceSource))
    allocate(iSource(max(nPointSource, nLineSource, nSurfaceSource)))
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'pointSources'
    if(verbose) print'(A)', 'Node    Source'
    do i = 1, nPointSource
       read(project,*) iNode(i), iSource(i)
       if(verbose) print'(I0,5X,I0)', iNode(i), iSource(i)
       call thermalAppl%addPointSource(iNode(i), iSource(i))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'lineSources'
    if(verbose) print'(A)', 'Node   Source'
    do i = 1, nLineSource
       read(project,*) iNode(i), iSource(i)
       if(verbose) print'(I0,5X,I0)', iNode(i), iSource(i)
    end do
    if(isQuadratic == 0) then
       do i = 1, nLineSource-1
          call thermalAppl%addLineSource(iNode(i:i+1), iSource(i))
       end do
    else if(isQuadratic == 1) then
       do i = 1, nLineSource-2, 2
          call thermalAppl%addLineSource(iNode(i:i+2), iSource(i))
       end do
    end if
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'surfaceSources'
    if(verbose) print'(A)', 'Element   Source'
    do i = 1, nSurfaceSource
       read(project,*) iElem(i), iSource(i)
       if(verbose) print'(I0,5X,I0)', iElem(i), iSource(i)
       call thermalAppl%addSurfaceSource(iElem(i), iSource(i))
    end do
  end subroutine readPointLineSurfaceSources

  subroutine readBoundaryConditions(thermalAppl)
    implicit none
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind)                   :: i, j, id, elemID, nPointID, iPoint
    integer(ikind), dimension(:), allocatable :: pointID
    real(rkind)                      :: value
    real(rkind)                      :: coef, temp
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichlet
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call thermalAppl%addDirichletPoint(id, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(isQuadratic == 0) then
       nPointID = 2
    else
       nPointID = 3
    end if
    allocate(pointID(nPointID))
    if(verbose) print'(/,A)', 'Normal Flux On Points conditions'
    if(verbose) print'(A)', 'Nodes     Value'
    do i = 1, nNormalFluxOP
       read(Project,*) iPoint, value
       if(verbose) print*, iPoint, value
       call thermalAppl%addNormalFluxPoint(iPoint, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Normal Flux On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nNormalFluxOL
       read(Project,*) elemID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       call thermalAppl%addNormalFluxLine(elemID, pointID, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Convection On Points conditions'
    if(verbose) print'(A)', 'Nodes     Coef     Temp'
    do i = 1, nConvectionOP
       read(Project,*) iPoint, coef, temp
       if(verbose) print*, iPoint, coef, temp
       call thermalAppl%addConvectionPoint(iPoint, coef, temp)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Convection On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Coef     Temp'
    do i = 1, nConvectionOL
       read(Project,*) elemID, (pointID(j),j=1,nPointID), coef, temp
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), coef, temp
       call thermalAppl%addConvectionLine(elemID, pointID, coef, temp)
    end do
    close(project)
  end subroutine readBoundaryConditions

  subroutine checknMaterial(nMaterial)
    implicit none
    integer(ikind), intent(inout) :: nMaterial
    if(nMaterial == 0) then
       print*, '*************************************************************'
       print*, '**  Material not asigned, auto asigning values equal to 1  **'
       print*, '*************************************************************'
       nMaterial = 1
       isMaterialAsigned = .false.
    end if
  end subroutine checknMaterial

  subroutine autoAsignMaterial(thermalAppl)
    implicit none
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind) :: i, iMat
    real(rkind) :: kx, ky
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         Kx            Ky    '
    do i = 1, nMaterial
       iMat = 1
       kx = 1
       ky = 1
       call thermalAppl%addMaterial(kx, ky)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X),A)', iMat, kx, ky, '*AUTO ASIGNED*'
    end do
  end subroutine autoAsignMaterial

  subroutine initElementsDefaultMat(thermalAppl)
    implicit none
    type(Thermal2DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind) :: i, j, iElem, iMat, iSource, iNode, Conectivities(8), auxInt
    character(len=13) :: type
    iMat = 1
    do i = 1, 28+nSourceOP+nSourceOL+nSourceOS
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index  |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       iMat = 1
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, iNode, (Conectivities(j),j=1,iNode)
       call thermalAppl%addElement(type, iNode, iMat, Conectivities)
    end do
  end subroutine initElementsDefaultMat
end module DataInputMOD
  
