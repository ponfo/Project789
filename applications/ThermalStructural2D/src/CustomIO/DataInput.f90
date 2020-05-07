module DataInputM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM
  use GeometryM
  use SourceM
  use ThermalMaterialM
  use ThermalElementM
  use Thermal2DApplicationM
  use ConvectionOnLineM
  use FluxOnLineM
  use PressumeM

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
  integer(ikind)               :: nFixTemp
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: nConvection
  integer(ikind)               :: nPressure
  integer(ikind)               :: nFixDisplX
  integer(ikind)               :: nFixDisplY
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nSourceOnPoints
  integer(ikind)               :: nSourceOnSurfaces
  integer(ikind)               :: nPointSource
  integer(ikind)               :: nSurfaceSource
  integer(ikind)               :: nLoadOnPoints
  integer(ikind)               :: nLoadOnSurfaces
  integer(ikind)               :: nPointLoad
  integer(ikind)               :: nSurfaceLoad
  character(100)               :: projectName
  character(100)               :: path
  character(100)               :: aux
  logical       , parameter    :: verbose = .false.
  logical                      :: isMaterialAsigned = .true.
  
  interface initFEM2D
     procedure :: initFEM2D
  end interface initFEM2D
  
contains
  
  subroutine initFEM2D(thermalStructAppl)
    implicit none
    type(ThermalStruct2DApplicationDT), intent(inout) :: thermalStructAppl
    print'(A)', 'Initializing ThermalStruct2D application'
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(thermalStructAppl)
    call debugLog('  Reading materials properties')
    call initMaterials(thermalStructAppl)
    call debugLog('  Reading elements')
    call initElements(thermalStructAppl)
    call debugLog('  Reading point and line Sources')
    call readPointLineSurfaceSources(thermalStructAppl)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(thermalStructAppl)
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
  
  subroutine initMesh(thermalStructAppl)
    implicit none
    type(ThermalStruct2DApplicationDT), intent(inout) :: thermalStructAppl
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
    read(project,*)  aux, nGauss    
    read(project,*)  aux, nFixTemp    
    read(project,*)  aux, nFixDisplX
    read(project,*)  aux, nFixDisplY
    read(project,*)  aux, nNormalFlux
    read(project,*)  aux, nConvection
    read(project,*)  aux, nPressure
    read(project,*)  aux, nSourceOnPoints
    read(project,*)  aux, nSourceOnSurfaces
    read(project,*)  aux, nPointSource
    read(project,*)  aux, nSurfaceSource
    read(project,*)  aux, nLoadOnPoints
    read(project,*)  aux, nLoadOnSurfaces
    read(project,*)  aux, nPointLoad
    read(project,*)  aux, nSurfaceLoad
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of Triangular elements..................: ', nTriangElem)
    call debugLog('    Number of Rectangular elements.................: ', nRectElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Fix Temperature conditions...........: ', nFixTemp)
    call debugLog('    Number of Fix Displacement X conditions........: ', nFixDisplX)
    call debugLog('    Number of Fix Displacement Y conditions........: ', nFixDisplY)
    call debugLog('    Number of NormalFluxOnLines conditions.........: ', nNormalFlux)   
    call debugLog('    Number of ConvectionOnLines conditions.........: ', nConvection)     
    call debugLog('    Number of Pressure conditions..................: ', nPressure)    
    call debugLog('    Number of Sources on points....................: ', nSourceOnPoints) 
    call debugLog('    Number of Sources on surfaces..................: ', nSourceOnSurfaces)
    call debugLog('    Number of points with pointSource..............: ', nPointSource)
    call debugLog('    Number of Surfaces with surfaceSource..........: ', nSurfaceSource)   
    call debugLog('    Number of Loads on points......................: ', nLoadOnPoints) 
    call debugLog('    Number of Loads on surfaces....................: ', nLoadOnSurfaces)
    call debugLog('    Number of points with pointLoad................: ', nPointLoad)
    call debugLog('    Number of Surfaces with surfaceLoad............: ', nSurfaceLoad)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)
    
    thermalStructAppl = thermalStruct2DApplication(            &
           nNode = nPoint                                      &
         , nElement = nTriangElem + nRectElem                  &
         , nNormalFlux = nNormalFlux                           &
         , nConvection = nConvection                           &
         , nPressure = nPressure                               &
         , nSource = nSourceOnPoints + nSourceOnSurfaces       &
         , nMaterial = nMaterial                               &
         , nGauss = nGauss                                     )
    
    do i = 1, 6
       read(project,*)
    end do
    if(verbose) print'(A)', 'Coordinates:'
    if(verbose) print'(1X,A)', 'NODES      X          Y'
    do i = 1, nPoint
       read(project,*) iPoint, x, y
       if(verbose) print'(3X,I0,3X,E10.3,3X,E10.3,3X,E10.3)', iPoint, x, y
       thermalStructAppl%node(iPoint) = node(iPoint, 1, x, y)
       call thermalStructAppl%node(iPoint)%assignDof(1, thermalStructAppl%model%dof(iPoint))
       call thermalStructAppl%model%addNode(iPoint, thermalStructAppl%node(iPoint))
    end do
  end subroutine initMesh
  
  subroutine initMaterials(thermalStructAppl)
    implicit none
    type(ThermalStruct2DApplicationDT), intent(inout) :: thermalStructAppl
    integer(ikind) :: i, iMat
    real(rkind) :: kx, ky
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         Kx            Ky    '
    do i = 1, nMaterial
       read(project,*) iMat, kx, ky
       thermalStructAppl%material(iMat) = thermalStructMaterial(kx, ky)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X))', iMat, kx, ky 
    end do
  end subroutine initMaterials
  
  subroutine initElements(thermalStructAppl)
    type(ThermalStruct2DApplicationDT), intent(inout) :: thermalStructAppl
    type(NodePtrDT), dimension(:), allocatable :: auxNode
    integer(ikind) :: i, j, iElem, iMat, nNode, Conectivities(8)
    character(len=13) :: type
    Conectivities = 0
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index  |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, nNode, (Conectivities(j),j=1,nNode)
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, nNode, (Conectivities(j),j=1,nNode)
       allocate(auxNode(nNode))
       do j = 1, nNode
          call auxNode(j)%associate(thermalStructAppl%node(conectivities(j)))
       end do
       thermalStructAppl%element(iElem) = thermalStructElement(iElem, auxNode, thermalStructAppl%material(iMat))
       call thermalStructAppl%model%addElement(i, thermalStructAppl%element(iElem))
       deallocate(auxNode)
    end do
  end subroutine initElements
  
  subroutine readPointLineSurfaceSources(thermalStructAppl)
    implicit none
    type(ThermalStruct2DApplicationDT), intent(inout) :: thermalStructAppl
    integer(ikind)                              :: i, countSource, auxInt
    integer(ikind)                              :: iNode, iElem, iSource
    character(150), dimension(1)                :: func
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'nSource'
    if(verbose) print'(A)', 'Source    Function'
    do i = 1, nSourceOnPoints+nSourceOnSurfaces
       read(project,*) iSource, func(1)
       thermalStructAppl%source(iSource) = source(2, 1, (/'x', 'y'/), func)
       if(verbose) print'(I0,5X,A)', iSource, func(1)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'pointSources'
    if(verbose) print'(A)', 'Node    Source'
    do i = 1, nPointSource
       read(project,*) iNode, iSource
       if(verbose) print'(I0,5X,I0)', iNode, iSource
       call thermalStructAppl%node(iNode)%assignSource(thermalStructAppl%source(iSource))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'surfaceSources'
    if(verbose) print'(A)', 'Element   Source'
    do i = 1, nSurfaceSource
       read(project,*) iElem, iSource
       if(verbose) print'(I0,5X,I0)', iElem, iSource
       call thermalStructAppl%element(iElem)%assignSource(thermalStructAppl%source(iSource))
    end do
  end subroutine readPointLineSurfaceSources

  subroutine readBoundaryConditions(thermalStructAppl)
    implicit none
    type(ThermalStruct2DApplicationDT), intent(inout)  :: thermalStructAppl
    integer(ikind)                               :: i, j, id, elemID, nPointID
    integer(ikind)                               :: iPoint, conditionCounter
    integer(ikind), dimension(:), allocatable    :: pointID
    real(rkind)                                  :: value
    real(rkind)                                  :: coef, temp
    type(ThermalStructElementDT)                       :: element
    type(NodePtrDT)  , dimension(:), allocatable :: node
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichlet
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call thermalStructAppl%node(id)%fixDof(1, value)
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
    allocate(node(nPointID))
    conditionCounter = 0
    if(verbose) print'(/,A)', 'Normal Flux On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nNormalFlux
       conditionCounter = conditionCounter + 1
       read(Project,*) elemID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       element = thermalStructAppl%element(elemID)
       do j = 1, nPointID
          node(j) = element%node(pointID(j))
       end do
       thermalStructAppl%normalFluxOL(i) = fluxOnLine(i, pointID, value, node, element%geometry)
       call thermalStructAppl%model%addCondition(conditionCounter, thermalStructAppl%normalFluxOL(i))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Convection On Lines conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Coef     Temp'
    do i = 1, nConvection
       conditionCounter = conditionCounter + 1
       read(Project,*) elemID, (pointID(j),j=1,nPointID), coef, temp
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), coef, temp
       element = thermalStructAppl%element(elemID)
       do j = 1, nPointID
          node(j) = element%node(pointID(j))
       end do
       thermalStructAppl%convectionOL(i) = convectionOnLine(i, pointID, coef, temp, node, element%geometry)
       call thermalStructAppl%model%addCondition(conditionCounter, thermalStructAppl%convectionOL(i))
    end do
    close(project)
  end subroutine readBoundaryConditions

end module DataInputM
  