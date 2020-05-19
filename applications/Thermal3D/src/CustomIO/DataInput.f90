module DataInputM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM
  use GeometryM
  use SourceM
  use ThermalMaterialM
  use ThermalElementM
  use Thermal3DApplicationM
  use ConvectionOnSurfaceM
  use FluxOnSurfaceM

  use MeshM
  use ModelM
  
  implicit none
  
  private
  public :: initFEM3D
  
  integer(ikind), parameter    :: projectData = 1
  integer(ikind), parameter    :: project = 2
  integer(ikind), parameter    :: functions = 4
  integer(ikind), dimension(8) :: date_time
  integer(ikind)               :: nElem
  integer(ikind)               :: nTetraElem
  integer(ikind)               :: nHexaElem
  integer(ikind)               :: nPoint
  integer(ikind)               :: iPoint
  integer(ikind)               :: nNormalFlux
  integer(ikind)               :: nDirichlet
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: nConvection
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nSourceOnPoints
  integer(ikind)               :: nSourceOnVolumes
  integer(ikind)               :: nPointSource
  integer(ikind)               :: nVolumeSource
  integer(ikind)               :: nnz
  character(100)               :: projectName
  character(100)               :: path
  character(100)               :: aux
  logical       , parameter    :: verbose = .false.
  logical                      :: isMaterialAsigned = .true.
  
  interface initFEM3D
     procedure :: initFEM3D
  end interface initFEM3D
  
contains
  
  subroutine initFEM3D(thermalAppl)
    implicit none
    type(Thermal3DApplicationDT), intent(inout) :: thermalAppl
    print'(A)', 'Initializing Thermal3D application'
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(thermalAppl)
    call debugLog('  Reading materials properties')
    call initMaterials(thermalAppl)
    call debugLog('  Reading elements')
    call initElements(thermalAppl)
    call debugLog('  Reading point and line Sources')
    call readSources(thermalAppl)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(thermalAppl)
    call debugLog('End loading data')
  end subroutine initFEM3D
  
  subroutine readProjectData
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(*(A))') projectName
    read(projectData, '(*(A))') path
    close(projectData)
    call debugLog('    Project name: ', trim(projectName))
    call debugLog('    Path: ', trim(path))
  end subroutine readProjectData
  
  subroutine initMesh(thermalAppl)
    implicit none
    type(Thermal3DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind) :: i
    real(rkind)    :: x, y, z
    open(project, file = trim(projectName)//'.dat')
    do i = 1, 9
       read(project,*)
    end do
    read(project,*)  aux, nElem
    read(project,*)  aux, nPoint
    read(project,*)  aux, isQuadratic
    read(project,*)  aux, nTetraElem
    read(project,*)  aux, nHexaElem
    read(project,*)  aux, nMaterial
    read(project,*)  aux, nGauss    
    read(project,*)  aux, nDirichlet
    read(project,*)  aux, nNormalFlux
    read(project,*)  aux, nConvection
    read(project,*)  aux, nSourceOnPoints
    read(project,*)  aux, nSourceOnVolumes
    read(project,*)  aux, nPointSource
    read(project,*)  aux, nVolumeSource
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of Tetrahedral elements.................: ', nTetraElem)
    call debugLog('    Number of Hexahedral elements..................: ', nHexaElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Dirichlet conditions.................: ', nDirichlet)
    call debugLog('    Number of NormalFluxOnSurfaces conditions......: ', nNormalFlux)   
    call debugLog('    Number of ConvectionOnSurfaces conditions......: ', nConvection)    
    call debugLog('    Number of Sources on points....................: ', nSourceOnPoints) 
    call debugLog('    Number of Sources on surfaces..................: ', nSourceOnVolumes)
    call debugLog('    Number of points with pointSource..............: ', nPointSource)
    call debugLog('    Number of Volumes with volumeSource............: ', nVolumeSource)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)

    if(isQuadratic == 0) then
       nnz = nElem*8*8
    else if (isQuadratic == 1) then
       nnz = nElem*20*20
    end if
    
    thermalAppl = thermal3DApplication(                        &
           nNode = nPoint                                      &
         , nElement = nElem                                    &
         , nNormalFlux = nNormalFlux                           &
         , nConvection = nConvection                           &
         , nSource = nSourceOnPoints + nSourceOnVolumes        &
         , nMaterial = nMaterial                               &
         , nGauss = nGauss                                     &
         , nnz = nnz                                           )
    
    do i = 1, 6
       read(project,*)
    end do
    if(verbose) print'(A)', 'Coordinates:'
    if(verbose) print'(1X,A)', 'NODES      X          Y          Z'
    do i = 1, nPoint
       read(project,*) iPoint, x, y, z
       if(verbose) print'(3X,I0,3X,E10.3,3X,E10.3,3X,E10.3,3X,E10.3)', iPoint, x, y, z
       thermalAppl%node(iPoint) = node(iPoint, 1, x, y, z)
       call thermalAppl%node(iPoint)%assignDof(1, thermalAppl%model%dof(iPoint))
       call thermalAppl%model%addNode(iPoint, thermalAppl%node(iPoint))
    end do
  end subroutine initMesh
  
  subroutine initMaterials(thermalAppl)
    implicit none
    type(Thermal3DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind) :: i, iMat
    real(rkind) :: kx, ky, kz
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         Kx            Ky            Kz'
    do i = 1, nMaterial
       read(project,*) iMat, kx, ky, kz
       thermalAppl%material(iMat) = thermalMaterial(kx, ky, kz)
       if(verbose) print'(4X,I0,7X,3(E10.3,3X))', iMat, kx, ky, kz
    end do
  end subroutine initMaterials
  
  subroutine initElements(thermalAppl)
    type(Thermal3DApplicationDT), intent(inout) :: thermalAppl
    type(NodePtrDT), dimension(:), allocatable :: auxNode
    integer(ikind) :: i, j, iElem, iMat, nNode, Conectivities(20)
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
          call auxNode(j)%associate(thermalAppl%node(conectivities(j)))
       end do
       thermalAppl%element(iElem) = thermalElement(iElem, auxNode, thermalAppl%material(iMat))
       call thermalAppl%model%addElement(i, thermalAppl%element(iElem))
       deallocate(auxNode)
    end do
  end subroutine initElements
  
  subroutine readSources(thermalAppl)
    implicit none
    type(Thermal3DApplicationDT), intent(inout) :: thermalAppl
    integer(ikind)                              :: i, countSource, auxInt
    integer(ikind)                              :: iNode, iElem, iSource
    character(150), dimension(1)                :: func
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'nSource'
    if(verbose) print'(A)', 'Source    Function'
    do i = 1, nSourceOnPoints+nSourceOnVolumes
       read(project,*) iSource, func(1)
       thermalAppl%source(iSource) = source(3, 1, (/'x', 'y', 'z'/), func)
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
       call thermalAppl%node(iNode)%assignSource(thermalAppl%source(iSource))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'surfaceSources'
    if(verbose) print'(A)', 'Element   Source'
    do i = 1, nVolumeSource
       read(project,*) iElem, iSource
       if(verbose) print'(I0,5X,I0)', iElem, iSource
       call thermalAppl%element(iElem)%assignSource(thermalAppl%source(iSource))
    end do
  end subroutine readSources

  subroutine readBoundaryConditions(thermalAppl)
    implicit none
    type(Thermal3DApplicationDT), intent(inout)  :: thermalAppl
    integer(ikind)                               :: i, j, id, elemID, nPointID
    integer(ikind)                               :: iPoint, conditionCounter
    integer(ikind), dimension(:), allocatable    :: pointID
    real(rkind)                                  :: value
    real(rkind)                                  :: coef, temp
    type(ThermalElementDT)                       :: element
    type(NodePtrDT)  , dimension(:), allocatable :: node
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichlet
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call thermalAppl%node(id)%fixDof(1, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    allocate(pointID(20))
    conditionCounter = 0
    if(verbose) print'(/,A)', 'Normal Flux On Surfaces conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nNormalFlux
       conditionCounter = conditionCounter + 1
       read(Project,*) elemID, nPointID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       element = thermalAppl%element(elemID)
       allocate(node(nPointID))
       do j = 1, nPointID
          node(j) = element%node(pointID(j))
       end do
       thermalAppl%normalFluxOS(i) = fluxOnSurface(i, pointID, value, node, element%geometry)
       call thermalAppl%model%addCondition(conditionCounter, thermalAppl%normalFluxOS(i))
       deallocate(node)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Convection On Surfaces conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Coef     Temp'
    do i = 1, nConvection
       conditionCounter = conditionCounter + 1
       read(Project,*) elemID, nPointID, (pointID(j),j=1,nPointID), coef, temp
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), coef, temp
       element = thermalAppl%element(elemID)
       allocate(node(nPointID))
       do j = 1, nPointID
          node(j) = element%node(pointID(j))
       end do
       thermalAppl%convectionOS(i) = convectionOnSurface(i, pointID, coef, temp, node, element%geometry)
       call thermalAppl%model%addCondition(conditionCounter, thermalAppl%convectionOS(i))
       deallocate(node)
    end do
    close(project)
  end subroutine readBoundaryConditions

end module DataInputM
  
