module DataInputM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM
  use GeometryM
  use SourceM
  use StructuralMaterialM
  use StructuralElementM
  use PressureM
  use Structural3DApplicationM

  use MeshM
  use ModelM
  
  implicit none
  
  private
  public :: initFEM3D
  
  integer(ikind), parameter    :: projectData = 1
  integer(ikind), parameter    :: project = 2
  integer(ikind), parameter    :: functions = 4
  integer(ikind), dimension(8) :: date_time
  integer(ikind), parameter    :: nDof = 3
  integer(ikind)               :: nElem
  integer(ikind)               :: nTetraElem
  integer(ikind)               :: nHexaElem
  integer(ikind)               :: nPoint
  integer(ikind)               :: iPoint
  integer(ikind)               :: nPressure
  integer(ikind)               :: nDirichletX
  integer(ikind)               :: nDirichletY
  integer(ikind)               :: nDirichletZ
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: nConvection
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nSourceOnPoints
  integer(ikind)               :: nSourceOnVolumes
  integer(ikind)               :: nPointSource
  integer(ikind)               :: nVolumeSource
  character(100)               :: projectName
  character(100)               :: aux
  logical       , parameter    :: verbose = .false.
  logical                      :: isMaterialAsigned = .true.
  
  interface initFEM3D
     procedure :: initFEM3D
  end interface initFEM3D
  
contains
  
  subroutine initFEM3D(structuralAppl)
    implicit none
    type(Structural3DApplicationDT), intent(inout) :: structuralAppl
    print'(A)', 'Initializing Structural3D application'
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(structuralAppl)
    call debugLog('  Reading materials properties')
    call initMaterials(structuralAppl)
    call debugLog('  Reading elements')
    call initElements(structuralAppl)
    call debugLog('  Reading point and line Sources')
    call readSources(structuralAppl)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(structuralAppl)
    call debugLog('End loading data')
  end subroutine initFEM3D
  
  subroutine readProjectData
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(*(A))') projectName
    close(projectData)
    call debugLog('    Project name: ', trim(projectName))
  end subroutine readProjectData
  
  subroutine initMesh(structuralAppl)
    implicit none
    type(Structural3DApplicationDT), intent(inout) :: structuralAppl
    integer(ikind) :: i
    integer(ikind) :: nnz
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
    read(project,*)  aux, nDirichletX
    read(project,*)  aux, nDirichletY
    read(project,*)  aux, nDirichletZ
    read(project,*)  aux, nPressure
    read(project,*)  aux, nSourceOnPoints
    read(project,*)  aux, nSourceOnVolumes
    read(project,*)  aux, nPointSource
    read(project,*)  aux, nVolumeSource
    
    if(verbose) print'(A,I0)','Number of Elements.............................: ', nElem
    if(verbose) print'(A,I0)','Are Elements Quadratic.........................: ', isQuadratic
    if(verbose) print'(A,I0)','Number of Tetrahedral elements.................: ', nTetraElem
    if(verbose) print'(A,I0)','Number of Hexahedral elements..................: ', nHexaElem
    if(verbose) print'(A,I0)','Number of Nodes................................: ', nPoint
    if(verbose) print'(A,I0)','Number of Dirichlet X conditions...............: ', nDirichletX   
    if(verbose) print'(A,I0)','Number of Dirichlet Y conditions...............: ', nDirichletY   
    if(verbose) print'(A,I0)','Number of Dirichlet Z conditions...............: ', nDirichletZ     
    if(verbose) print'(A,I0)','Number of Pressure conditions..................: ', nPressure 
    if(verbose) print'(A,I0)','Number of Loads on points......................: ', nSourceOnPoints 
    if(verbose) print'(A,I0)','Number of Loads on volumes.....................: ', nSourceOnVolumes
    if(verbose) print'(A,I0)','Number of points with pointSource..............: ', nPointSource
    if(verbose) print'(A,I0)','Number of volumes with volumeSource............: ', nVolumeSource
    if(verbose) print'(A,I0)','Number of Materials............................: ', nMaterial
    if(verbose) print'(A,I0)','Gauss cuadrature order.........................: ', nGauss
    
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of Tetrahedral elements.................: ', nTetraElem)
    call debugLog('    Number of Hexahedral elements..................: ', nHexaElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Dirichlet X conditions...............: ', nDirichletX)    
    call debugLog('    Number of Dirichlet Y conditions...............: ', nDirichletY)    
    call debugLog('    Number of Dirichlet Z conditions...............: ', nDirichletZ)      
    call debugLog('    Number of Pressure conditions..................: ', nPressure)    
    call debugLog('    Number of Sources on points....................: ', nSourceOnPoints) 
    call debugLog('    Number of Sources on surfaces..................: ', nSourceOnVolumes)
    call debugLog('    Number of points with pointSource..............: ', nPointSource)
    call debugLog('    Number of Volumes with volumeSource............: ', nVolumeSource)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)

    if(isQuadratic == 0) then
       nnz = nElem*8*8*3*3
    else if (isQuadratic == 1) then
       nnz = nElem*20*20*3*3
    end if
    
    structuralAppl = structural3DApplication(                  &
           nNode = nPoint                                      &
         , nElement = nElem                                    &
         , nPressure = nPressure                               &
         , nSource = nSourceOnPoints + nSourceOnVolumes        &
         , nMaterial = nMaterial                               &
         , nGauss = nGauss                                     &
         , nnz = nnz                                           )
    
    do i = 1, 6
       read(project,*)
    end do
    if(verbose) print'(A)', 'Coordinates:'
    if(verbose) print'(1X,A)', 'NODES      X          Y'
    do i = 1, nPoint
       read(project,*) iPoint, x, y, z
       if(verbose) print'(3X,I0,3X,E10.3,3X,E10.3,3X,E10.3)', iPoint, x, y, z
       structuralAppl%node(iPoint) = node(iPoint, 3, x, y, z)
       call structuralAppl%node(iPoint)%assignDof(1, structuralAppl%model%dof(iPoint*nDof-2))
       call structuralAppl%node(iPoint)%assignDof(2, structuralAppl%model%dof(iPoint*nDof-1))
       call structuralAppl%node(iPoint)%assignDof(3, structuralAppl%model%dof(iPoint*nDof))
       call structuralAppl%model%addNode(iPoint, structuralAppl%node(iPoint))
    end do
  end subroutine initMesh
  
  subroutine initMaterials(structuralAppl)
    implicit none
    type(Structural3DApplicationDT), intent(inout) :: structuralAppl
    integer(ikind) :: i, iMat
    real(rkind)    :: alpha, E, nu
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material         alpha        E        nu'
    do i = 1, nMaterial
       read(project,*) iMat, alpha, E, nu
       structuralAppl%material(iMat) = structuralMaterial(E, nu, alpha)
       if(verbose) print'(4X,I0,7X,5(E10.3,3X))', iMat, alpha, E, nu 
    end do
  end subroutine initMaterials
  
  subroutine initElements(structuralAppl)
    type(Structural3DApplicationDT), intent(inout) :: structuralAppl
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
          call auxNode(j)%associate(structuralAppl%node(conectivities(j)))
       end do
       structuralAppl%element(iElem) = structuralElement(iElem, auxNode, structuralAppl%material(iMat))
       call structuralAppl%model%addElement(i, structuralAppl%element(iElem))
       deallocate(auxNode)
    end do
  end subroutine initElements
  
  subroutine readSources(structuralAppl)
    implicit none
    type(Structural3DApplicationDT), intent(inout) :: structuralAppl
    integer(ikind)                              :: i, countSource, auxInt
    integer(ikind)                              :: iNode, iElem, iSource
    character(150), dimension(3)                :: func
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'nSource'
    if(verbose) print'(A)', 'Source    Function'
    do i = 1, nSourceOnPoints+nSourceOnVolumes
       read(project,*) iSource, func(1), func(2), func(3)
       structuralAppl%source(iSource) = source(3, 3, (/'x', 'y', 'z'/), func)
       if(verbose) print'(I0,5X,30A,30A,30A)', iSource, func(1), func(2), func(3)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'pointSources'
    if(verbose) print'(A)', 'Node    Load'
    do i = 1, nPointSource
       read(project,*) iNode, iSource
       if(verbose) print'(I0,5X,I0)', iNode, iSource
       call structuralAppl%node(iNode)%assignSource(structuralAppl%source(iSource))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'volumeSources'
    if(verbose) print'(A)', 'Element   Load'
    do i = 1, nVolumeSource
       read(project,*) iElem, iSource
       if(verbose) print'(I0,5X,I0)', iElem, iSource
       call structuralAppl%element(iElem)%assignSource(structuralAppl%source(iSource))
    end do
  end subroutine readSources

  subroutine readBoundaryConditions(structuralAppl)
    implicit none
    type(Structural3DApplicationDT), intent(inout)  :: structuralAppl
    integer(ikind)                               :: i, j, id, elemID, nPointID
    integer(ikind)                               :: iPoint, conditionCounter
    integer(ikind), dimension(:), allocatable    :: pointID
    real(rkind)                                  :: value
    real(rkind)                                  :: coef, temp
    type(StructuralElementDT)                    :: element
    type(NodePtrDT)  , dimension(:), allocatable :: node
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet X conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichletX
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call structuralAppl%node(id)%fixDof(1, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet Y conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichletY
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call structuralAppl%node(id)%fixDof(2, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Dirichlet Z conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDirichletZ
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call structuralAppl%node(id)%fixDof(3, value)
    end do
    do i = 1, 7
       read(project,*)
    end do
    allocate(pointID(20))
    conditionCounter = 0
    if(verbose) print'(/,A)', 'Pressure On Surfaces conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nPressure
       conditionCounter = conditionCounter + 1
       read(Project,*) elemID, nPointID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       element = structuralAppl%element(elemID)
       allocate(node(nPointID))
       do j = 1, nPointID
          node(j) = element%node(pointID(j))
       end do
       structuralAppl%pressure(i) = &
            pressure(i, pointID, value, node, element%geometry, element%material)
       call structuralAppl%model%addCondition(conditionCounter, structuralAppl%pressure(i))
       deallocate(node)
    end do
    close(project)
  end subroutine readBoundaryConditions

end module DataInputM
  
