module DataInputM
  
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM
  use GeometryM
  use SourceM
  use CFDMaterialM
  use CFDElementM
  use NormalVelocityM
  use CFDApplicationM

  use MeshM
  use ModelM
  
  implicit none
  
  private
  public :: initFEM2D
  
  integer(ikind), parameter    :: projectData = 1
  integer(ikind), parameter    :: project = 2
  integer(ikind), dimension(8) :: date_time
  integer(ikind), parameter    :: nDof = 4
  integer(ikind)               :: nElem
  integer(ikind)               :: nTriangElem
  integer(ikind)               :: nRectElem
  integer(ikind)               :: nPoint
  integer(ikind)               :: iPoint
  integer(ikind)               :: nNormalVelocity
  integer(ikind)               :: nTemperature
  integer(ikind)               :: nDensity
  integer(ikind)               :: nVelocity
  integer(ikind)               :: nMaterial
  integer(ikind)               :: nGauss
  integer(ikind)               :: isQuadratic
  integer(ikind)               :: nSourceOnPoints
  integer(ikind)               :: nSourceOnSurfaces
  integer(ikind)               :: nPointSource
  integer(ikind)               :: nSurfaceSource
  integer(ikind)               :: printStep
  real(rkind)                  :: t0
  real(rkind)                  :: errorTol
  integer(ikind)               :: maxIter
  real(rkind)                  :: fSafe
  real(rkind)                  :: constant
  real(rkind)                  :: R
  real(rkind)                  :: gamma
  real(rkind)                  :: mu
  real(rkind)                  :: k
  real(rkind)                  :: rho
  real(rkind)                  :: Vx
  real(rkind)                  :: Vy
  real(rkind)                  :: Vc
  real(rkind)                  :: T
  real(rkind)                  :: P
  real(rkind)                  :: mach
  real(rkind)                  :: Cv
  character(100)               :: projectName
  character(100)               :: aux
  logical       , parameter    :: verbose = .false.
  
  interface initFEM2D
     procedure :: initFEM2D
  end interface initFEM2D
  
contains
  
  subroutine initFEM2D(cfdAppl)
    implicit none
    type(CFDApplicationDT), intent(inout) :: cfdAppl
    print'(A)', 'Initializing CFD2D application'
    call initLog(.true., 'log.dat')
    call debugLog('  Reading project data')
    call readProjectData
    call debugLog('  Reading mesh data')
    call initMesh(cfdAppl)
    call debugLog('  Reading materials properties')
    call initMaterials(cfdAppl)
    call debugLog('  Reading elements')
    call initElements(cfdAppl)
    call debugLog('  Reading point and surface Sources')
    call readPointLineSurfaceSources(cfdAppl)
    call debugLog('  Reading Boundary Conditions')
    call readBoundaryConditions(cfdAppl)
    call debugLog('End loading data')
  end subroutine initFEM2D
  
  subroutine readProjectData
    implicit none
    open(projectData, file = 'projectData.dat')
    read(projectData, '(*(A))') projectName
    close(projectData)
    call debugLog('    Project name: ', trim(projectName))
  end subroutine readProjectData
  
  subroutine initMesh(cfdAppl)
    implicit none
    type(CFDApplicationDT), intent(inout) :: cfdAppl
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
    read(project,*)  aux, nVelocity    
    read(project,*)  aux, nDensity    
    read(project,*)  aux, nTemperature
    read(project,*)  aux, nNormalVelocity
    read(project,*)  aux, nSourceOnPoints
    read(project,*)  aux, nSourceOnSurfaces
    read(project,*)  aux, nPointSource
    read(project,*)  aux, nSurfaceSource
    read(project,*)  aux, printStep
    read(project,*)  aux, t0 
    read(project,*)  aux, errorTol 
    read(project,*)  aux, maxIter
    read(project,*)  aux, fSafe
    read(project,*)  aux, constant
    constant = 1.d0/constant
    
    if(verbose) print'(A,I0)','Number of Elements.............................: ', nElem
    if(verbose) print'(A,I0)','Are Elements Quadratic.........................: ', isQuadratic
    if(verbose) print'(A,I0)','Number of Triangular elements..................: ', nTriangElem
    if(verbose) print'(A,I0)','Number of Rectangular elements.................: ', nRectElem
    if(verbose) print'(A,I0)','Number of Nodes................................: ', nPoint
    if(verbose) print'(A,I0)','Number of Velocity conditions..................: ', nVelocity
    if(verbose) print'(A,I0)','Number of Density conditions...................: ', nDensity
    if(verbose) print'(A,I0)','Number of Temperature conditions...............: ', nTemperature   
    if(verbose) print'(A,I0)','Number of Normal Velocity conditions...........: ', nNormalVelocity 
    if(verbose) print'(A,I0)','Number of Loads on points......................: ', nSourceOnPoints 
    if(verbose) print'(A,I0)','Number of Loads on surfaces....................: ', nSourceOnSurfaces
    if(verbose) print'(A,I0)','Number of points with pointSource..............: ', nPointSource
    if(verbose) print'(A,I0)','Number of Surfaces with surfaceSource..........: ', nSurfaceSource
    if(verbose) print'(A,I0)','Number of Materials............................: ', nMaterial
    if(verbose) print'(A,I0)','Gauss cuadrature order.........................: ', nGauss
    if(verbose) print'(A,I0)','PrintStep......................................: ', printStep
    if(verbose) print'(A,I0)','Initial time...................................: ', t0
    if(verbose) print'(A,E14.7)','Error tolerance................................:', errorTol
    if(verbose) print'(A,I0)','Max Iterations.................................:', maxIter
    if(verbose) print'(A,E14.7)','Safety Factor..................................:', fSafe
    if(verbose) print'(A,E14.7)','Shock Capturing constant.......................:', constant
    
    call debugLog('    Number of Elements.............................: ', nElem)
    call debugLog('    Are Elements Quadratic.........................: ', isQuadratic)
    call debugLog('    Number of Triangular elements..................: ', nTriangElem)
    call debugLog('    Number of Rectangular elements.................: ', nRectElem)
    call debugLog('    Number of Nodes................................: ', nPoint)
    call debugLog('    Number of Velocity conditions..........---.....: ', nVelocity)
    call debugLog('    Number of Density conditions...................: ', nDensity)
    call debugLog('    Number of Temperature conditions...............: ', nTemperature)    
    call debugLog('    Number of Pressure conditions..................: ', nNormalVelocity) 
    call debugLog('    Number of Loads on points......................: ', nSourceOnPoints) 
    call debugLog('    Number of Loads on surfaces....................: ', nSourceOnSurfaces)
    call debugLog('    Number of points with pointSource..............: ', nPointSource)
    call debugLog('    Number of Surfaces with surfaceSource..........: ', nSurfaceSource)
    call debugLog('    Number of Materials............................: ', nMaterial)
    call debugLog('    Gauss cuadrature order.........................: ', nGauss)
    call debugLog('    PrintStep......................................: ', printStep)
    call debugLog('    Initial time...................................: ', t0)
    call debugLog('    Error tolerance................................: ', errorTol)
    call debugLog('    Max Iterations.................................: ', maxIter)
    call debugLog('    Safety Factor..................................: ', fSafe)
    call debugLog('    Shock Capturing constant.......................: ', constant)
    
    cfdAppl = cfdApplication(                  &
           nNode = nPoint                                      &
         , nElement = nTriangElem + nRectElem                  &
         , nNormalVelocity = nNormalVelocity                   &
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
       cfdAppl%node(iPoint) = node(iPoint, 4, x, y)
       call cfdAppl%node(iPoint)%assignDof(1, cfdAppl%model%dof(iPoint*nDof-3))
       call cfdAppl%node(iPoint)%assignDof(2, cfdAppl%model%dof(iPoint*nDof-2))
       call cfdAppl%node(iPoint)%assignDof(3, cfdAppl%model%dof(iPoint*nDof-1))
       call cfdAppl%node(iPoint)%assignDof(4, cfdAppl%model%dof(iPoint*nDof))
       call cfdAppl%model%addNode(iPoint, cfdAppl%node(iPoint))
    end do
  end subroutine initMesh

  subroutine initMaterials(cfdAppl)
    implicit none
    type(CFDApplicationDT), intent(inout) :: cfdAppl
    integer(ikind) :: i, iMat
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Material    R    gamma     mu      k     Vx_inf   Vy_inf     T     rho      Mach       Cv        P     '
    do i = 1, nMaterial
       read(project,*) iMat, R, gamma, mu, k, Vx, Vy, T, rho, mach, Cv, P
       if (T == 0.d0) T = P/(rho*R)
       if (P == 0.d0) P = rho*R*T
       if (rho == 0.d0) rho = P/(R*T)
       if ((Vx**2+Vy**2) == 0.d0) then
          Vx = sqrt(gamma*R*T)*mach
          Vy = 0.d0
       end if
       Vc = sqrt(gamma*R*T)
       call cfdAppl%setTransientValues(printStep, t0, errorTol, maxIter, fSafe, constant, R, Cv, Vc, gamma)
       cfdAppl%material(iMat) = cfdMaterial(R, gamma, mu, k, Vx, Vy, T, rho, mach, Cv, P, Vc)
       if(verbose) print'(4X,I0,7X,2(E10.3,3X))', iMat, R, gamma, mu, k, Vx, Vy, T, rho, mach, Cv, P
    end do
  end subroutine initMaterials
  
  subroutine initElements(cfdAppl)
    type(CFDApplicationDT), intent(inout) :: cfdAppl
    type(NodePtrDT), dimension(:), allocatable :: auxNode
    integer(ikind) :: i, j, iElem, iMat, nNode, Conectivities(8)
    character(len=13) :: type
    Conectivities = 0
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(A)', 'Element  |      Type      |  material index |  nNodes  |  connectivities'
    do i = 1, nElem
       read(project,*) iElem, type, iMat, nNode, (Conectivities(j),j=1,nNode)
       if(verbose) print'(I5,A15,I18,I14,5X,*(I5,X))', iElem, type, iMat, nNode, (Conectivities(j),j=1,nNode)
       allocate(auxNode(nNode))
       do j = 1, nNode
          call auxNode(j)%associate(cfdAppl%node(conectivities(j)))
       end do
       cfdAppl%element(iElem) = cfdElement(iElem, auxNode, cfdAppl%material(iMat))
       call cfdAppl%model%addElement(i, cfdAppl%element(iElem))
       deallocate(auxNode)
    end do
  end subroutine initElements
  
  subroutine readPointLineSurfaceSources(cfdAppl)
    implicit none
    type(CFDApplicationDT), intent(inout) :: cfdAppl
    integer(ikind)                        :: i, countSource, auxInt
    integer(ikind)                        :: iNode, iElem, iSource
    character(150), dimension(2)          :: func
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'nSource'
    if(verbose) print'(A)', 'Source    Function'
    do i = 1, nSourceOnPoints+nSourceOnSurfaces
       read(project,*) iSource, func(1), func(2)
       cfdAppl%source(iSource) = source(2, 2, (/'x', 'y'/), func)
       if(verbose) print'(I0,5X,30A,30A)', iSource, func(1), func(2)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'pointSources'
    if(verbose) print'(A)', 'Node    Load'
    do i = 1, nPointSource
       read(project,*) iNode, iSource
       if(verbose) print'(I0,5X,I0)', iNode, iSource
       call cfdAppl%node(iNode)%assignSource(cfdAppl%source(iSource))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'surfaceSources'
    if(verbose) print'(A)', 'Element   Load'
    do i = 1, nSurfaceSource
       read(project,*) iElem, iSource
       if(verbose) print'(I0,5X,I0)', iElem, iSource
       call cfdAppl%element(iElem)%assignSource(cfdAppl%source(iSource))
    end do
  end subroutine readPointLineSurfaceSources

  subroutine readBoundaryConditions(cfdAppl)
    implicit none
    type(CFDApplicationDT), intent(inout)        :: cfdAppl
    integer(ikind)                               :: i, j, id, elemID, nPointID
    integer(ikind)                               :: iPoint, conditionCounter
    integer(ikind), dimension(:), allocatable    :: pointID
    integer(ikind)                               :: value, valuex, valuey
    type(CFDElementDT)                           :: element
    type(NodePtrDT)  , dimension(:), allocatable :: node
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Density conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nDensity
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call cfdAppl%node(id)%fixDof(1, rho)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Velocity conditions'
    if(verbose) print'(A)', 'Node    Valuex    valuey'
    do i = 1, nVelocity
       read(Project,*) id, valuex, valuey
       if(verbose) print'(I0,5X,E10.3)', id, valuex, valuey
       call cfdAppl%node(id)%fixDof(2, rho*Vx)
       call cfdAppl%node(id)%fixDof(3, rho*Vy)
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(verbose) print'(/,A)', 'Temperature conditions'
    if(verbose) print'(A)', 'Node    Value'
    do i = 1, nTemperature
       read(Project,*) id, value
       if(verbose) print'(I0,5X,E10.3)', id, value
       call cfdAppl%node(id)%fixDof(4, rho*(P/((gamma-1.d0)*rho)+0.5d0*(Vx**2.d0+Vy**2.d0)))
    end do
    do i = 1, 7
       read(project,*)
    end do
    if(isQuadratic == 0) then
       nPointID = 2
    else if(isQuadratic == 1) then
       nPointID = 3
    end if
    allocate(pointID(nPointID))
    allocate(node(nPointID))
    conditionCounter = 0
    if(verbose) print'(/,A)', 'Normal Velocity conditions'
    if(verbose) print'(A)', 'Elem    Nodes     Value'
    do i = 1, nNormalVelocity
       conditionCounter = conditionCounter + 1
       read(Project,*) elemID, (pointID(j),j=1,nPointID), value
       if(verbose) print*, elemID, (pointID(j),j=1,nPointID), value
       element = cfdAppl%element(elemID)
       do j = 1, nPointID
          node(j) = element%node(pointID(j))
       end do
       cfdAppl%normalVelocity(i) = &
            normalVelocity(i, pointID, value, node, element%geometry, element%material)
       call cfdAppl%model%addCondition(conditionCounter, cfdAppl%normalVelocity(i))
    end do
    close(project)
  end subroutine readBoundaryConditions

end module DataInputM
  
