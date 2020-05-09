module ThermalStructuralElementM

  use UtilitiesM

  use Triangle2D3NodeM
  use Triangle2D6NodeM
  use Quadrilateral2D4NodeM
  use Quadrilateral2D8NodeM

  use IntegratorPtrM

  use LeftHandSideM

  use PointM
  use NodeM
  use NodePtrM

  use SourceM
  use SourcePtrM
  
  use ElementM

  use StructuralMaterialM

  implicit none

  private
  public :: ThermalStructuralElementDT, thermalStructuralElement, initGeometriesTS

  type, extends(ElementDT) :: ThermalStructuralElementDT
     class(StructuralMaterialDT), pointer :: material
   contains
     procedure, public  :: init
     procedure, public  :: calculateLHS
     procedure, public  :: calculateRHS
     procedure, public  :: calculateLocalSystem
     procedure, public  :: calculateResults
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
  end type ThermalStructuralElementDT

  interface thermalStructuralElement
     procedure :: constructor
  end interface thermalStructuralElement

  type(Triangle2D3NodeDT)     , target, save :: myTriangle2D3Node
  type(Triangle2D6NodeDT)     , target, save :: myTriangle2D6Node
  type(Quadrilateral2D4NodeDT), target, save :: myQuadrilateral2D4Node
  type(Quadrilateral2D8NodeDT), target, save :: myQuadrilateral2D8Node

contains

  type(ThermalStructuralElementDT) function constructor(id, node, material)
    implicit none
    integer(ikind)                          , intent(in) :: id
    type(NodePtrDT)           , dimension(:), intent(in) :: node
    type(StructuralMaterialDT), target      , intent(in) :: material
    call constructor%init(id, node, material)
  end function constructor

  subroutine init(this, id, node, material)
    implicit none
    class(ThermalStructuralElementDT)              , intent(inout) :: this
    integer(ikind)                          , intent(in)    :: id
    type(NodePtrDT)           , dimension(:), intent(in)    :: node
    type(StructuralMaterialDT), target      , intent(in)    :: material
    this%id = id
    this%node = node
    this%material => material
    if(size(node) == 3) then
       this%geometry => myTriangle2D3Node
    else if(size(node) == 4) then
       this%geometry => myQuadrilateral2D4Node
    else if(size(node) == 6) then
       this%geometry => myTriangle2D6Node
    else if(size(node) == 8) then
       this%geometry => myQuadrilateral2D8Node
    end if
    allocate(this%source(1))
  end subroutine init

  subroutine initGeometriesTS(nGauss)
    implicit none
    integer(ikind), intent(in) :: nGauss
    myTriangle2D3Node = triangle2D3Node(nGauss)
    myTriangle2D6Node = triangle2D6Node(nGauss)
    myQuadrilateral2D4Node = quadrilateral2D4Node(nGauss)
    myQuadrilateral2D8Node = quadrilateral2D8Node(nGauss)
  end subroutine initGeometriesTS

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ThermalStructuralElementDT)                            , intent(inout) :: this
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)            , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, ii, jj, k
    integer(ikind)                                                        :: nNode, nDof
    real(rkind)                                                           :: bi, bj, ci, cj
    real(rkind)                                                           :: d11, d12, d21, d22, d33
    real(rkind)                                                           :: thickness
    real(rkind)                                                           :: temp, strain
    real(rkind)                                                           :: val1, val2
    real(rkind)            , dimension(2,2)                               :: Kij
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    real(rkind)            , dimension(:,:)  , allocatable                :: valuedSource
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    nDof = 2
    d11 = this%material%d11
    d12 = this%material%d12
    d21 = this%material%d21
    d22 = this%material%d22
    d33 = this%material%d33
    thickness = this%material%thickness
    integrator = this%getIntegrator()
    lhs = leftHandSide(0, 0, nNode*nDof)
    allocate(rhs(nNode*nDof))
    allocate(nodalPoints(nNode))
    rhs = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    do i = 1, nNode
       do j = 1, nNode
          ii = nDof*i-1
          jj = nDof*j-1
          lhs%stiffness(ii,jj)     = 0._rkind
          lhs%stiffness(ii+1,jj)   = 0._rkind
          lhs%stiffness(ii,jj+1)   = 0._rkind
          lhs%stiffness(ii+1,jj+1) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             bi = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,i) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j)
             
             Kij(1,1) = bi*bj*d11 + ci*cj*d33
             Kij(1,2) = bi*cj*d12 + bj*ci*d33
             Kij(2,1) = ci*bj*d21 + bi*cj*d33
             Kij(2,2) = bi*bj*d33 + ci*cj*d22
             
             lhs%stiffness(ii,jj)     = &
                  lhs%stiffness(ii,jj)     + integrator%getWeight(k)*Kij(1,1)/jacobianDet(k)
             lhs%stiffness(ii,jj+1)   = &
                  lhs%stiffness(ii,jj+1)   + integrator%getWeight(k)*Kij(1,2)/jacobianDet(k)
             lhs%stiffness(ii+1,jj)   = &
                  lhs%stiffness(ii+1,jj)   + integrator%getWeight(k)*Kij(2,1)/jacobianDet(k)
             lhs%stiffness(ii+1,jj+1) = &
                  lhs%stiffness(ii+1,jj+1) + integrator%getWeight(k)*Kij(2,2)/jacobianDet(k)
          end do
       end do
       if(this%node(i)%hasSource()) then
          val1 = this%node(i)%ptr%source(2)%evaluate(1, (/this%node(i)%getx(), this%node(i)%gety()/))
          val2 = this%node(i)%ptr%source(2)%evaluate(2, (/this%node(i)%getx(), this%node(i)%gety()/))
          rhs(nDof*i-1) = rhs(nDof*i-1) + val1
          rhs(nDof*i)   = rhs(nDof*i)   + val2
       end if
       temp = this%node(i)%ptr%dof(1)%val
       !Deformación plana
       strain = (1+this%material%poissonCoef)*this%material%thermalCoef*(temp-this%material%stableTemp)
       !Tensión plana
       !strain = this%material%thermalCoef*(temp-this%material%stableTemp)
       val1 = 0._rkind
       val2 = 0._rkind
       do j = 1, integrator%getIntegTerms()
          bi = jacobian(j,2,2)*integrator%getDShapeFunc(j,1,i) &
               - jacobian(j,1,2)*integrator%getDShapeFunc(j,2,i)
          ci = jacobian(j,1,1)*integrator%getDShapeFunc(j,2,i) &
               - jacobian(j,2,1)*integrator%getDShapeFunc(j,1,i)
          val1 = val1 + integrator%getWeight(j)*bi*(d11*strain+d12*strain)*thickness
          val2 = val2 + integrator%getWeight(j)*ci*(d21*strain+d22*strain)*thickness
       end do
       rhs(i*nDof-1) = rhs(i*nDof-1) + val1
       rhs(i*nDof)   = rhs(i*nDof)   + val2
    end do
    lhs%stiffness = lhs%stiffness * thickness
    if(this%hasSource()) then
       allocate(valuedSource(2,integrator%getIntegTerms()))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       do i = 1, nNode
          val1 = 0._rkind
          val2 = 0._rkind
          do j = 1, integrator%getIntegTerms()
             val1 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(1,j)*jacobianDet(j)
             val2 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(2,j)*jacobianDet(j)
          end do
          rhs(i*nDof-1) = rhs(i*nDof-1) + val1
          rhs(i*nDof)   = rhs(i*nDof)   + val2
       end do
       deallocate(valuedSource)
    end if
    deallocate(jacobian)
    deallocate(jacobianDet)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ThermalStructuralElementDT)                            , intent(inout) :: this
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, ii, jj, k
    integer(ikind)                                                        :: nNode, nDof
    real(rkind)                                                           :: bi, bj, ci, cj
    real(rkind)            , dimension(2,2)                               :: Kij
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    nDof = 2
    integrator = this%getIntegrator()
    lhs = leftHandSide(0, 0, nNode*nDof)
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    do i = 1, nNode
       do j = 1, nNode
          ii = nDof*i-1
          jj = nDof*j-1
          lhs%stiffness(ii,jj)     = 0._rkind
          lhs%stiffness(ii+1,jj)   = 0._rkind
          lhs%stiffness(ii,jj+1)   = 0._rkind
          lhs%stiffness(ii+1,jj+1) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             bi = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,i) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j)
             
             Kij(1,1) = bi*bj*this%material%d11 + ci*cj*this%material%d33
             Kij(1,2) = bi*cj*this%material%d12 + bj*ci*this%material%d33
             Kij(2,1) = ci*bj*this%material%d21 + bi*cj*this%material%d33
             Kij(2,2) = bi*bj*this%material%d33 + ci*cj*this%material%d22
             
             lhs%stiffness(ii,jj)     = &
                  lhs%stiffness(ii,jj)     + integrator%getWeight(k)*Kij(1,1)/jacobianDet(k)
             lhs%stiffness(ii,jj+1)   = &
                  lhs%stiffness(ii,jj+1)   + integrator%getWeight(k)*Kij(1,2)/jacobianDet(k)
             lhs%stiffness(ii+1,jj)   = &
                  lhs%stiffness(ii+1,jj)   + integrator%getWeight(k)*Kij(2,1)/jacobianDet(k)
             lhs%stiffness(ii+1,jj+1) = &
                  lhs%stiffness(ii+1,jj+1) + integrator%getWeight(k)*Kij(2,2)/jacobianDet(k)
          end do
       end do
    end do
    lhs%stiffness = lhs%stiffness * this%material%thickness
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ThermalStructuralElementDT)                            , intent(inout) :: this
    real(rkind)            , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, nNode, nDof
    real(rkind)                                                           :: val1, val2
    real(rkind)                                                           :: bi, ci
    real(rkind)                                                           :: d11, d12, d21, d22, d33
    real(rkind)                                                           :: thickness
    real(rkind)                                                           :: temp, strain
    real(rkind)            , dimension(:,:)  , allocatable                :: valuedSource
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    nDof = 2
    allocate(rhs(nNode*nDof))
    integrator = this%getIntegrator()
    d11 = this%material%d11
    d12 = this%material%d12
    d21 = this%material%d21
    d22 = this%material%d22
    d33 = this%material%d33
    thickness = this%material%thickness
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    rhs = 0._rkind
    do i = 1, nNode
       temp = this%node(i)%ptr%dof(1)%val
       !Deformación plana
       strain = (1+this%material%poissonCoef)*this%material%thermalCoef*(temp-this%material%stableTemp)
       !Tensión plana
       !strain = this%material%thermalCoef*(temp-this%material%stableTemp)
       val1 = 0._rkind
       val2 = 0._rkind
       do j = 1, integrator%getIntegTerms()
          bi = jacobian(j,2,2)*integrator%getDShapeFunc(j,1,i) &
               - jacobian(j,1,2)*integrator%getDShapeFunc(j,2,i)
          ci = jacobian(j,1,1)*integrator%getDShapeFunc(j,2,i) &
               - jacobian(j,2,1)*integrator%getDShapeFunc(j,1,i)
          val1 = val1 + integrator%getWeight(j)*bi*(d11*strain+d12*strain)*thickness
          val2 = val2 + integrator%getWeight(j)*ci*(d21*strain+d22*strain)*thickness
       end do
       rhs(i*nDof-1) = rhs(i*nDof-1) + val1
       rhs(i*nDof)   = rhs(i*nDof)   + val2
    end do 
    do i = 1, nNode
       if(this%node(i)%hasSource()) then
          val1 = this%node(i)%ptr%source(2)%evaluate(1, (/this%node(i)%getx(), this%node(i)%gety()/))
          val2 = this%node(i)%ptr%source(2)%evaluate(2, (/this%node(i)%getx(), this%node(i)%gety()/))
          rhs(nDof*i-1) = rhs(nDof*i-1) + val1
          rhs(nDof*i)   = rhs(nDof*i)   + val2
       end if
    end do
    if(this%hasSource()) then
       allocate(valuedSource(2,integrator%getIntegTerms()))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       do i = 1, nNode
          val1 = 0._rkind
          val2 = 0._rkind
          do j = 1, integrator%getIntegTerms()
             val1 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(1,j)*jacobianDet(j)
             val2 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(2,j)*jacobianDet(j)
          end do
          rhs(i*nDof-1) = rhs(i*nDof-1) + val1
          rhs(i*nDof)   = rhs(i*nDof)   + val2
       end do
       deallocate(valuedSource)
       deallocate(jacobianDet)
    end if   
  end subroutine calculateRHS

  subroutine setupIntegration(this, integrator, valuedSource, jacobianDet)
    implicit none
    class(ThermalStructuralElementDT)                          , intent(inout) :: this
    type(IntegratorPtrDT)                               , intent(in)    :: integrator
    real(rkind), dimension(2,integrator%getIntegTerms()), intent(out)   :: valuedSource
    real(rkind), dimension(integrator%getIntegTerms())  , intent(out)   :: jacobianDet
    integer(ikind)                                                      :: i, nNode
    type(NodePtrDT), dimension(:), allocatable                          :: nodalPoints
    nNode = this%getnNode()
    allocate(nodalPoints(nNode))
    valuedSource = this%getValuedSource(integrator)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobianDet = this%geometry%jacobianDetAtGPoints(nodalPoints)
  end subroutine setupIntegration

  function getValuedSource(this, integrator)
    implicit none
    class(ThermalStructuralElementDT), intent(inout) :: this
    type(IntegratorPtrDT) , intent(in) :: integrator
    real(rkind), dimension(2,integrator%getIntegTerms()) :: getValuedSource
    integer(ikind) :: i, j, nNode
    real(rkind) :: x, y
    type(NodePtrDT), dimension(:), allocatable :: node
    nNode = this%getnNode()
    do i = 1, integrator%getIntegTerms()
       node = this%node
       x = 0
       y = 0
       do j = 1, nNode
          x = x + integrator%getShapeFunc(i,j)*node(j)%ptr%getx()
          y = y + integrator%getShapeFunc(i,j)*node(j)%ptr%gety()
       end do
       getValuedSource(1,i) = this%source(1)%evaluate(1, (/x,y/))
       getValuedSource(2,i) = this%source(1)%evaluate(2, (/x,y/))
    end do
  end function getValuedSource

  subroutine calculateResults(this, resultMat)
    implicit none
    class(ThermalStructuralElementDT)                     , intent(inout) :: this
    real(rkind)            , dimension(:,:,:), allocatable, intent(inout) :: resultMat
    integer(ikind)                                                        :: i, iGauss, nNode
    real(rkind)                                                           :: nsx, nsy !NormalStress
    real(rkind)                                                           :: shs      !SheatStress
    real(rkind)                                                           :: epx, epy !Strain
    real(rkind)                                                           :: bi, ci, xi, eta
    real(rkind)                                                           :: dNidx, dNidy, kx, ky
    real(rkind)                                                           :: d11, d12, d21, d22, d33
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    integrator = this%getIntegrator()
    nNode = this%getnNode()
    allocate(nodalPoints(nNode))
    allocate(resultMat(3,integrator%getIntegTerms(),2))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    do iGauss = 1, integrator%getIntegTerms()
       xi = integrator%getGPoint(iGauss,1)
       eta = integrator%getGPoint(iGauss,2)
       nsx = 0._rkind
       nsy = 0._rkind
       shs = 0._rkind
       epx = 0._rkind
       epy = 0._rkind
       do i = 1, nNode
          bi = jacobian(iGauss,2,2)*integrator%getDShapeFunc(iGauss,1,i) &
               - jacobian(iGauss,1,2)*integrator%getDShapeFunc(iGauss,2,i)
          ci = jacobian(iGauss,1,1)*integrator%getDShapeFunc(iGauss,2,i) &
               - jacobian(iGauss,2,1)*integrator%getDShapeFunc(iGauss,1,i)
          dNidx = bi/jacobianDet(iGauss)
          dNidy = ci/jacobianDet(iGauss)
          nsx = nsx + dNidx*this%node(i)%ptr%dof(2)%val
          nsy = nsy + dNidy*this%node(i)%ptr%dof(3)%val
          shs = shs + dNidx*this%node(i)%ptr%dof(3)%val + dNidy*this%node(i)%ptr%dof(2)%val
          epx = epx + dNidx*this%node(i)%ptr%dof(2)%val
          epy = epy + dNidy*this%node(i)%ptr%dof(3)%val
       end do
       d11 = this%material%d11
       d12 = this%material%d12
       d21 = this%material%d21
       d22 = this%material%d22
       d33 = this%material%d33
       resultMat(1,iGauss,1) = d11*nsx+d12*nsy
       resultMat(1,iGauss,2) = d21*nsx+d22*nsy
       resultMat(2,iGauss,1) = d33*shs
       resultMat(3,iGauss,1) = epx
       resultMat(3,iGauss,2) = epy
    end do
  end subroutine calculateResults

end module ThermalStructuralElementM
