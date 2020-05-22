module CFDElementM

  use UtilitiesM
  use SparseKit

  use Triangle2D3NodeM
  use Triangle2D6NodeM
  use Quadrilateral2D4NodeM
  use Quadrilateral2D8NodeM

  use IntegratorPtrM

  use LeftHandSideM
  use ProcessInfoM

  use PointM
  use NodeM
  use NodePtrM

  use SourceM
  use SourcePtrM

  use ElementM
  
  use CFDMaterialM

  implicit none

  private
  public :: CFDElementDT, cfdElement, initGeometries

  type, extends(ElementDT) :: CFDElementDT
     class(ThermalMaterialDT), pointer :: material
   contains
     procedure, public  :: init
     procedure, public  :: calculateLHS
     procedure, public  :: calculateRHS
     procedure, public  :: calculateLocalSystem
     procedure, public  :: calculateSystem
     procedure, public  :: calculateResults
     procedure, public  :: calculateDT
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
  end type CFDElementDT

  interface cfdElement
     procedure :: constructor
  end interface cfdElement

  type(Triangle2D3NodeDT)     , target, save :: myTriangle2D3Node
  type(Triangle2D6NodeDT)     , target, save :: myTriangle2D6Node
  type(Quadrilateral2D4NodeDT), target, save :: myQuadrilateral2D4Node
  type(Quadrilateral2D8NodeDT), target, save :: myQuadrilateral2D8Node

contains

  type(CFDElementDT) function constructor(id, node, material)
    implicit none
    integer(ikind)                          , intent(in) :: id
    type(NodePtrDT)           , dimension(:), intent(in) :: node
    class(ThermalMaterialDT), target      , intent(in) :: material
    call constructor%init(id, node, material)
  end function constructor

  subroutine init(this, id, node, material)
    implicit none
    class(CFDElementDT)                     , intent(inout) :: this
    integer(ikind)                          , intent(in)    :: id
    type(NodePtrDT)           , dimension(:), intent(in)    :: node
    class(ThermalMaterialDT), target        , intent(in)    :: material
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

  subroutine initGeometries(nGauss)
    implicit none
    integer(ikind), intent(in) :: nGauss
    myTriangle2D3Node = triangle2D3Node(nGauss)
    myTriangle2D6Node = triangle2D6Node(nGauss)
    myQuadrilateral2D4Node = quadrilateral2D4Node(nGauss)
    myQuadrilateral2D8Node = quadrilateral2D8Node(nGauss)
  end subroutine initGeometries

  subroutine calculateLocalSystem(this, processInfo, lhs, rhs)
    implicit none
    class(CFDElementDT)                                   , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)            , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, ii, jj, k
    integer(ikind)                                                        :: nNode, nDof
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    real(rkind)                                                           :: val1, val2
    real(rkind)                                                           :: val3, val4
    real(rkind)            , dimension(:,:)  , allocatable                :: valuedSource
    real(rkind)                                                           :: dt
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    nDof = this%node(1)%getnDof()
    integrator = this%getIntegrator()
    lhs = leftHandSide(nNode*nDof, 0, nNode*nDof)
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
          ii   = nDof*i-3
          jj   = nDof*j-3
          lhs%mass(ii  ,jj  ) = 0._rkind
          lhs%mass(ii  ,jj+1) = 0._rkind
          lhs%mass(ii  ,jj+2) = 0._rkind
          lhs%mass(ii  ,jj+3) = 0._rkind

          lhs%mass(ii+1,jj  ) = 0._rkind
          lhs%mass(ii+1,jj+1) = 0._rkind
          lhs%mass(ii+1,jj+2) = 0._rkind
          lhs%mass(ii+1,jj+3) = 0._rkind

          lhs%mass(ii+2,jj  ) = 0._rkind
          lhs%mass(ii+2,jj+1) = 0._rkind
          lhs%mass(ii+2,jj+2) = 0._rkind
          lhs%mass(ii+2,jj+3) = 0._rkind

          lhs%mass(ii+3,jj  ) = 0._rkind
          lhs%mass(ii+3,jj+1) = 0._rkind
          lhs%mass(ii+3,jj+2) = 0._rkind
          lhs%mass(ii+3,jj+3) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             lhs%mass(ii,jj  )   = &
                  lhs%mass(ii,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii,jj+1)   = &
                  lhs%mass(ii,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii,jj+2)   = &
                  lhs%mass(ii,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii,jj+3)   = &
                  lhs%mass(ii,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))

             lhs%mass(ii+1,jj)     = &
                  lhs%mass(ii+1,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+1,jj+1)   = &
                  lhs%mass(ii+1,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+1,jj+2)   = &
                  lhs%mass(ii+1,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+1,jj+32)  = &
                  lhs%mass(ii+1,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))

             lhs%mass(ii+2,jj)     = &
                  lhs%mass(ii+2,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+2,jj+1)   = &
                  lhs%mass(ii+2,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+2,jj+2)   = &
                  lhs%mass(ii+2,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+2,jj+3)   = &
                  lhs%mass(ii+2,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))

             lhs%mass(ii+3,jj)     = &
                  lhs%mass(ii+3,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+3,jj+1)   = &
                  lhs%mass(ii+3,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+3,jj+2)   = &
                  lhs%mass(ii+3,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+3,jj+3)   = &
                  lhs%mass(ii+3,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j)) 
          end do
       end do
       if(this%node(i)%hasSource()) then
          val1 = this%node(i)%ptr%source(1)%evaluate(1, (/this%node(i)%getx(), this%node(i)%gety()/))
          val2 = this%node(i)%ptr%source(1)%evaluate(2, (/this%node(i)%getx(), this%node(i)%gety()/))
          val3 = this%node(i)%ptr%source(1)%evaluate(3, (/this%node(i)%getx(), this%node(i)%gety()/))
          val4 = this%node(i)%ptr%source(1)%evaluate(4, (/this%node(i)%getx(), this%node(i)%gety()/))
          rhs(nDof*i-3) = rhs(nDof*i-3) + val1
          rhs(nDof*i-2) = rhs(nDof*i-2) + val2
          rhs(nDof*i-1) = rhs(nDof*i-1) + val3
          rhs(nDof*i)   = rhs(nDof*i)   + val4
       end if
    end do
    if(this%hasSource()) then
       allocate(valuedSource(4,integrator%getIntegTerms()))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       do i = 1, nNode
          val1 = 0._rkind
          val2 = 0._rkind
          val3 = 0._rkind
          val4 = 0._rkind
          do j = 1, integrator%getIntegTerms()
             val1 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(1,j)*jacobianDet(j)
             val2 = val2 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(2,j)*jacobianDet(j)
             val3 = val3 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(3,j)*jacobianDet(j)
             val4 = val4 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(4,j)*jacobianDet(j)
          end do
          rhs(i*nDof-3) = rhs(i*nDof-3) + val1
          rhs(i*nDof-2) = rhs(i*nDof-2) + val2
          rhs(i*nDof-1) = rhs(i*nDof-1) + val3
          rhs(i*nDof)   = rhs(i*nDof)   + val4
       end do
       deallocate(valuedSource)
    end if
    deallocate(jacobian)
    deallocate(jacobianDet)
    call this%calculateSystem(processInfo, lhs)
    call this%calculateDT(processInfo, dt)
    call processInfo%setMinimumDT(dt)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, processInfo, lhs)
    implicit none
    class(CFDElementDT)                                   , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, ii, jj, k
    integer(ikind)                                                        :: nNode, nDof
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    nDof = this%node(1)%getnDof()
    integrator = this%getIntegrator()
    lhs = leftHandSide(nNode*nDof, 0, 0)
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    do i = 1, nNode
       do j = 1, nNode
          ii   = nDof*i-3
          jj   = nDof*j-3
          lhs%mass(ii  ,jj  ) = 0._rkind
          lhs%mass(ii  ,jj+1) = 0._rkind
          lhs%mass(ii  ,jj+2) = 0._rkind
          lhs%mass(ii  ,jj+3) = 0._rkind

          lhs%mass(ii+1,jj  ) = 0._rkind
          lhs%mass(ii+1,jj+1) = 0._rkind
          lhs%mass(ii+1,jj+2) = 0._rkind
          lhs%mass(ii+1,jj+3) = 0._rkind

          lhs%mass(ii+2,jj  ) = 0._rkind
          lhs%mass(ii+2,jj+1) = 0._rkind
          lhs%mass(ii+2,jj+2) = 0._rkind
          lhs%mass(ii+2,jj+3) = 0._rkind

          lhs%mass(ii+3,jj  ) = 0._rkind
          lhs%mass(ii+3,jj+1) = 0._rkind
          lhs%mass(ii+3,jj+2) = 0._rkind
          lhs%mass(ii+3,jj+3) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             lhs%mass(ii,jj  )   = &
                  lhs%mass(ii,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii,jj+1)   = &
                  lhs%mass(ii,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii,jj+2)   = &
                  lhs%mass(ii,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii,jj+3)   = &
                  lhs%mass(ii,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))

             lhs%mass(ii+1,jj)     = &
                  lhs%mass(ii+1,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+1,jj+1)   = &
                  lhs%mass(ii+1,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+1,jj+2)   = &
                  lhs%mass(ii+1,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+1,jj+32)  = &
                  lhs%mass(ii+1,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))

             lhs%mass(ii+2,jj)     = &
                  lhs%mass(ii+2,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+2,jj+1)   = &
                  lhs%mass(ii+2,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+2,jj+2)   = &
                  lhs%mass(ii+2,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+2,jj+3)   = &
                  lhs%mass(ii+2,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))

             lhs%mass(ii+3,jj)     = &
                  lhs%mass(ii+3,jj)     + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+3,jj+1)   = &
                  lhs%mass(ii+3,jj+1)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+3,jj+2)   = &
                  lhs%mass(ii+3,jj+2)   + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))          
             lhs%mass(ii+3,jj+3)   = &
                  lhs%mass(ii+3,jj+3) + (integrator%getWeight(k)*jacobianDet(k)*2&
                  *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j)) 
          end do
       end do
    end do
  end subroutine calculateLHS

  subroutine calculateRHS(this, processInfo, rhs)
    implicit none
    class(CFDElementDT)                          , intent(inout) :: this
    type(ProcessInfoDT)                               , intent(inout) :: processInfo
    real(rkind)            , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                      :: i, j, nNode, nDof
    real(rkind)                                                         :: val1, val2, val3, val4
    real(rkind)            , dimension(:,:), allocatable                :: valuedSource
    real(rkind)            , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                               :: integrator
    nNode = this%getnNode()
    nDof = this%node(1)%getnDof()
    allocate(rhs(nNode*nDof))
    rhs = 0._rkind
    do i = 1, nNode
       if(this%node(i)%hasSource()) then
          val1 = this%node(i)%ptr%source(1)%evaluate(1, (/this%node(i)%getx(), this%node(i)%gety()/))
          val2 = this%node(i)%ptr%source(1)%evaluate(2, (/this%node(i)%getx(), this%node(i)%gety()/))
          val3 = this%node(i)%ptr%source(1)%evaluate(3, (/this%node(i)%getx(), this%node(i)%gety()/))
          val4 = this%node(i)%ptr%source(1)%evaluate(4, (/this%node(i)%getx(), this%node(i)%gety()/))
          rhs(nDof*i-3) = rhs(nDof*i-3) + val1
          rhs(nDof*i-2) = rhs(nDof*i-2) + val2
          rhs(nDof*i-1) = rhs(nDof*i-1) + val3
          rhs(nDof*i)   = rhs(nDof*i)   + val4
       end if
    end do
    if(this%hasSource()) then
       integrator = this%getIntegrator()
       allocate(valuedSource(4,integrator%getIntegTerms()))
       allocate(jacobianDet(integrator%getIntegTerms()))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       do i = 1, nNode
          val1 = 0._rkind
          val2 = 0._rkind
          val3 = 0._rkind
          val4 = 0._rkind
          do j = 1, integrator%getIntegTerms()
             val1 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(1,j)*jacobianDet(j)
             val2 = val2 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(2,j)*jacobianDet(j)
             val3 = val3 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(3,j)*jacobianDet(j)
             val4 = val4 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(4,j)*jacobianDet(j)
          end do
          rhs(i*nDof-3) = rhs(i*nDof-3) + val1
          rhs(i*nDof-2) = rhs(i*nDof-2) + val2
          rhs(i*nDof-1) = rhs(i*nDof-1) + val3
          rhs(i*nDof)   = rhs(i*nDof)   + val4
       end do
       deallocate(valuedSource)
       deallocate(jacobianDet)
    end if
  end subroutine calculateRHS

  subroutine setupIntegration(this, integrator, valuedSource, jacobianDet)
    implicit none
    class(CFDElementDT)                          , intent(inout) :: this
    type(IntegratorPtrDT)                               , intent(in)    :: integrator
    real(rkind), dimension(4,integrator%getIntegTerms()), intent(out)   :: valuedSource
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
    class(CFDElementDT), intent(inout) :: this
    type(IntegratorPtrDT) , intent(in) :: integrator
    real(rkind), dimension(4,integrator%getIntegTerms()) :: getValuedSource
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
       getValuedSource(3,i) = this%source(1)%evaluate(3, (/x,y/))
       getValuedSource(4,i) = this%source(1)%evaluate(4, (/x,y/))
    end do
  end function getValuedSource

  subroutine calculateResults(this, processInfo, resultMat)
    implicit none
    class(StructuralElementDT)                            , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
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
          nsx = nsx + dNidx*this%node(i)%ptr%dof(1)%val
          nsy = nsy + dNidy*this%node(i)%ptr%dof(2)%val
          shs = shs + dNidx*this%node(i)%ptr%dof(2)%val + dNidy*this%node(i)%ptr%dof(1)%val
          epx = epx + dNidx*this%node(i)%ptr%dof(1)%val
          epy = epy + dNidy*this%node(i)%ptr%dof(2)%val
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
  
  subroutine calculateSystem(this, processInfo, lhs)
    implicit none
    class(CFDElementDT)                                   , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)            , dimension(:,:,:), allocatable                :: resultMat
    integer(ikind)                                                        :: i, j, ii, jj, k
    integer(ikind)                                                        :: nNode, nDof
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(4,4)                               :: A1, A2, K11, K12, K21, K22
    real(rkind)            , dimension(:,:)  , allocatable                :: inverse
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    nDof = this%node(1)%getnDof()
    integrator = this%getIntegrator()
    allocate(nodalPoints(nNode))
    allocate(inverse(nNode*nDof,nNode*nDof))
    allocate(resultMat(nNode*nDof, nNode*nDof, 11))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    resultMat = 0.d0
    do i = 1, nNode
       do j = 1, nNode
          ii   = nDof*i-3
          jj   = nDof*j-3
          do k = 1, integrator%getIntegTerms()
             ! N * dN/dx
             resultMat(ii  ,jj  ,1) = resultMat(ii  ,jj  ,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
             resultMat(ii  ,jj+1,1) = resultMat(ii  ,jj+1,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii  ,jj+2,1) = resultMat(ii  ,jj+2,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii  ,jj+3,1) = resultMat(ii  ,jj+3,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

             resultMat(ii+1,jj  ,1) = resultMat(ii+1,jj  ,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
             resultMat(ii+1,jj+1,1) = resultMat(ii+1,jj+1,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii+1,jj+2,1) = resultMat(ii+1,jj+2,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii+1,jj+3,1) = resultMat(ii+1,jj+3,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

             resultMat(ii+2,jj  ,1) = resultMat(ii+2,jj  ,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii+2,jj+1,1) = resultMat(ii+2,jj+1,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii+2,jj+2,1) = resultMat(ii+2,jj+2,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
             resultMat(ii+2,jj+3,1) = resultMat(ii+2,jj+3,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

             resultMat(ii+3,jj  ,1) = resultMat(ii+3,jj  ,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii+3,jj+1,1) = resultMat(ii+3,jj+1,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
             resultMat(ii+3,jj+2,1) = resultMat(ii+3,jj+2,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
             resultMat(ii+3,jj+3,1) = resultMat(ii+3,jj+3,1) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

             !====================================================================
             ! N * dN/dy
             resultMat(ii  ,jj  ,2) = resultMat(ii  ,jj  ,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii  ,jj+1,2) = resultMat(ii  ,jj+1,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii  ,jj+2,2) = resultMat(ii  ,jj+2,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii  ,jj+3,2) = resultMat(ii  ,jj+3,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             resultMat(ii+1,jj  ,2) = resultMat(ii+1,jj  ,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
             resultMat(ii+1,jj+1,2) = resultMat(ii+1,jj+1,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+1,jj+2,2) = resultMat(ii+1,jj+2,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+1,jj+3,2) = resultMat(ii+1,jj+3,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             resultMat(ii+2,jj  ,2) = resultMat(ii+2,jj  ,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+2,jj+1,2) = resultMat(ii+2,jj+1,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+2,jj+2,2) = resultMat(ii+2,jj+2,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
             resultMat(ii+2,jj+3,2) = resultMat(ii+2,jj+3,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             resultMat(ii+3,jj  ,2) = resultMat(ii+3,jj  ,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+3,jj+1,2) = resultMat(ii+3,jj+1,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
             resultMat(ii+3,jj+2,2) = resultMat(ii+3,jj+2,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+3,jj+3,2) = resultMat(ii+3,jj+3,2) &
                  +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             !====================================================================
             ! dN/dx**2
             resultMat(ii  ,jj  ,3) = resultMat(ii  ,jj  ,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii  ,jj+1,3) = resultMat(ii  ,jj+1,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii  ,jj+2,3) = resultMat(ii  ,jj+2,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii  ,jj+3,3) = resultMat(ii  ,jj+3,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

             resultMat(ii+1,jj  ,3) = resultMat(ii+1,jj  ,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)          
             resultMat(ii+1,jj+1,3) = resultMat(ii+1,jj+1,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii+1,jj+2,3) = resultMat(ii+1,jj+2,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii+1,jj+3,3) = resultMat(ii+1,jj+3,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

             resultMat(ii+2,jj  ,3) = resultMat(ii+2,jj  ,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii+2,jj+1,3) = resultMat(ii+2,jj+1,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii+2,jj+2,3) = resultMat(ii+2,jj+2,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)          
             resultMat(ii+2,jj+3,3) = resultMat(ii+2,jj+3,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

             resultMat(ii+3,jj  ,3) = resultMat(ii+3,jj  ,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii+3,jj+1,3) = resultMat(ii+3,jj+1,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)          
             resultMat(ii+3,jj+2,3) = resultMat(ii+3,jj+2,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
             resultMat(ii+3,jj+3,3) = resultMat(ii+3,jj+3,3) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

             !====================================================================
             ! dN/dy**2
             resultMat(ii  ,jj  ,4) = resultMat(ii  ,jj  ,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)
             resultMat(ii  ,jj+1,4) = resultMat(ii  ,jj+1,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)
             resultMat(ii  ,jj+2,4) = resultMat(ii  ,jj+2,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)
             resultMat(ii  ,jj+3,4) = resultMat(ii  ,jj+3,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

             resultMat(ii+1,jj  ,4) = resultMat(ii+1,jj  ,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)          
             resultMat(ii+1,jj+1,4) = resultMat(ii+1,jj+1,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
             resultMat(ii+1,jj+2,4) = resultMat(ii+1,jj+2,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
             resultMat(ii+1,jj+3,4) = resultMat(ii+1,jj+3,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

             resultMat(ii+2,jj  ,4) = resultMat(ii+2,jj  ,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
             resultMat(ii+2,jj+1,4) = resultMat(ii+2,jj+1,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
             resultMat(ii+2,jj+2,4) = resultMat(ii+2,jj+2,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)          
             resultMat(ii+2,jj+3,4) = resultMat(ii+2,jj+3,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

             resultMat(ii+3,jj  ,4) = resultMat(ii+3,jj  ,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
             resultMat(ii+3,jj+1,4) = resultMat(ii+3,jj+1,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)          
             resultMat(ii+3,jj+2,4) = resultMat(ii+3,jj+2,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
             resultMat(ii+3,jj+3,4) = resultMat(ii+3,jj+3,4) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

             !====================================================================
             ! dN/dx * dN/dy
             resultMat(ii  ,jj  ,5) = resultMat(ii  ,jj  ,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii  ,jj+1,5) = resultMat(ii  ,jj+1,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii  ,jj+2,5) = resultMat(ii  ,jj+2,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii  ,jj+3,5) = resultMat(ii  ,jj+3,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             resultMat(ii+1,jj  ,5) = resultMat(ii+1,jj  ,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
             resultMat(ii+1,jj+1,5) = resultMat(ii+1,jj+1,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+1,jj+2,5) = resultMat(ii+1,jj+2,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+1,jj+3,5) = resultMat(ii+1,jj+3,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             resultMat(ii+2,jj  ,5) = resultMat(ii+2,jj  ,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+2,jj+1,5) = resultMat(ii+2,jj+1,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+2,jj+2,5) = resultMat(ii+2,jj+2,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
             resultMat(ii+2,jj+3,5) = resultMat(ii+2,jj+3,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

             resultMat(ii+3,jj  ,5) = resultMat(ii+3,jj  ,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+3,jj+1,5) = resultMat(ii+3,jj+1,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
             resultMat(ii+3,jj+2,5) = resultMat(ii+3,jj+2,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
             resultMat(ii+3,jj+3,5) = resultMat(ii+3,jj+3,5) &
                  +(integrator%getWeight(k)&
                  *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))
             *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))
          end do
       end do
    end do
    gamma = this%material%gamma
    mu    = this%material%mu
    R     = this%material%R
    k     = this%material%k
    do i = 1, nNode
       v(1) = (this%node(i)%ptr%dof(2)%val/this%node(i)%ptr%dof(1)%val)
       v(2) = (this%node(i)%ptr%dof(3)%val/this%node(i)%ptr%dof(1)%val)
       E    = (this%node(i)%ptr%dof(4)%val/this%node(i)%ptr%dof(1)%val)
       A1  = MA1(v,gamma,E)
       A2  = MA2(v,gamma,E)
       K11 = MK11(v,mu,R,gamma,E,k)
       K12 = MK12(v,mu,R,gamma,E,k)
       K21 = MK21(v,mu,R,gamma,E,k)
       K22 = MK22(v,mu,R,gamma,E,k)
       do j = 1, nNode
          ii   = nDof*i-3
          jj   = nDof*j-3
          resultMat(ii:ii+3,jj:jj+3,6 ) = resultMat(ii:ii+3,jj:jj+3,6 ) +  A1(:,:)
          resultMat(ii:ii+3,jj:jj+3,7 ) = resultMat(ii:ii+3,jj:jj+3,7 ) +  A2(:,:)
          resultMat(ii:ii+3,jj:jj+3,8 ) = resultMat(ii:ii+3,jj:jj+3,8 ) + K11(:,:)
          resultMat(ii:ii+3,jj:jj+3,9 ) = resultMat(ii:ii+3,jj:jj+3,9 ) + K12(:,:)
          resultMat(ii:ii+3,jj:jj+3,10) = resultMat(ii:ii+3,jj:jj+3,10) + K21(:,:)
          resultMat(ii:ii+3,jj:jj+3,11) = resultMat(ii:ii+3,jj:jj+3,11) + K22(:,:)
          
       end do
    end do
    inverse = 0._rkind
    do i = 1, nNode*nDof
       do j = 1, nNode*nDof
          inverse(i,i) = inverse(i,i) + inverse(i,j)
       end do
       inverse(i,i) = 1._rkind/inverse(i,i)
    end do
    tau  = calculateTau()
    mu_a = calculateMu_a()
    lhs%stiffness =                      matmul(resultMat(:,:,1),resultMat(:,:,6)) &
         +                              matmul(resultMat(:,:,2),resultMat(:,:,7))  &
         + tau*(matmul(resultMat(:,:,3),matmul(resultMat(:,:,6),resultMat(:,:,6))) &
         +      matmul(resultMat(:,:,5),matmul(resultMat(:,:,6),resultMat(:,:,7))) &
         +      matmul(resultMat(:,:,5),matmul(resultMat(:,:,7),resultMat(:,:,6))) &
         +      matmul(resultMat(:,:,4),matmul(resultMat(:,:,7),resultMat(:,:,7))) &
         -      matmul(inverse         ,(matmul(resultMat(:,:,1),resultMat(:,:,6)) &
         +                              matmul(resultMat(:,:,2),resultMat(:,:,7))))&
         + mu_a*(resultMat(:,:,3)+resultMat(:,:,4)+2._rkind*resultMat(:,:,5))      &
         +                              matmul(resultMat(:,:,3),resultMat(:,:,8))  &
         +                              matmul(resultMat(:,:,5),resultMat(:,:,9))  &
         +                              matmul(resultMat(:,:,5),resultMat(:,:,10)) &
         +                              matmul(resultMat(:,:,4),resultMat(:,:,11))          
  end subroutine calculateResults

  subroutine calculateDT(this, processInfo, dt)
    implicit none
    class(CFDElementDT)          , intent(inout) :: this
    type(ProcessInfoDT)          , intent(inout) :: processInfo
    real(rkind), dimension(:,:,:), allocatable   :: jacobian
    real(rkind), dimension(:)    , allocatable   :: jacobianDet
    type(IntegratorPtrDT)                        :: integrator
    real(rkind)                  , intent(inout) :: dt
    integer(ikind)                               :: i, nNode, nDof
    real(rkind)                                  :: area, dt_min, alpha, deltaTU   
    real(rkind)                                  :: val1, val2, V, dt_elem, deltaTC    
    real(rkind)                                  :: Vx, Vy, Vxmax, Vymax, fSafe, cota 
    type(NodePtrDT), dimension(:), allocatable   :: nodalPoints
    fSafe = processInfo%getConstants(1)
    dt_min = 1.d20
    Vxmax = 0.d0
    Vymax = 0.d0
    nNode = this%getnNode()
    nDof = this%node(1)%getnDof()
    integrator = this%getIntegrator()
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    area = this%geometry%getLenght(nodalPoints)
    do i = 1, nNode
       val1 = 0._rkind
       val2 = 0._rkind
       do j = 1, integrator%getIntegTerms()
          val1 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *this%node(i)%ptr%dof(2)%val/this%node(i)%ptr%dof(1)%val*jacobianDet(j)
          val2 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *this%node(i)%ptr%dof(3)%val/this%node(i)%ptr%dof(1)%val*jacobianDet(j)
       end do
       Vx = Vx + val1
       Vy = Vy + val2
    end do
    if (Vx .gt. Vxmax) Vxmax = Vx
    if (Vy .gt. Vymax) Vymax = Vy
    V       = sqrt((Vxmax**2+Vymax**2))
    alpha   = min(V*sqrt(jacobianDet(1))/(2._rkind*0.001d0)/3._rkind,1._rkind)
    deltaTU = 1._rkind/(4._rkind*0.001d0/sqrt(jacobianDet(1))**2._rkind+alpha*V/sqrt(jacobianDet(1)))
    deltaTC = 1._rkind/(4._rkind*0.001d0/sqrt(jacobianDet(1))**2._rkind)
    dt_elem  = fsafe/(1._rkind/deltaTC+1._rkind/deltaTU)
    dt = dt_elem
    IF (dt_elem .LT. dt_min) dt_min = dt_elem
    cota = 10.d0*dt_min
    !EL VALOR 10 ESTA PUESTO A OJO #MODIFICAR SI ES NECESARIO#
    if(dt > cota) dt = cota
  end subroutine calculateDT

  function calculateTau()
    implicit none
  end function calculateTau

  function calculateMu_a()
    implicit none
  end function calculateMu_a
  
  function MA1(v,gamma,E)
    implicit none
    real(rkind), dimension(4,4) :: A1
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: gamma
    real(rkind)                 :: E
    A1(1,1) = 0._rkind
    A1(1,2) = 1._rkind
    A1(1,3) = 0._rkind
    A1(1,4) = 0._rkind
    
    A1(2,1) = (((gamma-1._rkind)/2._rkind)*(v(1)**2+v(2)**2))-v(1)**2
    A1(2,2) = (3._rkind-gamma)*v(1)
    A1(2,3) = -(gamma-1._rkind)*v(2)
    A1(2,4) = (gamma-1._rkind)
    
    A1(3,1) = -v(1)*v(2)
    A1(3,2) = v(2)
    A1(3,3) = v(1)
    A1(3,4) = 0._rkind
    
    A1(4,1) = (((gamma-1._rkind)*(v(1)**2+v(2)**2))-(gamma*E))*v(1)
    A1(4,2) = (gamma*E)-((gamma-1._rkind)*(v(1)**2+v(2)**2)/2._rkind)&
         -((gamma-1._rkind)*v(1)*v(2))
    A1(4,3) = -(gamma-1._rkind)*v(1)*v(2)
    A1(4,4) = gamma*v(1)
  end function MA1

  function MA2(v,gamma,E)
    implicit none
    real(rkind), dimension(4,4) :: A2
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: gamma
    real(rkind)                 :: E
    A2(1,1) = 0._rkind
    A2(1,2) = 0._rkind
    A2(1,3) = 1._rkind
    A2(1,4) = 0._rkind
    
    A2(2,1) = -v(1)*v(2)
    A2(2,2) = v(2)
    A2(2,3) = v(1)
    A2(2,4) = 0._rkind
    
    A2(3,1) = (((gamma-1._rkind)/2._rkind)*(v(1)**2+v(2)**2))-v(2)**2
    A2(3,2) = -(gamma-1._rkind)*v(1)
    A2(3,3) = (3._rkind-gamma)*v(2)
    A2(3,4) = (gamma-1._rkind)
    
    A2(4,1) = (((gamma-1._rkind)*(v(1)**2+v(2)**2))-(gamma*E))*v(2)
    A2(4,2) = -(gamma-1._rkind)*v(1)*v(2)
    A2(4,3) = (gamma*E)-((gamma-1._rkind)*(v(1)**2+v(2)**2)/2._rkind)&
         -((gamma-1._rkind)*v(2)**2)
    A2(4,4) = gamma*v(2)
  end function MA2

  function MK11(v,mu,R,gamma,E,k)
    implicit none
    real(rkind), dimension(4,4) :: K11
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: mu
    real(rkind)                 :: Cv
    real(rkind)                 :: E
    real(rkind)                 :: k
    real(rkind)                 :: R
    real(rkind)                 :: gamma 
    Cv = R/(gamma-1)
    K11(1,1) = 0._rkind
    K11(1,2) = 0._rkind
    K11(1,3) = 0._rkind 
    K11(1,4) = 0._rkind
    
    K11(2,1) = -(4._rkind/3._rkind)*mu*v(1)
    K11(2,2) = (4._rkind/3._rkind)*mu
    K11(2,3) = 0._rkind
    K11(2,4) = 0._rkind
    
    K11(3,1) = -mu*v(2)
    K11(3,2) = 0._rkind
    K11(3,3) = mu
    K11(3,4) = 0._rkind
    
    K11(4,1) = ((k/Cv)*(((v(1)**2+v(2)**2)/2._rkind)    &
         -(E-((v(1)**2+v(2)**2)/2._rkind))))              &
         -(mu*(v(1)**2+v(2)**2))-(mu*v(1)*v(1)/3._rkind)
    K11(4,2) = ((mu/3._rkind)+mu-(k/Cv))*v(1)
    K11(4,3) = (mu-(k/Cv))*v(2)
    K11(4,4) = (k/Cv)
  end function MK11

  function MK12(v,mu,R,gamma,E,k)
    implicit none
    real(rkind), dimension(4,4) :: K12
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: mu
    real(rkind)                 :: Cv
    real(rkind)                 :: E
    real(rkind)                 :: k
    real(rkind)                 :: R
    real(rkind)                 :: gamma 
    Cv = R/(gamma-1)
    K12(1,1) = 0._rkind
    K12(1,2) = 0._rkind
    K12(1,3) = 0._rkind
    K12(1,4) = 0._rkind
    
    K12(2,1) = 2._rkind*mu*v(2)/3._rkind
    K12(2,2) = 0._rkind
    K12(2,3) = -2._rkind*mu/3._rkind
    K12(2,4) = 0._rkind
    
    K12(3,1) = -mu*v(1)
    K12(3,2) = mu
    K12(3,3) = 0._rkind
    K12(3,4) = 0._rkind
    
    K12(4,1) = -mu*v(1)*v(2)/3._rkind
    K12(4,2) = mu*v(2)
    K12(4,3) = -2._rkind*mu*v(1)/3._rkind
    K12(4,4) = 0._rkind
  end function MK12

  function MK21(v,mu,R,gamma,E,k)
    implicit none
    real(rkind), dimension(4,4) :: K21
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: mu
    real(rkind)                 :: Cv
    real(rkind)                 :: E
    real(rkind)                 :: k
    real(rkind)                 :: R
    real(rkind)                 :: gamma 
    Cv = R/(gamma-1)
    K21(1,1) = 0._rkind
    K21(1,2) = 0._rkind
    K21(1,3) = 0._rkind
    K21(1,4) = 0._rkind
    
    K21(2,1) = -mu*v(2)
    K21(2,2) = 0._rkind
    K21(2,3) = mu
    K21(2,4) = 0._rkind
    
    K21(3,1) = (2._rkind*mu*v(1)/3._rkind)-(2._rkind*mu/3._rkind)
    K21(3,2) = 0._rkind
    K21(3,3) = 0._rkind
    K21(3,4) = 0._rkind
    
    K21(4,1) = -mu*v(1)*v(2)/3._rkind
    K21(4,2) = -2._rkind*mu*v(2)/3._rkind
    K21(4,3) = mu*v(1)
    K21(4,4) = 0._rkind
  end function MK21

  function MK22(v,mu,R,gamma,E,k)
    implicit none
    real(rkind), dimension(4,4) :: K22
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: mu
    real(rkind)                 :: Cv
    real(rkind)                 :: E
    real(rkind)                 :: k
    real(rkind)                 :: R
    real(rkind)                 :: gamma 
    Cv = R/(gamma-1)
    K22(1,1) = 0._rkind
    K22(1,2) = 0._rkind
    K22(1,3) = 0._rkind
    K22(1,4) = 0._rkind
    
    K22(2,1) = -mu*v(1)
    K22(2,2) = mu
    K22(2,3) = 0._rkind
    K22(2,4) = 0._rkind
    
    K22(3,1) = -(4._rkind/3._rkind)*mu*v(2)
    K22(3,2) = 0._rkind
    K22(3,3) = (4._rkind/3._rkind)*mu
    K22(3,4) = 0._rkind
    
    K22(4,1) = ((k/Cv)*(((v(1)**2+v(2)**2)/2._rkind)    &
         -(E-((v(1)**2+v(2)**2)/2._rkind))))              &
         -(mu*(v(1)**2+v(2)**2))-(mu*v(2)*v(2)/3._rkind)
    K22(4,2) = (mu-(k/Cv))*v(1)
    K22(4,3) = ((mu/3._rkind)+mu-(k/Cv))*v(2)
    K22(4,4) = (k/Cv)
  end function MK22
  
end module CFDElementM

!!$    integer(ikind)                    :: i
!!$    integer(ikind)                    :: n
!!$    n = model%getnNode()
!!$    allocate(model%results%density(n)       , model%results%velocity(n,n) )
!!$    allocate(model%results%internalEnergy(n), model%results%temperature(n))
!!$    allocate(model%results%pressure(n)      , model%results%mach(n)       )
!!$    do i = 1, size(model%dof),4
!!$       model%results%density(i)        = model%dof(i)
!!$       model%results%velocity(i,i)     = (/model%dof(i+1),model%dof(i+2)/)/model%dof(i)
!!$       model%results%internalEnergy(i) = model%dof(i+3)/model%dof(i)
!!$       model%results%temperature(i)    = &
!!$            (model%gamma-1)/R*(model%dof(i+3)-0.5*(model%dof(i+1)**2+model%dof(i+2)**2))
!!$       model%results%pressure(i)       = &
!!$            (model%gamma-1)*model%dof(i)*(model%dof(i+3)-0.5*(model%dof(i+1)**2+model%dof(i+2)**2))
!!$       model%results%mach(i)           = sqrt(model%dof(i+1)**2+model%dof(i+2)**2)         &
!!$            /sqrt(model%gamma*model%R                                                   &
!!$            *(model%gamma-1)/R*(model%dof(i+3)-0.5*(model%dof(i+1)**2+model%dof(i+2)**2)))
!!$    end do
