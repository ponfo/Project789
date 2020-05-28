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
     class(CFDMaterialDT), pointer :: material
   contains
     procedure, public  :: init
     procedure, public  :: calculateLHS
     procedure, public  :: calculateRHS
     procedure, public  :: calculateLocalSystem
     procedure, public  :: calculateResults
     procedure, public  :: calculateDt
     procedure, public  :: calculateTauNu
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
  real(rkind), dimension(:,:,:), allocatable, save :: localResultMat

contains

  type(CFDElementDT) function constructor(id, node, material)
    implicit none
    integer(ikind)                          , intent(in) :: id
    type(NodePtrDT)           , dimension(:), intent(in) :: node
    class(CFDMaterialDT), target            , intent(in) :: material
    call constructor%init(id, node, material)
  end function constructor

  subroutine init(this, id, node, material)
    implicit none
    class(CFDElementDT)                     , intent(inout) :: this
    integer(ikind)                          , intent(in)    :: id
    type(NodePtrDT)           , dimension(:), intent(in)    :: node
    class(CFDMaterialDT), target            , intent(in)    :: material
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
             lhs%mass(ii+1,jj+2)  = &
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
    call this%calculateDt(processInfo)
    deallocate(jacobian)
    deallocate(jacobianDet)
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
    class(CFDElementDT)                                   , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    real(rkind)            , dimension(:,:,:), allocatable, intent(inout) :: resultMat
    integer(ikind)                                                        :: i, j, ii, jj, k
    integer(ikind)                                                        :: nNode, nDof
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:,:)  , allocatable                :: auxMatrix
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    if (processInfo%getStep()==0) then
       nNode = this%getnNode()
       nDof = this%node(1)%getnDof()
       integrator = this%getIntegrator()
       allocate(nodalPoints(nNode))
       if (allocated(localResultMat)) deallocate(localResultMat)
       allocate(localResultMat(nNode*nDof, nNode*nDof, 6))
       do i = 1, nNode
          nodalPoints(i) = this%node(i)
       end do
       jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
       jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
       localResultMat = 0.d0
       do i = 1, nNode
          do j = 1, nNode
             ii   = nDof*i-3
             jj   = nDof*j-3
             do k = 1, integrator%getIntegTerms()
                ! N * dN/dx
                localResultMat(ii  ,jj  ,1) = localResultMat(ii  ,jj  ,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
                localResultMat(ii  ,jj+1,1) = localResultMat(ii  ,jj+1,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii  ,jj+2,1) = localResultMat(ii  ,jj+2,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii  ,jj+3,1) = localResultMat(ii  ,jj+3,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

                localResultMat(ii+1,jj  ,1) = localResultMat(ii+1,jj  ,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
                localResultMat(ii+1,jj+1,1) = localResultMat(ii+1,jj+1,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii+1,jj+2,1) = localResultMat(ii+1,jj+2,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii+1,jj+3,1) = localResultMat(ii+1,jj+3,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

                localResultMat(ii+2,jj  ,1) = localResultMat(ii+2,jj  ,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii+2,jj+1,1) = localResultMat(ii+2,jj+1,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii+2,jj+2,1) = localResultMat(ii+2,jj+2,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
                localResultMat(ii+2,jj+3,1) = localResultMat(ii+2,jj+3,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

                localResultMat(ii+3,jj  ,1) = localResultMat(ii+3,jj  ,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii+3,jj+1,1) = localResultMat(ii+3,jj+1,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))          
                localResultMat(ii+3,jj+2,1) = localResultMat(ii+3,jj+2,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k))           
                localResultMat(ii+3,jj+3,1) = localResultMat(ii+3,jj+3,1) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))/jacobianDet(k)) 

                !====================================================================
                ! N * dN/dy
                localResultMat(ii  ,jj  ,2) = localResultMat(ii  ,jj  ,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii  ,jj+1,2) = localResultMat(ii  ,jj+1,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii  ,jj+2,2) = localResultMat(ii  ,jj+2,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii  ,jj+3,2) = localResultMat(ii  ,jj+3,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                localResultMat(ii+1,jj  ,2) = localResultMat(ii+1,jj  ,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
                localResultMat(ii+1,jj+1,2) = localResultMat(ii+1,jj+1,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+1,jj+2,2) = localResultMat(ii+1,jj+2,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+1,jj+3,2) = localResultMat(ii+1,jj+3,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                localResultMat(ii+2,jj  ,2) = localResultMat(ii+2,jj  ,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+2,jj+1,2) = localResultMat(ii+2,jj+1,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+2,jj+2,2) = localResultMat(ii+2,jj+2,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
                localResultMat(ii+2,jj+3,2) = localResultMat(ii+2,jj+3,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                localResultMat(ii+3,jj  ,2) = localResultMat(ii+3,jj  ,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+3,jj+1,2) = localResultMat(ii+3,jj+1,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
                localResultMat(ii+3,jj+2,2) = localResultMat(ii+3,jj+2,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+3,jj+3,2) = localResultMat(ii+3,jj+3,2) &
                     +(integrator%getWeight(k)*integrator%getShapeFunc(k,i)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                !====================================================================
                ! dN/dx**2
                localResultMat(ii  ,jj  ,3) = localResultMat(ii  ,jj  ,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii  ,jj+1,3) = localResultMat(ii  ,jj+1,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii  ,jj+2,3) = localResultMat(ii  ,jj+2,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii  ,jj+3,3) = localResultMat(ii  ,jj+3,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

                localResultMat(ii+1,jj  ,3) = localResultMat(ii+1,jj  ,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)          
                localResultMat(ii+1,jj+1,3) = localResultMat(ii+1,jj+1,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii+1,jj+2,3) = localResultMat(ii+1,jj+2,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii+1,jj+3,3) = localResultMat(ii+1,jj+3,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

                localResultMat(ii+2,jj  ,3) = localResultMat(ii+2,jj  ,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii+2,jj+1,3) = localResultMat(ii+2,jj+1,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii+2,jj+2,3) = localResultMat(ii+2,jj+2,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)          
                localResultMat(ii+2,jj+3,3) = localResultMat(ii+2,jj+3,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

                localResultMat(ii+3,jj  ,3) = localResultMat(ii+3,jj  ,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii+3,jj+1,3) = localResultMat(ii+3,jj+1,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)          
                localResultMat(ii+3,jj+2,3) = localResultMat(ii+3,jj+2,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2)           
                localResultMat(ii+3,jj+3,3) = localResultMat(ii+3,jj+3,3) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))**2/jacobianDet(k)**2) 

                !====================================================================
                ! dN/dy**2
                localResultMat(ii  ,jj  ,4) = localResultMat(ii  ,jj  ,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)
                localResultMat(ii  ,jj+1,4) = localResultMat(ii  ,jj+1,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)
                localResultMat(ii  ,jj+2,4) = localResultMat(ii  ,jj+2,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)
                localResultMat(ii  ,jj+3,4) = localResultMat(ii  ,jj+3,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

                localResultMat(ii+1,jj  ,4) = localResultMat(ii+1,jj  ,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)          
                localResultMat(ii+1,jj+1,4) = localResultMat(ii+1,jj+1,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
                localResultMat(ii+1,jj+2,4) = localResultMat(ii+1,jj+2,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
                localResultMat(ii+1,jj+3,4) = localResultMat(ii+1,jj+3,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

                localResultMat(ii+2,jj  ,4) = localResultMat(ii+2,jj  ,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
                localResultMat(ii+2,jj+1,4) = localResultMat(ii+2,jj+1,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
                localResultMat(ii+2,jj+2,4) = localResultMat(ii+2,jj+2,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)          
                localResultMat(ii+2,jj+3,4) = localResultMat(ii+2,jj+3,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

                localResultMat(ii+3,jj  ,4) = localResultMat(ii+3,jj  ,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
                localResultMat(ii+3,jj+1,4) = localResultMat(ii+3,jj+1,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)          
                localResultMat(ii+3,jj+2,4) = localResultMat(ii+3,jj+2,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2)           
                localResultMat(ii+3,jj+3,4) = localResultMat(ii+3,jj+3,4) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))**2/jacobianDet(k)**2) 

                !====================================================================
                ! dN/dx * dN/dy
                localResultMat(ii  ,jj  ,5) = localResultMat(ii  ,jj  ,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii  ,jj+1,5) = localResultMat(ii  ,jj+1,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii  ,jj+2,5) = localResultMat(ii  ,jj+2,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii  ,jj+3,5) = localResultMat(ii  ,jj+3,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                localResultMat(ii+1,jj  ,5) = localResultMat(ii+1,jj  ,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
                localResultMat(ii+1,jj+1,5) = localResultMat(ii+1,jj+1,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+1,jj+2,5) = localResultMat(ii+1,jj+2,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+1,jj+3,5) = localResultMat(ii+1,jj+3,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                localResultMat(ii+2,jj  ,5) = localResultMat(ii+2,jj  ,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+2,jj+1,5) = localResultMat(ii+2,jj+1,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+2,jj+2,5) = localResultMat(ii+2,jj+2,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
                localResultMat(ii+2,jj+3,5) = localResultMat(ii+2,jj+3,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k)) 

                localResultMat(ii+3,jj  ,5) = localResultMat(ii+3,jj  ,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+3,jj+1,5) = localResultMat(ii+3,jj+1,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))          
                localResultMat(ii+3,jj+2,5) = localResultMat(ii+3,jj+2,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))           
                localResultMat(ii+3,jj+3,5) = localResultMat(ii+3,jj+3,5) &
                     +(integrator%getWeight(k)&
                     *(jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                     - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i))&
                     *(jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                     - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j))/jacobianDet(k))
             end do
          end do
       end do
    end if
    call this%calculateTauNu(processInfo, auxMatrix)
    localResultMat(1:5,1:5,6) = auxMatrix(:,:)
    resultMat = localResultMat
    call this%calculateDt(processInfo)
  end subroutine calculateResults

  subroutine calculateDt(this, processInfo)
    implicit none
    class(CFDElementDT), intent(inout)         :: this
    class(ProcessInfoDT), intent(inout)        :: processInfo
    integer(ikind)                             :: i, j, nNode, nodeID, nDof, nCond, iCond
    integer(ikind)                             :: iNode, iNodeID , nElem, iElem
    type(NodePtrDT)                            :: node
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:)    , allocatable :: jacobianDet, constants, dof
    type(IntegratorPtrDT)                      :: integrator
    real(rkind)                                :: dt_min, alpha, deltaTU, dt, rho, val3   
    real(rkind)                                :: V, dt_elem, deltaTC    
    real(rkind)                                :: Vx, Vy, Vxmax, Vymax, fSafe, cota 
    type(NodePtrDT), dimension(:), allocatable :: nodalPoints
    constants = processInfo%getConstants()
    nNode = this%getnNode()
    fSafe = constants(1)
    dof = constants(13:size(constants))
    dt_min = 1.d20
    Vxmax = 0.d0
    Vymax = 0.d0
    integrator = this%getIntegrator()
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%getNode(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    do i = 1, nNode
       Vx  = 0.d0
       Vy  = 0.d0
       rho = 0.d0
       nodeID = this%getNodeID(i)
       do j = 1, integrator%getIntegTerms()
          Vx = rho + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4-2)*jacobianDet(j)
          Vy = rho + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4-1)*jacobianDet(j)
          rho = rho + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4-3)*jacobianDet(j)
       end do
       if (rho .ne. 0.d0) then
          Vx  = Vx/rho
          Vy  = Vy/rho
       end if
       if (Vx .gt. Vxmax) Vxmax = Vx
       if (Vy .gt. Vymax) Vymax = Vy
    end do
    V       = sqrt((Vxmax**2+Vymax**2))*0.5d0
    alpha   = min(V*sqrt(jacobianDet(1))/(2._rkind*0.001d0)/3._rkind,1._rkind)
    deltaTU = 1._rkind/(4._rkind*0.001d0/sqrt(jacobianDet(1))**2._rkind+alpha*V/sqrt(jacobianDet(1)))
    deltaTC = 1._rkind/(4._rkind*0.001d0/sqrt(jacobianDet(1))**2._rkind)
    dt_elem  = fsafe/(1._rkind/deltaTC+1._rkind/deltaTU)
    dt = dt_elem
    IF (dt_elem .LT. dt_min) dt_min = dt_elem
    cota = dt_min
    if(dt > cota) dt = cota
    call processInfo%setMinimumDt(dt)
  end subroutine calculateDt

  subroutine calculateTauNu(this, processInfo, matrix)
    implicit none
    class(CFDElementDT), intent(inout)         :: this
    class(ProcessInfoDT), intent(inout)        :: processInfo
    real(rkind), dimension(:,:), allocatable, intent(inout) :: matrix
    integer(ikind)                             :: i, j, nNode, nodeID, nDof, nCond, iCond
    integer(ikind)                             :: iNode, iNodeID , nElem, iElem
    type(NodePtrDT)                            :: node
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:)    , allocatable :: jacobianDet, constants, dof
    type(IntegratorPtrDT)                      :: integrator
    type(NodePtrDT), dimension(:), allocatable :: nodalPoints
    real(rkind)                                :: rho, H_JGN, H_RGN, H_RGNE, val3, val4, nu  
    real(rkind)                                :: val1, val2, V, deltaTC, dRho, dV, dT, dtime
    real(rkind)                                :: dxRho, dxV, dxT, dyRho, dyV, dyT, T, Tinf, dNidy
    real(rkind)                                :: Vx, Vy, Vxmax, Vymax, gamma, cota, Cv, Vc, dNidx
    real(rkind)                                :: VxT, VyT, VxJ, VyJ, VxV, VyV, fmu, smu, bi, ci 
    real(rkind)                                :: cte, Rhoinf, Tau, T_SUGN1, T_SUGN2, T_SUGN3
    real(rkind)                                :: term1, term2, TR1, Tau_SUNG3, Tau_SUNG3_E
    real(rkind)                                :: H_RGN1, H_RGN2, resumen
    constants = processInfo%getConstants()
    nNode = this%getnNode()
    if (allocated(nodalPoints)) deallocate(nodalPoints)
    allocate(nodalPoints(nNode))
    if (allocated(matrix)) deallocate(matrix)
    allocate(matrix(5,5))
    matrix = 0.d0
    cte = constants(2)
    Cv = constants(4)
    Vc = constants(5)
    dtime = processInfo%getDt()
    Tinf = this%material%T_inf
    Rhoinf = this%material%Rho
    dof = constants(13:size(constants))
    gamma = constants(6)
    Tau = 0.d0
    H_RGNE = 1.D-10
    H_RGN  = 1.D-10
    H_JGN  = 1.D-10
    integrator = this%getIntegrator()
    do i = 1, nNode
       nodalPoints(i) = this%getNode(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    rho = 0.d0
    Vx  = 0.d0
    Vy  = 0.d0
    T   = 0.d0
    do i = 1, nNode
       val1 = 0._rkind
       val2 = 0._rkind
       val3 = 0._rkind
       val4 = 0._rkind
       nodeID = this%getNodeID(i)
       do j = 1, integrator%getIntegTerms()
          val1 = val1 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4-3)*jacobianDet(j)
          val2 = val2 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4-2)*jacobianDet(j)
          val3 = val3 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4-1)*jacobianDet(j)
          val4 = val4 + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
               *dof(nodeID*4)*jacobianDet(j)
       end do
       rho = rho + val1
       Vx  = Vx  + val2
       Vy  = Vy  + val3
       T   = T   + val4
    end do
    if (rho .ne. 0.d0) then
       rho = rho/nNode
       Vx  = Vx/(rho*nNode)
       Vy  = Vy/(rho*nNode)
       T   = ((T/(rho*nNode))-0.5d0*(Vx**2+Vy**2))/Cv
    end if
    V = sqrt((Vx**2+Vy**2))
    dxRho = 0._rkind
    dxT   = 0._rkind
    dxV   = 0._rkind
    dyRho = 0._rkind
    dyT   = 0._rkind
    dyV   = 0._rkind
    do i = 1, nNode
       nodeID = this%getNodeID(i)
       do j = 1, integrator%getIntegTerms()
          bi = jacobian(j,2,2)*integrator%getDShapeFunc(j,1,i) &
               - jacobian(j,1,2)*integrator%getDShapeFunc(j,2,i)
          ci = jacobian(j,1,1)*integrator%getDShapeFunc(j,2,i) &
               - jacobian(j,2,1)*integrator%getDShapeFunc(j,1,i)
          dNidx = bi/jacobianDet(j)
          dNidy = ci/jacobianDet(j)
          dxRho = dxRho + dNidx*rho
          dxT   = dxT   + dNidy*T
          dxV   = dxV   + dNidx*V
          dyRho = dyRho + dNidy*rho
          dyT   = dyT   + dNidx*T
          dyV   = dyV   + dNidy*V
       end do
    end do
    dRho =sqrt(dxRho**2+dyRho**2)+1.d-20
    dT   =sqrt(dxT**2+dyT**2)+1.d-20
    dV   = sqrt(dxV**2+dyV**2)+1.d-20
    if (dT .ne. 0.d0) then
       VxT = dxT/dT
       VyT = dyT/dT
    else
       VxT = dxT
       VyT = dyT
    end if
    if (dRho .ne. 0.d0) then
       VxJ = dxRho/dRho
       VyJ = dyRho/dRho
    else
       VxJ = dxRho
       VyJ = dyRho
    end if
    if (dV .ne. 0.d0) then
       VxV = dxV/dV
       VyV = dyV/dV
    else
       VxV = dxV
       VyV = dyV
    end if
    smu = 110.d0
    fmu = 0.41685d0*(abs(T)/Tinf)**1.5d0*(Tinf+smu)/(abs(T)+smu)
    term1  = 0._rkind
    term2  = 0._rkind
    H_RGN1 = 0._rkind
    H_RGN2 = 0._rkind
    do j = 1, nNode
       nodeID = this%getNodeID(j)
       do i = 1, integrator%getIntegTerms()
          bi = jacobian(i,2,2)*integrator%getDShapeFunc(i,1,j) &
               - jacobian(i,1,2)*integrator%getDShapeFunc(i,2,j)
          ci = jacobian(i,1,1)*integrator%getDShapeFunc(i,2,j) &
               - jacobian(i,2,1)*integrator%getDShapeFunc(i,1,j)
          dNidx  = bi/jacobianDet(i)
          dNidy  = ci/jacobianDet(i)
          term1  = abs(dNidx*Vx  + dNidy*Vy)
          term2  = abs(dNidx*VxJ + dNidy*VyJ)
          H_RGN1 = abs(dNidx*VxT + dNidy*VyT)
          H_RGN2 = abs(dNidx*VxV + dNidy*VyV)
       end do
       Tau    = Tau+term1+term2*Vc
       H_RGNE = H_RGNE+H_RGN1
       H_RGN  = H_RGN+H_RGN2
       H_JGN  = H_JGN+term2
    end do
    if (Tau .ne. 0.d0) then
       Tau = 1.d0/Tau
    end if
    if (H_RGNE .ne. 0.d0) then
       H_RGNE = 2.d0/H_RGNE
    end if
    if (H_RGN .ne. 0.d0) then
       H_RGN  = 2.d0/H_RGN
    end if
    if (H_RGN .gt. 1.d3) H_RGN=0.d0
    if (H_JGN .ne. 0.d0) then
       H_JGN=2.D0/H_JGN
    end if
    if (H_JGN .gt. 1.d10) H_JGN=0.d0
    if (Rho == 0.d0) then
       TR1 = 0.d0
    else
       TR1 = dRho*H_JGN/Rho
    end if
    nu  = (sqrt(abs(TR1))+TR1**2)*Vc*2*(H_JGN/2.d0)
    if (Tau .ne. 0.d0) then
       resumen = ((1.d0/Tau)**2.d0 +(2.d0/dtime)**2.d0)**(-0.5d0)
    end if
    T_SUGN1 = Resumen
    T_SUGN2 = Resumen
    T_SUGN3 = Resumen
    if (fmu.ne.0.d0) then
       TAU_SUNG3   = H_RGN**2.d0/(4.d0*fmu/Rhoinf)
       TAU_SUNG3_E = H_RGNE**2.d0/(4.d0*fmu/Rhoinf)
       T_SUGN2     = (((1.d0/Tau)**2.d0 +(2.d0/dtime)**2.d0)&
            +1.d0/TAU_SUNG3**2.d0)**(-0.5d0)
       T_SUGN3     = (((1.d0/Tau)**2.d0 +(2.d0/dtime)**2.d0)&
            +1.d0/TAU_SUNG3_E**2.d0)**(-0.5d0)
    end if
    matrix(1,1) = T_SUGN1
    matrix(2,2) = T_SUGN2
    matrix(3,3) = T_SUGN2
    matrix(4,4) = T_SUGN3
    matrix(5,1) = nu*cte
  end subroutine calculateTauNu
  
end module CFDElementM
