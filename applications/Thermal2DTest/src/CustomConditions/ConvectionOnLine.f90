module ConvectionOnLineM
  use UtilitiesM

  use PointM
  use NodeM
  use NodePtrM
  use GeometryM

  use IntegratorPtrM
  
  use ConditionM

  implicit none

  private
  public :: ConvectionOnLineDT, convectionOnLine

  type, extends(ConditionDT) :: ConvectionOnLineDT
     integer(ikind), dimension(:), allocatable :: nodeIDList
     real(rkind)                               :: coef
     real(rkind)                               :: temp
   contains
     procedure, public :: init

     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type ConvectionOnLineDT

  interface convectionOnLine
     procedure :: constructor
  end interface convectionOnLine

contains

  type(ConvectionOnLineDT) function constructor(id, nodeIDList, coef, temp, node, geometry)
    implicit none
    integer(ikind)                 , intent(in) :: id
    integer(ikind)   , dimension(:), intent(in) :: nodeIDList
    real(rkind)                    , intent(in) :: coef
    real(rkind)                    , intent(in) :: temp
    type(NodePtrDT)  , dimension(:), intent(in) :: node
    class(GeometryDT), pointer     , intent(in) :: geometry
    call constructor%init(id, nodeIDList, coef, temp, node, geometry)
  end function constructor

  subroutine init(this, id, nodeIDList, coef, temp, node, geometry)
    implicit none
    class(ConvectionOnLineDT)      , intent(inout) :: this
    integer(ikind)                 , intent(in)    :: id
    integer(ikind)   , dimension(:), intent(in)    :: nodeIDList
    real(rkind)                    , intent(in)    :: coef
    real(rkind)                    , intent(in)    :: temp
    type(NodePtrDT)  , dimension(:), intent(in)    :: node
    class(GeometryDT), pointer     , intent(in)    :: geometry
    this%id = id
    this%affectsLHS = .true.
    this%affectsRHS = .true.
    this%nodeIDList = nodeIDList
    this%coef = coef
    this%temp = temp
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, k
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(PointDT)                                                         :: gaussPoint
    type(NodePtrDT)          , dimension(:)   , allocatable               :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%integrator
    allocate(lhs(nNode,nNode))
    allocate(rhs(nNode))
    allocate(jacobianDet(integrator%ptr%integTerms))
    allocate(nodalPoints(nNode))
    gaussPoint = point(0._rkind, 0._rkind)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    do i = 1, integrator%ptr%integTerms
       call gaussPoint%updatePoint(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%geometry%jacobianDet(this%nodeIDList, gaussPoint, nodalPoints)
    end do
    do i = 1, nNode
       do j = 1, nNode
          lhs(i,j) = 0._rkind
          int = 0._rkind
          do k = 1, integrator%ptr%integTerms
             int = int + integrator%ptr%weight(k)*integrator%ptr%shapeFunc(k,this%nodeIDList(i)) &
                  * integrator%ptr%shapeFunc(k,this%nodeIDList(j))*this%coef
          end do
          lhs(i,j) = lhs(i,j) + int
       end do
    end do
    do i = 1, nNode
       int = 0._rkind
       do j = 1, integrator%ptr%integTerms
          int = int + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,this%nodeIDList(i)) &
               * this%coef*this%temp*jacobianDet(j)
       end do
    end do
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, k
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(PointDT)                                                         :: gaussPoint
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%integrator
    allocate(lhs(nNode,nNode))
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    do i = 1, integrator%ptr%integTerms
       call gaussPoint%updatePoint(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%geometry%jacobianDet(this%nodeIDList, gaussPoint, nodalPoints)
    end do
    do i = 1, nNode
       do j = 1, nNode
          lhs(i,j) = 0._rkind
          int = 0._rkind
          do k = 1, integrator%ptr%integTerms
             int = int + integrator%ptr%weight(k)*integrator%ptr%shapeFunc(k,nodalPoints(i)%getID()) &
                  * integrator%ptr%shapeFunc(k,nodalPoints(j)%getID())*this%coef
          end do
          lhs(i,j) = lhs(i,j) + int
       end do
    end do
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(PointDT)                                                         :: gaussPoint
    type(NodePTrDT)          , dimension(:)    , allocatable              :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%integrator
    allocate(rhs(nNode))
    allocate(jacobianDet(integrator%ptr%integTerms))
    allocate(nodalPoints(nNode))
    gaussPoint = point(0._rkind, 0._rkind)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    do i = 1, integrator%ptr%integTerms
       call gaussPoint%updatePoint(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%geometry%jacobianDet(this%nodeIDList, gaussPoint, nodalPoints)
    end do
    do i = 1, nNode
       int = 0._rkind
       do j = 1, integrator%ptr%integTerms
          int = int + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,nodalPoints(i)%getID()) &
               * this%coef*this%temp*jacobianDet(j)
       end do
    end do
  end subroutine calculateRHS

end module ConvectionOnLineM
