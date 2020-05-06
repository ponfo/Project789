module Triangle2D6NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use Line2D3NodeM
  use NodeM
  use NodePtrM

  use IntegratorM

  use GeometryM

  implicit none

  private
  public :: Triangle2D6NodeDT, triangle2D6Node

  type, extends(GeometryDT) :: Triangle2D6NodeDT
   contains
     procedure, public  :: init
     procedure, public  :: getLenght
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: jacobianAllNodes
     procedure, public  :: jacobianSomeNodes
     procedure, public  :: jacobianAtGPoints
     procedure, public  :: jacobianDetFromCoordAllNodes
     procedure, public  :: jacobianDetFromCoordSomeNodes
     procedure, public  :: jacobianDetFromJacobian
     procedure, public  :: jacobianDetAtGPointsFromCoord
     procedure, public  :: jacobianDetAtGPointsFromJacobian
     procedure, private :: valueShapeFuncAtGPoints
  end type Triangle2D6NodeDT

  interface triangle2D6Node
     procedure :: constructor
  end interface triangle2D6Node

  integer(ikind), parameter :: NNODE = 6

contains

  type(Triangle2D6NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Triangle2D6NodeDT), intent(inout) :: this
    integer(ikind)          , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%dim = 2
    this%integrator = integrator(gaussOrder, 'triangle')
    call this%valueShapeFuncAtGPoints()
    this%boundaryGeometry = line2D3Node(gaussOrder)
  end subroutine init

  real(rkind) function getLenght(this, node)
    implicit none
    class(Triangle2D6NodeDT)                       , intent(inout) :: this
    class(NodePtrDT)        , dimension(this%nNode), intent(in)    :: node
    getLenght = node(1)%getx() * (node(2)%gety()-node(3)%gety()) &
         +      node(2)%getx() * (node(3)%gety()-node(1)%gety()) &
         +      node(3)%getx() * (node(1)%gety()-node(2)%gety())
  end function getLenght

  function shapeFunc(this, point)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)         :: this
    class(PointDT)          , intent(in)            :: point
    real(rkind)             , dimension(this%nNode) :: shapeFunc
    !Corners
    shapeFunc(1) = (1-point%getx()-point%gety())*(1-2*point%getx()-2*point%gety())
    shapeFunc(2) = point%getx()*(2*point%getx()-1)
    shapeFunc(3) = point%gety()*(2*point%gety()-1)
    !Sides
    shapeFunc(4) = 4*point%getx()*(1-point%getx()-point%gety())
    shapeFunc(5) = 4*point%getx()*point%gety()
    shapeFunc(6) = 4*point%gety()*(1-point%getx()-point%gety())
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)                   :: this
    class(PointDT)          , intent(in)                      :: point
    real(rkind)             , dimension(this%dim, this%nNode) :: dShapeFunc
    !Corners
    dShapeFunc(1,1) = 4*point%getx()+4*point%gety()-3
    dShapeFunc(2,1) = 4*point%gety()+4*point%getx()-3
    dShapeFunc(1,2) = 4*point%getx()-1
    dShapeFunc(2,2) = 0.d0
    dShapeFunc(1,3) = 0.d0
    dShapeFunc(2,3) = 4*point%gety()-1
    !Sides
    dShapeFunc(1,4) = -8*point%getx()-4*point%gety()+4
    dShapeFunc(2,4) = -4*point%getx()
    dShapeFunc(1,5) = 4*point%gety()
    dShapeFunc(2,5) = 4*point%getx()
    dShapeFunc(1,6) = -4*point%gety()
    dShapeFunc(2,6) = -8*point%gety()-4*point%getx()+4
  end function dShapeFunc

  function jacobianAllNodes(this, pointToValue, node)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)                     :: this
    class(PointDT)          , intent(in)                        :: pointToValue
    class(NodePtrDT)        , dimension(this%nNode), intent(in) :: node
    real(rkind)             , dimension(this%dim, this%dim)      :: jacobianAllNodes
    integer(ikind)                                               :: i
    real(rkind)             , dimension(2, NNODE)                :: dsf
    jacobianAllNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, this%nNode
       jacobianAllNodes(1,1) = jacobianAllNodes(1,1) + dsf(1,i)*node(i)%getx() !dx/d(xi)
       jacobianAllNodes(1,2) = jacobianAllNodes(1,2) + dsf(1,i)*node(i)%gety() !dy/d(xi)
       jacobianAllNodes(2,1) = jacobianAllNodes(2,1) + dsf(2,i)*node(i)%getx() !dx/d(eta)
       jacobianAllNodes(2,2) = jacobianAllNodes(2,2) + dsf(2,i)*node(i)%gety() !dy/d(eta)
    end do
  end function jacobianAllNodes

  function jacobianSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)                     :: this
    integer(ikind)          , dimension(:)         , intent(in) :: indexList
    class(PointDT)          , intent(in)                        :: pointToValue
    class(NodePtrDT)        , dimension(:)         , intent(in) :: node
    real(rkind)             , dimension(this%dim, this%dim)      :: jacobianSomeNodes
    integer(ikind)                                               :: i
    real(rkind)             , dimension(2, NNODE)                :: dsf
    jacobianSomeNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, size(node)
       jacobianSomeNodes(1,1) = jacobianSomeNodes(1,1) + dsf(1,indexList(i))*node(i)%getx() !dx/d(xi)
       jacobianSomeNodes(1,2) = jacobianSomeNodes(1,2) + dsf(1,indexList(i))*node(i)%gety() !dy/d(xi)
       jacobianSomeNodes(2,1) = jacobianSomeNodes(2,1) + dsf(2,indexList(i))*node(i)%getx() !dx/d(eta)
       jacobianSomeNodes(2,2) = jacobianSomeNodes(2,2) + dsf(2,indexList(i))*node(i)%gety() !dy/d(eta)
    end do
  end function jacobianSomeNodes

  function jacobianAtGPoints(this, node)
    implicit none
    class(Triangle2D6NodeDT)                            , intent(inout) :: this
    class(NodePtrDT), dimension(this%nNode)             , intent(in)    :: node
    real(rkind), dimension(this%integrator%integTerms,this%dim,this%dim) :: jacobianAtGPoints
    integer(ikind)                                                      :: i, j
    real(rkind), dimension(this%integrator%integTerms,2,NNODE)          :: dsf
    jacobianAtGPoints = 0._rkind
    dsf = this%integrator%dShapeFunc
    do i = 1, this%integrator%integTerms
       do j = 1, this%nNode
          jacobianAtGPoints(i,1,1) = jacobianAtGPoints(i,1,1) + dsf(i,1,j)*node(j)%getx()
          jacobianAtGPoints(i,1,2) = jacobianAtGPoints(i,1,2) + dsf(i,1,j)*node(j)%gety()
          jacobianAtGPoints(i,2,1) = jacobianAtGPoints(i,2,1) + dsf(i,2,j)*node(j)%getx()
          jacobianAtGPoints(i,2,2) = jacobianAtGPoints(i,2,2) + dsf(i,2,j)*node(j)%gety()
       end do
    end do
  end function jacobianAtGPoints

  real(rkind) function jacobianDetFromCoordAllNodes(this, pointToValue, node)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)                     :: this
    class(PointDT)          , intent(in)                        :: pointToValue
    class(NodePtrDT)        , dimension(this%nNode), intent(in) :: node
    real(rkind)             , dimension(2,2)                    :: jacobian
    jacobian = this%jacobian(pointToValue, node)
    jacobianDetFromCoordAllNodes = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoordAllNodes

  real(rkind) function jacobianDetFromCoordSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)            :: this
    integer(ikind)          , dimension(:), intent(in) :: indexList
    class(PointDT)          , intent(in)               :: pointToValue
    class(NodePtrDT)        , dimension(:), intent(in) :: node
    real(rkind)             , dimension(2,2)           :: jacobian
    jacobian = this%jacobian(indexList, pointToValue, node)
    jacobianDetFromCoordSomeNodes = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoordSomeNodes

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Triangle2D6NodeDT)   , intent(inout) :: this
    real(rkind), dimension(:,:), intent(in) :: jacobian
    jacobianDetFromJacobian = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromJacobian

  function jacobianDetAtGPointsFromCoord(this, node)
    implicit none
    class(Triangle2D6NodeDT)               , intent(inout)      :: this
    class(NodePtrDT), dimension(this%nNode), intent(in)         :: node
    real(rkind)     , dimension(this%integrator%integTerms)     :: jacobianDetAtGPointsFromCoord
    integer(ikind)                                              :: i
    real(rkind)     , dimension(this%integrator%integTerms,2,2) :: jacobian
    jacobian = this%jacobianAtGPoints(node)
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromCoord(i) = &
            jacobian(i,1,1)*jacobian(i,2,2)-jacobian(i,1,2)*jacobian(i,2,1)
    end do
  end function jacobianDetAtGPointsFromCoord

  function jacobianDetAtGPointsFromJacobian(this, jacobian)
    implicit none
    class(Triangle2D6NodeDT)          , intent(inout) :: this
    real(rkind), dimension(:,:,:)        , intent(in) :: jacobian
    real(rkind), dimension(this%integrator%integTerms) :: jacobianDetAtGPointsFromJacobian
    integer(ikind)                                     :: i
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromJacobian(i) = &
            jacobian(i,1,1)*jacobian(i,2,2)-jacobian(i,1,2)*jacobian(i,2,1)
    end do
  end function jacobianDetAtGPointsFromJacobian

  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Triangle2D6NodeDT), intent(inout) :: this
    integer(ikind)                          :: i
    integer(ikind)                          :: integTerms
    real(rkind)                             :: x
    real(rkind)                             :: y
    integTerms = this%integrator%integTerms
    allocate(this%integrator%shapeFunc(integTerms, NNODE))
    allocate(this%integrator%dShapeFunc(integTerms, 2, NNODE))
    do i = 1, integTerms
       x = this%integrator%gPoint(i,1)
       y = this%integrator%gPoint(i,2)
       this%integrator%shapeFunc(i, 1:NNODE) = this%shapeFunc(point(x, y))
       this%integrator%dShapeFunc(i, 1:2, 1:NNODE) = this%dShapeFunc(point(x, y))
    end do
  end subroutine valueShapeFuncAtGPoints

end module Triangle2D6NodeM
