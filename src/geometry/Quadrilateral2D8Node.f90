module Quadrilateral2D8NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use NodeM
  use NodePtrM

  use IntegratorM
  
  use GeometryM

  implicit none
  
  private
  public :: Quadrilateral2D8NodeDT, quadrilateral2D8Node

  type, extends(GeometryDT) :: Quadrilateral2D8NodeDT
   contains
     procedure, public  :: init
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: jacobianAllNodes
     procedure, public  :: jacobianSomeNodes
     procedure, public  :: jacobianDetFromCoordAllNodes
     procedure, public  :: jacobianDetFromCoordSomeNodes
     procedure, public  :: jacobianDetFromJacobian
     procedure, private :: valueShapeFuncAtGPoints
  end type Quadrilateral2D8NodeDT

  interface quadrilateral2D8Node
     procedure :: constructor
  end interface quadrilateral2D8Node

  integer(ikind), parameter :: NNODE = 8

contains

  type(Quadrilateral2D8NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout) :: this
    integer(ikind)               , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%integrator = integrator(gaussOrder, 'quadrilateral')
    call this%valueShapeFuncAtGPoints()
  end subroutine init

  function shapeFunc(this, point)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout)         :: this
    class(PointDT)               , intent(in)            :: point
    real(rkind)                  , dimension(this%nNode) :: shapeFunc
    !Corners
    shapeFunc(1) = (1.d0/4.d0)*(1-point%gety())*(1-point%getx())*(-1-point%getx()-point%gety())
    shapeFunc(2) = (1.d0/4.d0)*(1-point%gety())*(1+point%getx())*(-1-point%gety()+point%getx())
    shapeFunc(3) = (1.d0/4.d0)*(1+point%gety())*(1+point%getx())*(-1+point%getx()+point%gety())
    shapeFunc(4) = (1.d0/4.d0)*(1+point%gety())*(1-point%getx())*(-1-point%getx()+point%gety())
    !Sides
    shapefunc(5) = (1.d0/2.d0)*(1-point%gety())*(1-point%getx()*point%getx())
    shapeFunc(6) = (1.d0/2.d0)*(1+point%getx())*(1-point%gety()*point%gety())
    shapeFunc(7) = (1.d0/2.d0)*(1+point%gety())*(1-point%getx()*point%getx())
    shapefunc(8) = (1.d0/2.d0)*(1-point%getx())*(1-point%gety()*point%gety())
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout)                   :: this
    class(PointDT)               , intent(in)                      :: point
    real(rkind)                  , dimension(this%dim, this%nNode) :: dShapeFunc
    !Corners
    dShapeFunc(1,1) = point%getx()/4.d0+point%gety()/4.d0-(point%getx()*point%gety())/4.d0 &
         -(point%gety()/4.d0-1.d0/4.d0)*(point%getx()+point%gety()+1)-1.d0/4.d0
    dShapeFunc(2,1) = point%getx()/4.d0+point%gety()/4.d0-(point%getx()*point%gety())/4.d0 &
         -(point%getx()/4.d0-1.d0/4.d0)*(point%getx()+point%gety()+1)-1.d0/4.d0
    dShapeFunc(1,2) = (point%gety()/4.d0-1.d0/4.d0)*(point%gety()-point%getx()+1) &
         -(point%gety()/4.d0-1.d0/4.d0)*(point%getx()+1)
    dShapeFunc(2,2) = (point%gety()/4.d0-1.d0/4.d0)*(point%getx()+1) &
         +((point%getx()+1)*(point%gety()-point%getx()+1))/4.d0
    dShapeFunc(1,3) = (point%gety()/4.d0+1.d0/4.d0)*(point%getx()+1) &
         +(point%gety()/4.d0+1.d0/4.d0)*(point%getx()+point%gety()-1)
    dShapeFunc(2,3) = (point%gety()/4.d0+1.d0/4.d0)*(point%getx()+1) &
         +((point%getx()+1)*(point%getx()+point%gety()-1))/4.d0
    dShapeFunc(1,4) = (point%gety()/4.d0+1.d0/4.d0)*(point%getx() &
         -point%gety()+1)+(point%gety()/4.d0+1/4.d0)*(point%getx()-1)
    dShapeFunc(2,4) = ((point%getx()-1)*(point%getx()-point%gety()+1))/4.d0 &
         -(point%gety()/4.d0+1.d0/4.d0)*(point%getx()-1)
    !Sides
    dShapeFunc(1,5) = 2*point%getx()*(point%gety()/2.d0-1.d0/2.d0)
    dShapeFunc(2,5) = (point%getx()**2)/2.d0-1.d0/2.d0
    dShapeFunc(1,6) = 1.d0/2.d0-(point%gety()**2)/2.d0
    dShapeFunc(2,6) = -2*point%gety()*(point%getx()/2.d0+1.d0/2.d0)
    dShapeFunc(1,7) = -2*point%getx()*(point%gety()/2.d0+1.d0/2.d0)
    dShapeFunc(2,7) = 1.d0/2.d0-(point%getx()**2)/2.d0
    dShapeFunc(1,8) = (point%gety()**2)/2.d0-1.d0/2.d0
    dShapeFunc(2,8) = 2*point%gety()*(point%getx()/2.d0-1.d0/2.d0)
  end function dShapeFunc

  function jacobianAllNodes(this, pointToValue, node)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout)                     :: this
    class(PointDT)               , intent(in)                        :: pointToValue
    class(NodePtrDT)             , dimension(this%nNode), intent(in) :: node
    real(rkind)                  , dimension(this%dim, this%dim)     :: jacobianAllNodes
    integer(ikind)                                                   :: i
    real(rkind)                 , dimension(2, NNODE)                :: dsf
    jacobianAllNodes = 0.d0
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
    class(Quadrilateral2D8NodeDT), intent(inout)                     :: this
    integer(ikind)               , dimension(:)         , intent(in) :: indexList
    class(PointDT)               , intent(in)                        :: pointToValue
    class(NodePtrDT)             , dimension(:)         , intent(in) :: node
    real(rkind)                  , dimension(this%dim, this%dim)      :: jacobianSomeNodes
    integer(ikind)                                                   :: i
    real(rkind)                 , dimension(2, NNODE)                :: dsf
    jacobianSomeNodes = 0.d0
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, size(node)
       jacobianSomeNodes(1,1) = jacobianSomeNodes(1,1) + dsf(1,indexList(i))*node(i)%getx() !dx/d(xi)
       jacobianSomeNodes(1,2) = jacobianSomeNodes(1,2) + dsf(1,indexList(i))*node(i)%gety() !dy/d(xi)
       jacobianSomeNodes(2,1) = jacobianSomeNodes(2,1) + dsf(2,indexList(i))*node(i)%getx() !dx/d(eta)
       jacobianSomeNodes(2,2) = jacobianSomeNodes(2,2) + dsf(2,indexList(i))*node(i)%gety() !dy/d(eta)
    end do
  end function jacobianSomeNodes

  real(rkind) function jacobianDetFromCoordAllNodes(this, pointToValue, node)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout)                     :: this
    class(PointDT)               , intent(in)                        :: pointToValue
    class(NodePtrDT)             , dimension(this%nNode), intent(in) :: node
    real(rkind)                  , dimension(2,2)                    :: jacobian
    jacobian = this%jacobian(pointToValue, node)
    jacobianDetFromCoordAllNodes = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoordAllNodes

  real(rkind) function jacobianDetFromCoordSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout)            :: this
    integer(ikind)               , dimension(:), intent(in) :: indexList
    class(PointDT)               , intent(in)               :: pointToValue
    class(NodePtrDT)             , dimension(:), intent(in) :: node
    real(rkind)                  , dimension(2,2)           :: jacobian
    jacobian = this%jacobian(indexList, pointToValue, node)
    jacobianDetFromCoordSomeNodes = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoordSomeNodes

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout) :: this
    real(rkind)  , dimension(:,:), intent(in)    :: jacobian
    jacobianDetFromJacobian = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromJacobian

  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Quadrilateral2D8NodeDT), intent(inout) :: this
    integer(ikind)                               :: i
    integer(ikind)                               :: integTerms
    real(rkind)                                  :: x
    real(rkind)                                  :: y
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

end module Quadrilateral2D8NodeM
