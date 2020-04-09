module Triangle2D3NodeM
  use UtilitiesM
  use DebuggerM

  use PointM

  use IntegratorM
  
  use GeometryM

  implicit none
  
  private
  public :: Triangle2D3NodeDT, triangle2D3Node

  type, extends(GeometryDT) :: Triangle2D3NodeDT
     procedure, public  :: init
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: jacobian
     generic  , public  :: jacobianDet => jacobianDetFromCoord, jacobianDetFromJacobian
     procedure, private :: jacobianDetFromCoord
     procedure, private :: jacobianDetFromJacobian
  end type Triangle2D3NodeDT

  interface triangle2D3Node
     procedure :: constructor
  end interface triangle2D3Node

  integer(ikind), parameter :: NNODE = 3

contains

  type(Triangle2D3NodeDT) function constructor(node, integrator)
    implicit none
    type(NodePtrDT), dimension(NNODE), intent(in) :: node
    type(IntegratorDT)               , intent(in) :: integrator
    call constructor%init(node, integrator)
  end function constructor

  subroutine init(this, node, integrator)
    implicit none
    class(Triangle2D3NodeDT)         , intent(inout) :: this
    type(NodePtrDT), dimension(NNODE), intent(in)    :: node
    type(IntegratorDT)               , intent(in)    :: integrator
    this%nNode = NNODE
    this%node(1:NNODE) = node(1:NNODE)
    this%integrator = integrator
  end subroutine init

  function shapeFunc(this, u, v)
    implicit none
    class(Triangle2D3NodeDT), intent(inout)    :: this
    real(rkind)             , intent(in)       :: u
    real(rkind)             , intent(in)       :: v
    real(rkind)             , dimension(NNODE) :: shapeFunc
    shapeFunc(1) = 1 - u - v
    shapeFunc(2) = u
    shapeFunc(3) = v
  end function shapeFunc

  function dShapeFunc(this, u, v)
    implicit none
    class(Triangle2D3NodeDT), intent(inout)       :: this
    real(rkind)             , intent(in)          :: u
    real(rkind)             , intent(in)          :: v
    real(rkind)             , dimension(2, NNODE) :: dShapeFunc
    dShapeFunc(1,1) = -1
    dShapeFunc(1,2) = 1
    dShapeFunc(1,3) = 0
    dShapeFunc(2,1) = -1
    dShapeFunc(2,2) = 0
    dShapeFunc(2,3) = 1
  end function dShapeFunc

  function jacobian(this, u, v, point)
    implicit none
    class(Triangle2D3NodeDT), intent(inout)                :: this
    real(rkind)             , intent(in)                   :: u
    real(rkind)             , intent(in)                   :: v
    class(PointDT)          , dimension(NNODE), intent(in) :: point
    real(rkind)             , dimension(2,2)               :: jacobian
    integer(ikind)                                         :: i
    real(rkind)             , dimension(2, NNODE)          :: dsf
    jacobian = 0.d0
    dsf = this%dShapeFunc(u,v)
    do i = 1, NNODE
       jacobian(1,1) = jacobian(1,1) + dsf(1,i)*point(i)%getx() !dx/d(xi)
       jacobian(1,2) = jacobian(1,2) + dsf(1,i)*point(i)%gety() !dy/d(xi)
       jacobian(2,1) = jacobian(2,1) + dsf(2,i)*point(i)%getx() !dx/d(eta)
       jacobian(2,2) = jacobian(2,2) + dsf(2,i)*point(i)%gety() !dy/d(eta)
    end do
  end function jacobian

  real(rkind) function jacobianDetFromCoord(this, u, v, point)
    implicit none
    class(Triangle2D3NodeDT), intent(inout)                :: this
    real(rkind)             , intent(in)                   :: u
    real(rkind)             , intent(in)                   :: v
    class(PointDT)          , dimension(NNODE), intent(in) :: point
    real(rkind)                                            :: jacobianDetFromCoord
    real(rkind)             , dimension(2,2)               :: jacobian
    jacobian = this%jacobian(u,v,point)
    jacobianDet = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoord

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Triangle2D3NodeDT), intent(inout) :: this
    real(rkind), dimension(2,2), intent(in) :: jacobian
    real(rkind)                             :: jacobianDetFromJacobian
    jacobianDet = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromJacobian
  
end module Triangle2D3NodeM
