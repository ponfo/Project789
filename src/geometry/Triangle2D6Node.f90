module Triangle2D6NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
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
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: jacobian
     generic  , public  :: jacobianDet => jacobianDetFromCoord, jacobianDetFromJacobian
     procedure, private :: jacobianDetFromCoord
     procedure, private :: jacobianDetFromJacobian
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
    this%integrator = integrator(gaussOrder, 'triangle')
    call this%valueShapeFuncAtGPoints()
  end subroutine init

  function shapeFunc(this, u, v)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)    :: this
    real(rkind)             , intent(in)       :: u
    real(rkind)             , intent(in)       :: v
    real(rkind)             , dimension(NNODE) :: shapeFunc
    !Corners
    shapeFunc(1) = (1-u-v)*(1-2*u-2*v)
    shapeFunc(2) = u*(2*u-1)
    shapeFunc(3) = v*(2*v-1)
    !Sides
    shapeFunc(4) = 4*u*(1-u-v)
    shapeFunc(5) = 4*u*v
    shapeFunc(6) = 4*v*(1-u-v)
  end function shapeFunc

  function dShapeFunc(this, u, v)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)       :: this
    real(rkind)             , intent(in)          :: u
    real(rkind)             , intent(in)          :: v
    real(rkind)             , dimension(2, NNODE) :: dShapeFunc
    !Corners
    dShapeFunc(1,1) = 4*u+4*v-3
    dShapeFunc(2,1) = 4*v+4*u-3
    dShapeFunc(1,2) = 4*u-1
    dShapeFunc(2,2) = 0.d0
    dShapeFunc(1,3) = 0.d0
    dShapeFunc(2,3) = 4*v-1
    !Sides
    dShapeFunc(1,4) = -8*u-4*v+4
    dShapeFunc(2,4) = -4*u
    dShapeFunc(1,5) = 4*v
    dShapeFunc(2,5) = 4*u
    dShapeFunc(1,6) = -4*v
    dShapeFunc(2,6) = -8*v-4*u+4
  end function dShapeFunc

  function jacobian(this, u, v, point)
    implicit none
    class(Triangle2D6NodeDT), intent(inout)                   :: this
    real(rkind)             , intent(in)                      :: u
    real(rkind)             , intent(in)                      :: v
    class(PointDT)          , dimension(NNODE), intent(inout) :: point
    real(rkind)             , dimension(2,2)                  :: jacobian
    integer(ikind)                                            :: i
    real(rkind)             , dimension(2, NNODE)             :: dsf
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
    class(Triangle2D6NodeDT), intent(inout)                   :: this
    real(rkind)             , intent(in)                      :: u
    real(rkind)             , intent(in)                      :: v
    class(PointDT)          , dimension(NNODE), intent(inout) :: point
    real(rkind)             , dimension(2,2)                  :: jacobian
    jacobian = this%jacobian(u,v,point)
    jacobianDetFromCoord = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoord

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Triangle2D6NodeDT), intent(inout) :: this
    real(rkind), dimension(2,2), intent(in) :: jacobian
    jacobianDetFromJacobian = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromJacobian

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
       this%integrator%shapeFunc(i, 1:NNODE) = this%shapeFunc(x, y)
       this%integrator%dShapeFunc(i, 1:2, 1:NNODE) = this%dShapeFunc(x, y)
    end do
  end subroutine valueShapeFuncAtGPoints

end module Triangle2D6NodeM
