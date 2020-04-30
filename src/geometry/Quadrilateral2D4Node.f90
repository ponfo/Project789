module Quadrilateral2D4NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use NodeM
  use NodePtrM

  use IntegratorM
  
  use GeometryM

  implicit none
  
  private
  public :: Quadrilateral2D4NodeDT, quadrilateral2D4Node

  type, extends(GeometryDT) :: Quadrilateral2D4NodeDT
   contains
     procedure, public  :: init
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: jacobian
     procedure, public  :: jacobianDetFromCoord
     procedure, public  :: jacobianDetFromJacobian
     procedure, private :: valueShapeFuncAtGPoints
  end type Quadrilateral2D4NodeDT

  interface quadrilateral2D4Node
     procedure :: constructor
  end interface quadrilateral2D4Node

  integer(ikind), parameter :: NNODE = 4

contains

  type(Quadrilateral2D4NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout) :: this
    integer(ikind)               , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%integrator = integrator(gaussOrder, 'quadrilateral')
    call this%valueShapeFuncAtGPoints()
  end subroutine init

  function shapeFunc(this, point)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)         :: this
    class(PointDT)               , intent(in)            :: point
    real(rkind)                  , dimension(this%nNode) :: shapeFunc
    shapeFunc(1) = (1.d0/4.d0)*(1-point%getx())*(1-point%gety())
    shapeFunc(2) = (1.d0/4.d0)*(1+point%getx())*(1-point%gety())
    shapeFunc(3) = (1.d0/4.d0)*(1+point%getx())*(1+point%gety())
    shapeFunc(4) = (1.d0/4.d0)*(1-point%getx())*(1+point%gety())
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)                   :: this
    class(PointDT)               , intent(in)                      :: point
    real(rkind)                  , dimension(this%dim, this%nNode) :: dShapeFunc
    dShapeFunc(1,1) = -(1.d0/4.d0)*(1-point%gety())
    dShapeFunc(2,1) = -(1.d0/4.d0)*(1-point%getx())
    dShapeFunc(1,2) =  (1.d0/4.d0)*(1-point%gety())
    dShapeFunc(2,2) = -(1.d0/4.d0)*(1+point%getx())
    dShapeFunc(1,3) =  (1.d0/4.d0)*(1+point%gety())
    dShapeFunc(2,3) =  (1.d0/4.d0)*(1+point%getx())
    dShapeFunc(1,4) = -(1.d0/4.d0)*(1+point%gety())
    dShapeFunc(2,4) =  (1.d0/4.d0)*(1-point%getx())
  end function dShapeFunc

  function jacobian(this, pointToValue, nodalPoints)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)                     :: this
    class(PointDT)               , intent(in)                        :: pointToValue
    class(PointDT)               , dimension(this%nNode), intent(in) :: nodalPoints
    real(rkind)                  , dimension(this%dim, this%dim)      :: jacobian
    integer(ikind)                                                    :: i
    real(rkind)                  , dimension(2, NNODE)                :: dsf
    jacobian = 0.d0
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, NNODE
       jacobian(1,1) = jacobian(1,1) + dsf(1,i)*nodalPoints(i)%getx() !dx/d(xi)
       jacobian(1,2) = jacobian(1,2) + dsf(1,i)*nodalPoints(i)%gety() !dy/d(xi)
       jacobian(2,1) = jacobian(2,1) + dsf(2,i)*nodalPoints(i)%getx() !dx/d(eta)
       jacobian(2,2) = jacobian(2,2) + dsf(2,i)*nodalPoints(i)%gety() !dy/d(eta)
    end do
  end function jacobian

  real(rkind) function jacobianDetFromCoord(this, pointToValue, nodalPoints)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)                     :: this
    class(PointDT)               , intent(in)                        :: pointToValue
    class(PointDT)               , dimension(this%nNode), intent(in) :: nodalPoints
    real(rkind)                  , dimension(2,2)                    :: jacobian
    jacobian = this%jacobian(pointToValue, nodalPoints)
    jacobianDetFromCoord = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoord

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout) :: this
    real(rkind)  , dimension(:,:), intent(in)    :: jacobian
    jacobianDetFromJacobian = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromJacobian

  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout) :: this
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

end module Quadrilateral2D4NodeM
