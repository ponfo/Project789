module Quadrilateral2D4NodeM
  use UtilitiesM
  use DebuggerM
  
  use NodeM
  use NodePtrM

  use IntegratorM
  
  use GeometryM

  implicit none
  
  private
  public :: Quadrilateral2D4NodeDT, quadrilateral2D4Node

  type, extends(GeometryDT) :: Quadrilateral2D4NodeDT
     procedure, public  :: init
     procedure, public  :: shapeFunc
     procedure, public  :: dShapeFunc
     procedure, public  :: jacobian
     generic  , public  :: jacobianDet => jacobianDetFromCoord, jacobianDetFromJacobian
     procedure, private :: jacobianDetFromCoord
     procedure, private :: jacobianDetFromJacobian
  end type Quadrilateral2D4NodeDT

  interface quadrilateral2D4Node
     procedure :: constructor
  end interface quadrilateral2D4Node

  integer(ikind), parameter :: NNODE = 4

contains

  type(Quadrilateral2D4NodeDT) function constructor(node, integrator)
    implicit none
    type(NodePtrDT), dimension(NNODE), intent(in) :: node
    type(IntegratorDT)               , intent(in) :: integrator
    call constructor%init(node, integrator)
  end function constructor

  subroutine init(this, node, integrator)
    implicit none
    class(Quadrilateral2D4NodeDT)         , intent(inout) :: this
    type(NodePtrDT)     , dimension(NNODE), intent(in)    :: node
    type(IntegratorDT)                    , intent(in)    :: integrator
    this%nNode = NNODE
    this%node(1:NNODE) = node(1:NNODE)
    this%integrator = integrator
  end subroutine init

  function shapeFunc(this, u, v)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)    :: this
    real(rkind)                  , intent(in)       :: u
    real(rkind)                  , intent(in)       :: v
    real(rkind)                  , dimension(NNODE) :: shapeFunc
    shapeFunc(1) = (1.d0/4.d0)*(1-u)*(1-v)
    shapeFunc(2) = (1.d0/4.d0)*(1+u)*(1-v)
    shapeFunc(3) = (1.d0/4.d0)*(1+u)*(1+v)
    shapeFunc(4) = (1.d0/4.d0)*(1-u)*(1+v)
  end function shapeFunc

  function dShapeFunc(this, u, v)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)       :: this
    real(rkind)                  , intent(in)          :: u
    real(rkind)                  , intent(in)          :: v
    real(rkind)                  , dimension(2, NNODE) :: dShapeFunc
    dShapeFunc(1,1) = -(1.d0/4.d0)*(1-v)
    dShapeFunc(2,1) = -(1.d0/4.d0)*(1-u)
    dShapeFunc(1,2) =  (1.d0/4.d0)*(1-v)
    dShapeFunc(2,2) = -(1.d0/4.d0)*(1+u)
    dShapeFunc(1,3) =  (1.d0/4.d0)*(1+v)
    dShapeFunc(2,3) =  (1.d0/4.d0)*(1+u)
    dShapeFunc(1,4) = -(1.d0/4.d0)*(1+v)
    dShapeFunc(2,4) =  (1.d0/4.d0)*(1-u)
  end function dShapeFunc

  function jacobian(this, u, v)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)       :: this
    real(rkind)                  , intent(in)          :: u
    real(rkind)                  , intent(in)          :: v
    real(rkind)                  , dimension(2,2)      :: jacobian
    integer(ikind)                                     :: i
    real(rkind)                  , dimension(2, NNODE) :: dsf
    jacobian = 0.d0
    dsf = this%dShapeFunc(u,v)
    do i = 1, NNODE
       jacobian(1,1) = jacobian(1,1) + dsf(1,i)*this%node(i)%getx() !dx/d(xi)
       jacobian(1,2) = jacobian(1,2) + dsf(1,i)*this%node(i)%gety() !dy/d(xi)
       jacobian(2,1) = jacobian(2,1) + dsf(2,i)*this%node(i)%getx() !dx/d(eta)
       jacobian(2,2) = jacobian(2,2) + dsf(2,i)*this%node(i)%gety() !dy/d(eta)
    end do
  end function jacobian

  real(rkind) function jacobianDetFromCoord(this, u, v)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout)  :: this
    real(rkind)                  , intent(in)     :: u
    real(rkind)                  , intent(in)     :: v
    real(rkind)                                   :: jacobianDetFromCoord
    real(rkind)                  , dimension(2,2) :: jacobian
    jacobian = this%jacobian(u,v)
    jacobianDet = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromCoord

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Quadrilateral2D4NodeDT), intent(inout) :: this
    real(rkind)  , dimension(2,2), intent(in)    :: jacobian
    real(rkind)                                  :: jacobianDetFromJacobian
    jacobianDet = jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
  end function jacobianDetFromJacobian

end module Quadrilateral2D4NodeM
