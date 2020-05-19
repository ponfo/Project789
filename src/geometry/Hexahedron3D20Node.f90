module Hexahedron3D20NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use Quadrilateral3D8NodeM
  use NodePtrM

  use IntegratorM

  use GeometryM

  implicit none

  private
  public :: Hexahedron3D20NodeDT, hexahedron3D20Node

  type, extends(GeometryDT) :: Hexahedron3D20NodeDT
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
  end type Hexahedron3D20NodeDT

  interface hexahedron3D20Node
     procedure :: constructor
  end interface hexahedron3D20Node

  integer(ikind), parameter :: NNODE = 20

contains

  type(Hexahedron3D20NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout) :: this
    integer(ikind)             , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%dim = 3
    this%integrator = integrator(gaussOrder, 'hexahedron')
    call this%valueShapeFuncAtGPoints()
    this%boundaryGeometry = quadrilateral3D8Node(gaussOrder)
  end subroutine init

  real(rkind) function getLenght(this, node)
    implicit none
    class(Hexahedron3D20NodeDT)                       , intent(inout) :: this
    class(NodePtrDT)           , dimension(this%nNode), intent(in)    :: node
    getLenght = 0._rkind !Pendiente
  end function getLenght

  function shapeFunc(this, point)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout)         :: this
    class(PointDT)             , intent(in)            :: point
    real(rkind)                , dimension(this%nNode) :: shapeFunc
    real(rkind)                                        :: u, v, w
    u = point%getx()
    v = point%gety()
    w = point%getz()
    shapeFunc(1)  = -(1._rkind/8._rkind)*(1-u)*(1-v)*(1-w)*(2+u+v+w)
    shapeFunc(2)  = -(1._rkind/8._rkind)*(1+u)*(1-v)*(1-w)*(2-u+v+w)
    shapeFunc(3)  = -(1._rkind/8._rkind)*(1+u)*(1+v)*(1-w)*(2-u-v+w)
    shapeFunc(4)  = -(1._rkind/8._rkind)*(1-u)*(1+v)*(1-w)*(2+u-v+w)
    shapeFunc(5)  = -(1._rkind/8._rkind)*(1-u)*(1-v)*(1+w)*(2+u+v-w)
    shapeFunc(6)  = -(1._rkind/8._rkind)*(1+u)*(1-v)*(1+w)*(2-u+v-w)
    shapeFunc(7)  = -(1._rkind/8._rkind)*(1+u)*(1+v)*(1+w)*(2-u-v-w)
    shapeFunc(8)  = -(1._rkind/8._rkind)*(1-u)*(1+v)*(1+w)*(2+u-v-w)
    
    shapeFunc(9)  =  (1._rkind/4._rkind)*(1-u)*(1+u)*(1-v)*(1-w)
    shapeFunc(10) =  (1._rkind/4._rkind)*(1+u)*(1-v)*(1+v)*(1-w)
    shapeFunc(11) =  (1._rkind/4._rkind)*(1+u)*(1-u)*(1+v)*(1-w)
    shapeFunc(12) =  (1._rkind/4._rkind)*(1-u)*(1-v)*(1+v)*(1-w)
    
    shapeFunc(13) =  (1._rkind/4._rkind)*(1-u)*(1-v)*(1-w)*(1+w)
    shapeFunc(14) =  (1._rkind/4._rkind)*(1+u)*(1-v)*(1-w)*(1+w)
    shapeFunc(15) =  (1._rkind/4._rkind)*(1+u)*(1+v)*(1-w)*(1+w)
    shapeFunc(16) =  (1._rkind/4._rkind)*(1-u)*(1+v)*(1-w)*(1+w)

    shapeFunc(17) =  (1._rkind/4._rkind)*(1-u)*(1+u)*(1-v)*(1+w)
    shapeFunc(18) =  (1._rkind/4._rkind)*(1+u)*(1-v)*(1+v)*(1+w)
    shapeFunc(19) =  (1._rkind/4._rkind)*(1-u)*(1+u)*(1+v)*(1+w)
    shapeFunc(20) =  (1._rkind/4._rkind)*(1-u)*(1-v)*(1+v)*(1+w)
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout)                  :: this
    class(PointDT)             , intent(in)                     :: point
    real(rkind)                , dimension(this%dim, this%nNode) :: dShapeFunc
    real(rkind)                                                  :: u, v, w
    u = point%getx()
    v = point%gety()
    w = point%getz()
    dShapeFunc(1,1)  =  (1._rkind/8._rkind)*(v-1)*(w-1)*(2*u+v+w+1)
    dShapeFunc(2,1)  =  (1._rkind/8._rkind)*(u-1)*(w-1)*(u+2*v+w+1)
    dShapeFunc(3,1)  =  (1._rkind/8._rkind)*(u-1)*(v-1)*(u+v+2*w+1)
    dShapeFunc(1,2)  = -(1._rkind/8._rkind)*(v-1)*(w-1)*(-2*u+v+w+1)
    dShapeFunc(2,2)  =  (1._rkind/8._rkind)*(u+1)*(w-1)*(u-2*v-w-1)
    dShapeFunc(3,2)  =  (1._rkind/8._rkind)*(u+1)*(v-1)*(u-v-2*w-1)
    dShapeFunc(1,3)  = -(1._rkind/8._rkind)*(v+1)*(w-1)*(2*u+v-w-1)
    dShapeFunc(2,3)  = -(1._rkind/8._rkind)*(u+1)*(w-1)*(u+2*v-w-1)
    dShapeFunc(3,3)  = -(1._rkind/8._rkind)*(u+1)*(v+1)*(u+v-2*w-1)
    dShapeFunc(1,4)  =  (1._rkind/8._rkind)*(v+1)*(w-1)*(-2*u+v-w-1)
    dShapeFunc(2,4)  = -(1._rkind/8._rkind)*(u-1)*(w-1)*(u-2*v+w+1)
    dShapeFunc(3,4)  = -(1._rkind/8._rkind)*(u-1)*(v+1)*(u-v+2*w+1)
    dShapeFunc(1,5)  = -(1._rkind/8._rkind)*(v-1)*(w+1)*(2*u+v-w+1)
    dShapeFunc(2,5)  = -(1._rkind/8._rkind)*(u-1)*(w+1)*(u+2*v-w+1)
    dShapeFunc(3,5)  = -(1._rkind/8._rkind)*(u-1)*(v-1)*(u+v-2*w+1)
    dShapeFunc(1,6)  =  (1._rkind/8._rkind)*(v-1)*(w+1)*(-2*u+v-w+1)
    dShapeFunc(2,6)  = -(1._rkind/8._rkind)*(u+1)*(w+1)*(u-2*v+w-1)
    dShapeFunc(3,6)  = -(1._rkind/8._rkind)*(u+1)*(v-1)*(u-v+2*w-1)
    dShapeFunc(1,7)  =  (1._rkind/8._rkind)*(v+1)*(w+1)*(2*u+v+w-1)
    dShapeFunc(2,7)  =  (1._rkind/8._rkind)*(u+1)*(w+1)*(u+2*v+w-1)
    dShapeFunc(3,7)  =  (1._rkind/8._rkind)*(u+1)*(v+1)*(u+v+2*w-1)
    dShapeFunc(1,8)  = -(1._rkind/8._rkind)*(v+1)*(w+1)*(-2*u+v+w-1)
    dShapeFunc(2,8)  =  (1._rkind/8._rkind)*(u-1)*(w+1)*(u-2*v-w+1)
    dShapeFunc(3,8)  =  (1._rkind/8._rkind)*(u-1)*(v+1)*(u-v-2*w+1)
    dShapeFunc(1,9)  = -(1._rkind/2._rkind)*u*(v-1)*(w-1)
    dShapeFunc(2,9)  = -(1._rkind/4._rkind)*(u-1)*(u+1)*(w-1)
    dShapeFunc(3,9)  = -(1._rkind/4._rkind)*(u-1)*(u+1)*(v-1)
    dShapeFunc(1,10) =  (1._rkind/4._rkind)*(v-1)*(v+1)*(w-1)
    dShapeFunc(2,10) =  (1._rkind/2._rkind)*(u+1)*v*(w-1)
    dShapeFunc(3,10) =  (1._rkind/4._rkind)*(u+1)*(v-1)*(v+1)
    dShapeFunc(1,11) =  (1._rkind/2._rkind)*u*(v+1)*(w-1)
    dShapeFunc(2,11) =  (1._rkind/4._rkind)*(u-1)*(u+1)*(w-1)
    dShapeFunc(3,11) =  (1._rkind/4._rkind)*(u-1)*(u+1)*(v+1)
    dShapeFunc(1,12) = -(1._rkind/4._rkind)*(v-1)*(v+1)*(w-1)
    dShapeFunc(2,12) = -(1._rkind/2._rkind)*(u-1)*v*(w-1)
    dShapeFunc(3,12) = -(1._rkind/4._rkind)*(u-1)*(v-1)*(v+1)
    
    dShapeFunc(1,17) =  (1._rkind/2._rkind)*u*(v-1)*(w+1)
    dShapeFunc(2,17) =  (1._rkind/4._rkind)*(u-1)*(u+1)*(w+1)
    dShapeFunc(3,17) =  (1._rkind/4._rkind)*(u-1)*(u+1)*(v-1)
    dShapeFunc(1,18) = -(1._rkind/4._rkind)*(v-1)*(v+1)*(w+1)
    dShapeFunc(2,18) = -(1._rkind/2._rkind)*(u+1)*v*(w+1)
    dShapeFunc(3,18) = -(1._rkind/4._rkind)*(u+1)*(v-1)*(v+1)
    dShapeFunc(1,19) = -(1._rkind/2._rkind)*u*(v+1)*(w+1)
    dShapeFunc(2,19) = -(1._rkind/4._rkind)*(u-1)*(u+1)*(w+1)
    dShapeFunc(3,19) = -(1._rkind/4._rkind)*(u-1)*(u+1)*(v+1)
    dShapeFunc(1,20) =  (1._rkind/4._rkind)*(v-1)*(v+1)*(w+1)
    dShapeFunc(2,20) =  (1._rkind/2._rkind)*(u-1)*v*(w+1)
    dShapeFunc(3,20) =  (1._rkind/4._rkind)*(u-1)*(v-1)*(v+1)
    
    dShapeFunc(1,13) = -(1._rkind/4._rkind)*(v-1)*(w-1)*(w+1)
    dShapeFunc(2,13) = -(1._rkind/4._rkind)*(u-1)*(w-1)*(w+1)
    dShapeFunc(3,13) = -(1._rkind/2._rkind)*(u-1)*(v-1)*w
    dShapeFunc(1,14) =  (1._rkind/4._rkind)*(v-1)*(w-1)*(w+1)
    dShapeFunc(2,14) =  (1._rkind/4._rkind)*(u+1)*(w-1)*(w+1)
    dShapeFunc(3,14) =  (1._rkind/2._rkind)*(u+1)*(v-1)*w
    dShapeFunc(1,15) = -(1._rkind/4._rkind)*(v+1)*(w-1)*(w+1)
    dShapeFunc(2,15) = -(1._rkind/4._rkind)*(u+1)*(w-1)*(w+1)
    dShapeFunc(3,15) = -(1._rkind/2._rkind)*(u+1)*(v+1)*w
    dShapeFunc(1,16) =  (1._rkind/4._rkind)*(v+1)*(w-1)*(w+1)
    dShapeFunc(2,16) =  (1._rkind/4._rkind)*(u-1)*(w-1)*(w+1)
    dShapeFunc(3,16) =  (1._rkind/2._rkind)*(u-1)*(v+1)*w
  end function dShapeFunc

  function jacobianAllNodes(this, pointToValue, node)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout)                     :: this
    class(PointDT)             , intent(in)                        :: pointToValue
    class(NodePtrDT)           , dimension(this%nNode), intent(in) :: node
    real(rkind)                , dimension(this%dim, this%dim)     :: jacobianAllNodes
    integer(ikind)                                                 :: i
    real(rkind)                , dimension(3, NNODE)               :: dsf
    jacobianAllNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, this%nNode
       jacobianAllNodes(1,1) = jacobianAllNodes(1,1) + dsf(1,i)*node(i)%getx() !dx/d(xi)
       jacobianAllNodes(1,2) = jacobianAllNodes(1,2) + dsf(1,i)*node(i)%gety() !dy/d(xi)
       jacobianAllNodes(1,3) = jacobianAllNodes(1,3) + dsf(1,i)*node(i)%getz() !dz/d(xi)
       jacobianAllNodes(2,1) = jacobianAllNodes(2,1) + dsf(2,i)*node(i)%getx() !dx/d(eta)
       jacobianAllNodes(2,2) = jacobianAllNodes(2,2) + dsf(2,i)*node(i)%gety() !dy/d(eta)
       jacobianAllNodes(2,3) = jacobianAllNodes(2,3) + dsf(2,i)*node(i)%getz() !dz/d(eta)
       jacobianAllNodes(3,1) = jacobianAllNodes(3,1) + dsf(3,i)*node(i)%getx() !dx/d(zeta)
       jacobianAllNodes(3,2) = jacobianAllNodes(3,2) + dsf(3,i)*node(i)%gety() !dy/d(zeta)
       jacobianAllNodes(3,3) = jacobianAllNodes(3,3) + dsf(3,i)*node(i)%getz() !dz/d(zeta)
    end do
  end function jacobianAllNodes

  function jacobianSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout)                     :: this
    integer(ikind)             , dimension(:)         , intent(in) :: indexList
    class(PointDT)             , intent(in)                        :: pointToValue
    class(NodePtrDT)           , dimension(:)         , intent(in) :: node
    real(rkind)                , dimension(this%dim, this%dim)      :: jacobianSomeNodes
    integer(ikind)                                                  :: i
    real(rkind)                , dimension(3, NNODE)                :: dsf
    jacobianSomeNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, size(node)
       jacobianSomeNodes(1,1) = jacobianSomeNodes(1,1) + dsf(1,indexList(i))*node(i)%getx() !dx/d(xi)
       jacobianSomeNodes(1,2) = jacobianSomeNodes(1,2) + dsf(1,indexList(i))*node(i)%gety() !dy/d(xi)
       jacobianSomeNodes(1,3) = jacobianSomeNodes(1,3) + dsf(1,indexList(i))*node(i)%getz() !dz/d(xi)
       jacobianSomeNodes(2,1) = jacobianSomeNodes(2,1) + dsf(2,indexList(i))*node(i)%getx() !dx/d(eta)
       jacobianSomeNodes(2,2) = jacobianSomeNodes(2,2) + dsf(2,indexList(i))*node(i)%gety() !dy/d(eta)
       jacobianSomeNodes(2,3) = jacobianSomeNodes(2,3) + dsf(2,indexList(i))*node(i)%getz() !dz/d(eta)
       jacobianSomeNodes(3,1) = jacobianSomeNodes(3,1) + dsf(3,indexList(i))*node(i)%getx() !dx/d(zeta)
       jacobianSomeNodes(3,2) = jacobianSomeNodes(3,2) + dsf(3,indexList(i))*node(i)%gety() !dy/d(zeta)
       jacobianSomeNodes(3,3) = jacobianSomeNodes(3,3) + dsf(3,indexList(i))*node(i)%getz() !dz/d(zeta)
    end do
  end function jacobianSomeNodes

  function jacobianAtGPoints(this, node)
    implicit none
    class(Hexahedron3D20NodeDT)                         , intent(inout) :: this
    class(NodePtrDT), dimension(this%nNode)             , intent(in)    :: node
    real(rkind), dimension(this%integrator%integTerms,this%dim,this%dim) :: jacobianAtGPoints
    integer(ikind)                                                      :: i, j
    real(rkind), dimension(this%integrator%integTerms,3,NNODE)          :: dsf
    jacobianAtGPoints = 0._rkind
    dsf = this%integrator%dShapeFunc
    do i = 1, this%integrator%integTerms
       do j = 1, this%nNode
          jacobianAtGPoints(i,1,1) = jacobianAtGPoints(i,1,1) + dsf(i,1,j)*node(j)%getx() !dx/d(xi)
          jacobianAtGPoints(i,1,2) = jacobianAtGPoints(i,1,2) + dsf(i,1,j)*node(j)%gety() !dy/d(xi)
          jacobianAtGPoints(i,1,3) = jacobianAtGPoints(i,1,3) + dsf(i,1,j)*node(j)%getz() !dz/d(xi)
          jacobianAtGPoints(i,2,1) = jacobianAtGPoints(i,2,1) + dsf(i,2,j)*node(j)%getx() !dx/d(eta)
          jacobianAtGPoints(i,2,2) = jacobianAtGPoints(i,2,2) + dsf(i,2,j)*node(j)%gety() !dy/d(eta)
          jacobianAtGPoints(i,2,3) = jacobianAtGPoints(i,2,3) + dsf(i,2,j)*node(j)%getz() !dz/d(eta)
          jacobianAtGPoints(i,3,1) = jacobianAtGPoints(i,3,1) + dsf(i,3,j)*node(j)%getx() !dx/d(zeta)
          jacobianAtGPoints(i,3,2) = jacobianAtGPoints(i,3,2) + dsf(i,3,j)*node(j)%gety() !dy/d(zeta)
          jacobianAtGPoints(i,3,3) = jacobianAtGPoints(i,3,3) + dsf(i,3,j)*node(j)%getz() !dz/d(zeta)
       end do
    end do
  end function jacobianAtGPoints

  real(rkind) function jacobianDetFromCoordAllNodes(this, pointToValue, node)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout)                   :: this
    class(PointDT)          , intent(in)                        :: pointToValue
    class(NodePtrDT)        , dimension(this%nNode), intent(in) :: node
    real(rkind)             , dimension(3,3)                    :: jacobian
    jacobian = this%jacobian(pointToValue, node)
    jacobianDetFromCoordAllNodes = jacobian(1,1)*jacobian(2,2)*jacobian(3,3) &
         + jacobian(1,2)*jacobian(2,3)*jacobian(3,1)        &
         + jacobian(1,3)*jacobian(2,1)*jacobian(3,2)        &
         - jacobian(3,1)*jacobian(2,2)*jacobian(1,3)        &
         - jacobian(2,1)*jacobian(1,2)*jacobian(3,3)        &
         - jacobian(1,1)*jacobian(3,2)*jacobian(2,3)
  end function jacobianDetFromCoordAllNodes

  real(rkind) function jacobianDetFromCoordSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout)            :: this
    integer(ikind)          , dimension(:), intent(in) :: indexList
    class(PointDT)          , intent(in)               :: pointToValue
    class(NodePtrDT)        , dimension(:), intent(in) :: node
    real(rkind)             , dimension(3,3)           :: jacobian
    jacobian = this%jacobian(indexList, pointToValue, node)
    jacobianDetFromCoordSomeNodes = jacobian(1,1)*jacobian(2,2)*jacobian(3,3) &
         + jacobian(1,2)*jacobian(2,3)*jacobian(3,1)        &
         + jacobian(1,3)*jacobian(2,1)*jacobian(3,2)        &
         - jacobian(3,1)*jacobian(2,2)*jacobian(1,3)        &
         - jacobian(2,1)*jacobian(1,2)*jacobian(3,3)        &
         - jacobian(1,1)*jacobian(3,2)*jacobian(2,3)
  end function jacobianDetFromCoordSomeNodes

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Hexahedron3D20NodeDT)   , intent(inout) :: this
    real(rkind), dimension(:,:), intent(in)    :: jacobian
    jacobianDetFromJacobian = jacobian(1,1)*jacobian(2,2)*jacobian(3,3) &
         + jacobian(1,2)*jacobian(2,3)*jacobian(3,1)        &
         + jacobian(1,3)*jacobian(2,1)*jacobian(3,2)        &
         - jacobian(3,1)*jacobian(2,2)*jacobian(1,3)        &
         - jacobian(2,1)*jacobian(1,2)*jacobian(3,3)        &
         - jacobian(1,1)*jacobian(3,2)*jacobian(2,3)
  end function jacobianDetFromJacobian

  function jacobianDetAtGPointsFromCoord(this, node)
    implicit none
    class(Hexahedron3D20NodeDT)               , intent(inout)      :: this
    class(NodePtrDT), dimension(this%nNode), intent(in)         :: node
    real(rkind)     , dimension(this%integrator%integTerms)     :: jacobianDetAtGPointsFromCoord
    integer(ikind)                                              :: i
    real(rkind)     , dimension(this%integrator%integTerms,3,3) :: jacobian
    jacobian = this%jacobianAtGPoints(node)
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromCoord(i) = jacobian(i,1,1)*jacobian(i,2,2)*jacobian(i,3,3) &
         + jacobian(i,1,2)*jacobian(i,2,3)*jacobian(i,3,1)        &
         + jacobian(i,1,3)*jacobian(i,2,1)*jacobian(i,3,2)        &
         - jacobian(i,3,1)*jacobian(i,2,2)*jacobian(i,1,3)        &
         - jacobian(i,2,1)*jacobian(i,1,2)*jacobian(i,3,3)        &
         - jacobian(i,1,1)*jacobian(i,3,2)*jacobian(i,2,3)
    end do
  end function jacobianDetAtGPointsFromCoord

  function jacobianDetAtGPointsFromJacobian(this, jacobian)
    implicit none
    class(Hexahedron3D20NodeDT)          , intent(inout) :: this
    real(rkind), dimension(:,:,:)        , intent(in) :: jacobian
    real(rkind), dimension(this%integrator%integTerms) :: jacobianDetAtGPointsFromJacobian
    integer(ikind)                                     :: i
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromJacobian(i) = jacobian(i,1,1)*jacobian(i,2,2)*jacobian(i,3,3) &
         + jacobian(i,1,2)*jacobian(i,2,3)*jacobian(i,3,1)        &
         + jacobian(i,1,3)*jacobian(i,2,1)*jacobian(i,3,2)        &
         - jacobian(i,3,1)*jacobian(i,2,2)*jacobian(i,1,3)        &
         - jacobian(i,2,1)*jacobian(i,1,2)*jacobian(i,3,3)        &
         - jacobian(i,1,1)*jacobian(i,3,2)*jacobian(i,2,3)
    end do
  end function jacobianDetAtGPointsFromJacobian

  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Hexahedron3D20NodeDT), intent(inout) :: this
    integer(ikind)                             :: i
    integer(ikind)                             :: integTerms
    real(rkind)                                :: x
    real(rkind)                                :: y
    real(rkind)                                :: z
    integTerms = this%integrator%integTerms
    allocate(this%integrator%shapeFunc(integTerms, NNODE))
    allocate(this%integrator%dShapeFunc(integTerms, 3, NNODE))
    do i = 1, integTerms
       x = this%integrator%gPoint(i,1)
       y = this%integrator%gPoint(i,2)
       z = this%integrator%gPoint(i,3)
       this%integrator%shapeFunc(i, 1:NNODE) = this%shapeFunc(point(x, y, z))
       this%integrator%dShapeFunc(i, 1:3, 1:NNODE) = this%dShapeFunc(point(x, y, z))
    end do
  end subroutine valueShapeFuncAtGPoints

end module Hexahedron3D20NodeM
