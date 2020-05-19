module Quadrilateral3D8NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use NodeM
  use NodePtrM
  use Line3D3NodeM

  use IntegratorM
  
  use GeometryM

  implicit none
  
  private
  public :: Quadrilateral3D8NodeDT, quadrilateral3D8Node

  type, extends(GeometryDT) :: Quadrilateral3D8NodeDT
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
  end type Quadrilateral3D8NodeDT

  interface quadrilateral3D8Node
     procedure :: constructor
  end interface quadrilateral3D8Node

  integer(ikind), parameter :: NNODE = 8

contains

  type(Quadrilateral3D8NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout) :: this
    integer(ikind)               , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%dim = 3
    this%integrator = integrator(gaussOrder, 'quadrilateral')
    call this%valueShapeFuncAtGPoints()
    this%boundaryGeometry = line3D3Node(gaussOrder)
  end subroutine init

  real(rkind) function getLenght(this, node)
    implicit none
    class(Quadrilateral3D8NodeDT)                       , intent(inout) :: this
    class(NodePtrDT)             , dimension(this%nNode), intent(in)    :: node
    integer(ikind)                                                      :: i
    real(rkind)                  , dimension(:), allocatable            :: jacobianDet
    jacobianDet = this%jacobianDetAtGPoints(node)
    getLenght = 0._rkind
    do i = 1, this%integrator%getIntegTerms()
       getLenght = getLenght + this%integrator%getWeight(i)*jacobianDet(i)
    end do
  end function getLenght

  function shapeFunc(this, point)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout)         :: this
    class(PointDT)               , intent(in)            :: point
    real(rkind)                  , dimension(this%nNode) :: shapeFunc
    !Corners
    shapeFunc(1) = (1._rkind/4._rkind)*(1-point%gety())*(1-point%getx())*(-1-point%getx()-point%gety())
    shapeFunc(2) = (1._rkind/4._rkind)*(1-point%gety())*(1+point%getx())*(-1-point%gety()+point%getx())
    shapeFunc(3) = (1._rkind/4._rkind)*(1+point%gety())*(1+point%getx())*(-1+point%getx()+point%gety())
    shapeFunc(4) = (1._rkind/4._rkind)*(1+point%gety())*(1-point%getx())*(-1-point%getx()+point%gety())
    !Sides
    shapefunc(5) = (1._rkind/2._rkind)*(1-point%gety())*(1-point%getx()*point%getx())
    shapeFunc(6) = (1._rkind/2._rkind)*(1+point%getx())*(1-point%gety()*point%gety())
    shapeFunc(7) = (1._rkind/2._rkind)*(1+point%gety())*(1-point%getx()*point%getx())
    shapefunc(8) = (1._rkind/2._rkind)*(1-point%getx())*(1-point%gety()*point%gety())
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout)                   :: this
    class(PointDT)               , intent(in)                      :: point
    real(rkind)                  , dimension(this%dim, this%nNode) :: dShapeFunc
    !Corners
    dShapeFunc(1,1) = point%getx()/4._rkind+point%gety()/4._rkind                       &
         -(point%getx()*point%gety())/4._rkind-(point%gety()/4._rkind-1._rkind/4._rkind) &
         *(point%getx()+point%gety()+1)-1._rkind/4._rkind
    dShapeFunc(2,1) = point%getx()/4._rkind+point%gety()/4._rkind                       &
         -(point%getx()*point%gety())/4._rkind-(point%getx()/4._rkind-1._rkind/4._rkind) &
         *(point%getx()+point%gety()+1)-1._rkind/4._rkind
    dShapeFunc(1,2) = (point%gety()/4._rkind-1._rkind/4._rkind)*(point%gety()-point%getx()+1) &
         -(point%gety()/4._rkind-1._rkind/4._rkind)*(point%getx()+1)
    dShapeFunc(2,2) = (point%gety()/4._rkind-1._rkind/4._rkind)*(point%getx()+1) &
         +((point%getx()+1)*(point%gety()-point%getx()+1))/4._rkind
    dShapeFunc(1,3) = (point%gety()/4._rkind+1._rkind/4._rkind)*(point%getx()+1) &
         +(point%gety()/4._rkind+1._rkind/4._rkind)*(point%getx()+point%gety()-1)
    dShapeFunc(2,3) = (point%gety()/4._rkind+1._rkind/4._rkind)*(point%getx()+1) &
         +((point%getx()+1)*(point%getx()+point%gety()-1))/4._rkind
    dShapeFunc(1,4) = (point%gety()/4._rkind+1._rkind/4._rkind)*(point%getx() &
         -point%gety()+1)+(point%gety()/4._rkind+1/4._rkind)*(point%getx()-1)
    dShapeFunc(2,4) = ((point%getx()-1)*(point%getx()-point%gety()+1))/4._rkind &
         -(point%gety()/4._rkind+1._rkind/4._rkind)*(point%getx()-1)
    !Sides
    dShapeFunc(1,5) = 2*point%getx()*(point%gety()/2._rkind-1._rkind/2._rkind)
    dShapeFunc(2,5) = (point%getx()**2)/2._rkind-1._rkind/2._rkind
    dShapeFunc(1,6) = 1._rkind/2._rkind-(point%gety()**2)/2._rkind
    dShapeFunc(2,6) = -2*point%gety()*(point%getx()/2._rkind+1._rkind/2._rkind)
    dShapeFunc(1,7) = -2*point%getx()*(point%gety()/2._rkind+1._rkind/2._rkind)
    dShapeFunc(2,7) = 1._rkind/2._rkind-(point%getx()**2)/2._rkind
    dShapeFunc(1,8) = (point%gety()**2)/2._rkind-1._rkind/2._rkind
    dShapeFunc(2,8) = 2*point%gety()*(point%getx()/2._rkind-1._rkind/2._rkind)
  end function dShapeFunc

  function jacobianAllNodes(this, pointToValue, node)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout)                     :: this
    class(PointDT)               , intent(in)                        :: pointToValue
    class(NodePtrDT)             , dimension(this%nNode), intent(in) :: node
    real(rkind)                  , dimension(this%dim, this%dim)     :: jacobianAllNodes
    integer(ikind)                                                   :: i
    real(rkind)                 , dimension(3, NNODE)                :: dsf
    jacobianAllNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, this%nNode
       jacobianAllNodes(1,1) = jacobianAllNodes(1,1) + dsf(1,i)*node(i)%getx() !dx/d(xi)
       jacobianAllNodes(1,2) = jacobianAllNodes(1,2) + dsf(1,i)*node(i)%gety() !dy/d(xi)
       jacobianAllNodes(1,3) = jacobianAllNodes(1,3) + dsf(1,i)*node(i)%getz() !dz/d(xi)
       jacobianAllNodes(2,1) = jacobianAllNodes(2,1) + dsf(2,i)*node(i)%getx() !dx/d(eta)
       jacobianAllNodes(2,2) = jacobianAllNodes(2,2) + dsf(2,i)*node(i)%gety() !dy/d(eta)
       jacobianAllNodes(2,3) = jacobianAllNodes(2,3) + dsf(2,i)*node(i)%getz() !dz/d(eta)
    end do
  end function jacobianAllNodes

  function jacobianSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout)                     :: this
    integer(ikind)               , dimension(:)         , intent(in) :: indexList
    class(PointDT)               , intent(in)                        :: pointToValue
    class(NodePtrDT)             , dimension(:)         , intent(in) :: node
    real(rkind)                  , dimension(this%dim, this%dim)      :: jacobianSomeNodes
    integer(ikind)                                                   :: i
    real(rkind)                 , dimension(3, NNODE)                :: dsf
    jacobianSomeNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, size(node)
       jacobianSomeNodes(1,1) = jacobianSomeNodes(1,1) + dsf(1,indexList(i))*node(i)%getx() !dx/d(xi)
       jacobianSomeNodes(1,2) = jacobianSomeNodes(1,2) + dsf(1,indexList(i))*node(i)%gety() !dy/d(xi)
       jacobianSomeNodes(1,3) = jacobianSomeNodes(1,3) + dsf(1,indexList(i))*node(i)%getz() !dz/d(xi)
       jacobianSomeNodes(2,1) = jacobianSomeNodes(2,1) + dsf(2,indexList(i))*node(i)%getx() !dx/d(eta)
       jacobianSomeNodes(2,2) = jacobianSomeNodes(2,2) + dsf(2,indexList(i))*node(i)%gety() !dy/d(eta)
       jacobianSomeNodes(2,3) = jacobianSomeNodes(2,3) + dsf(2,indexList(i))*node(i)%getz() !dz/d(eta)
    end do
  end function jacobianSomeNodes

  function jacobianAtGPoints(this, node)
    implicit none
    class(Quadrilateral3D8NodeDT)                       , intent(inout) :: this
    class(NodePtrDT), dimension(this%nNode)             , intent(in)    :: node
    real(rkind), dimension(this%integrator%integTerms,this%dim,this%dim) :: jacobianAtGPoints
    integer(ikind)                                                      :: i, j
    real(rkind), dimension(this%integrator%integTerms,3,NNODE)          :: dsf
    jacobianAtGPoints = 0._rkind
    dsf = this%integrator%dShapeFunc
    do i = 1, this%integrator%integTerms
       do j = 1, this%nNode
          jacobianAtGPoints(i,1,1) = jacobianAtGPoints(i,1,1) + dsf(i,1,j)*node(j)%getx()
          jacobianAtGPoints(i,1,2) = jacobianAtGPoints(i,1,2) + dsf(i,1,j)*node(j)%gety()
          jacobianAtGPoints(i,1,3) = jacobianAtGPoints(i,1,3) + dsf(i,1,j)*node(j)%getz()
          jacobianAtGPoints(i,2,1) = jacobianAtGPoints(i,2,1) + dsf(i,2,j)*node(j)%getx()
          jacobianAtGPoints(i,2,2) = jacobianAtGPoints(i,2,2) + dsf(i,2,j)*node(j)%gety()
          jacobianAtGPoints(i,2,3) = jacobianAtGPoints(i,2,3) + dsf(i,2,j)*node(j)%getz()
       end do
    end do
  end function jacobianAtGPoints

  real(rkind) function jacobianDetFromCoordAllNodes(this, pointToValue, node)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout)                :: this
    class(PointDT)          , intent(in)                        :: pointToValue
    class(NodePtrDT)        , dimension(this%nNode), intent(in) :: node
    real(rkind)             , dimension(3,3)                    :: jacobian
    real(rkind)                                                 :: x, y, z
    jacobian = this%jacobian(pointToValue, node)
    x = jacobian(1,2)*jacobian(2,3)-jacobian(2,2)*jacobian(1,3)
    y = jacobian(2,1)*jacobian(1,3)-jacobian(1,1)*jacobian(2,3)
    z = jacobian(1,1)*jacobian(3,2)-jacobian(2,1)*jacobian(1,2)
    jacobianDetFromCoordAllNodes = sqrt(x*x+y*y+z*z)
  end function jacobianDetFromCoordAllNodes

  real(rkind) function jacobianDetFromCoordSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout)       :: this
    integer(ikind)          , dimension(:), intent(in) :: indexList
    class(PointDT)          , intent(in)               :: pointToValue
    class(NodePtrDT)        , dimension(:), intent(in) :: node
    real(rkind)             , dimension(3,3)           :: jacobian
    real(rkind)                                        :: x, y, z
    jacobian = this%jacobian(indexList, pointToValue, node)
    x = jacobian(1,2)*jacobian(2,3)-jacobian(2,2)*jacobian(1,3)
    y = jacobian(2,1)*jacobian(1,3)-jacobian(1,1)*jacobian(2,3)
    z = jacobian(1,1)*jacobian(3,2)-jacobian(2,1)*jacobian(1,2)
    jacobianDetFromCoordSomeNodes = sqrt(x*x+y*y+z*z)
  end function jacobianDetFromCoordSomeNodes

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout) :: this
    real(rkind), dimension(:,:)  , intent(in)    :: jacobian
    real(rkind)                                  :: x, y, z
    x = jacobian(1,2)*jacobian(2,3)-jacobian(2,2)*jacobian(1,3)
    y = jacobian(2,1)*jacobian(1,3)-jacobian(1,1)*jacobian(2,3)
    z = jacobian(1,1)*jacobian(3,2)-jacobian(2,1)*jacobian(1,2)
    jacobianDetFromJacobian = sqrt(x*x+y*y+z*z)
  end function jacobianDetFromJacobian

  function jacobianDetAtGPointsFromCoord(this, node)
    implicit none
    class(Quadrilateral3D8NodeDT)          , intent(inout)      :: this
    class(NodePtrDT), dimension(this%nNode), intent(in)         :: node
    real(rkind)     , dimension(this%integrator%integTerms)     :: jacobianDetAtGPointsFromCoord
    integer(ikind)                                              :: i
    real(rkind)     , dimension(this%integrator%integTerms,3,3) :: jacobian
    real(rkind)                                                 :: x, y, z
    jacobian = this%jacobianAtGPoints(node)
    do i = 1, this%integrator%integTerms
       x = jacobian(i,1,2)*jacobian(i,2,3)-jacobian(i,2,2)*jacobian(i,1,3)
       y = jacobian(i,2,1)*jacobian(i,1,3)-jacobian(i,1,1)*jacobian(i,2,3)
       z = jacobian(i,1,1)*jacobian(i,3,2)-jacobian(i,2,1)*jacobian(i,1,2)
       jacobianDetAtGPointsFromCoord(i) = sqrt(x*x+y*y+z*z)
    end do
  end function jacobianDetAtGPointsFromCoord

  function jacobianDetAtGPointsFromJacobian(this, jacobian)
    implicit none
    class(Quadrilateral3D8NodeDT)      , intent(inout) :: this
    real(rkind), dimension(:,:,:)        , intent(in) :: jacobian
    real(rkind), dimension(this%integrator%integTerms) :: jacobianDetAtGPointsFromJacobian
    integer(ikind)                                     :: i
    real(rkind)                                        :: x, y, z
    do i = 1, this%integrator%integTerms
       x = jacobian(i,1,2)*jacobian(i,2,3)-jacobian(i,2,2)*jacobian(i,1,3)
       y = jacobian(i,2,1)*jacobian(i,1,3)-jacobian(i,1,1)*jacobian(i,2,3)
       z = jacobian(i,1,1)*jacobian(i,3,2)-jacobian(i,2,1)*jacobian(i,1,2)
       jacobianDetAtGPointsFromJacobian(i) = sqrt(x*x+y*y+z*z)
    end do
  end function jacobianDetAtGPointsFromJacobian
  
  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Quadrilateral3D8NodeDT), intent(inout) :: this
    integer(ikind)                               :: i
    integer(ikind)                               :: integTerms
    real(rkind)                                  :: x
    real(rkind)                                  :: y
    integTerms = this%integrator%integTerms
    allocate(this%integrator%shapeFunc(integTerms, NNODE))
    allocate(this%integrator%dShapeFunc(integTerms, 3, NNODE))
    do i = 1, integTerms
       x = this%integrator%gPoint(i,1)
       y = this%integrator%gPoint(i,2)
       this%integrator%shapeFunc(i, 1:NNODE) = this%shapeFunc(point(x, y))
       this%integrator%dShapeFunc(i, 1:3, 1:NNODE) = this%dShapeFunc(point(x, y))
    end do
  end subroutine valueShapeFuncAtGPoints

end module Quadrilateral3D8NodeM
