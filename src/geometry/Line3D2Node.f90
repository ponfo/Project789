module Line3D2NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use NodePtrM

  use IntegratorM

  use GeometryM

  implicit none

  private
  public :: Line3D2NodeDT, line3D2Node

  type, extends(GeometryDT) :: Line3D2NodeDT
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
  end type Line3D2NodeDT

  interface line3D2Node
     procedure :: constructor
  end interface line3D2Node

  integer(ikind), parameter :: NNODE = 2

contains

  type(Line3D2NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Line3D2NodeDT), intent(inout) :: this
    integer(ikind)      , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%dim = 3
    this%integrator = integrator(gaussOrder, 'line')
    call this%valueShapeFuncAtGPoints()
  end subroutine init

  real(rkind) function getLenght(this, node)
    implicit none
    class(Line3D2NodeDT)                       , intent(inout) :: this
    class(NodePtrDT)    , dimension(this%nNode), intent(in)    :: node
    getLenght = sqrt((node(1)%getx()-node(2)%getx())**2 &
         + (node(1)%gety()-node(2)%gety())**2           &
         + (node(1)%getz()-node(2)%getz())**2           )
  end function getLenght

  function shapeFunc(this, point)
    implicit none
    class(Line3D2NodeDT), intent(inout)         :: this
    class(PointDT)      , intent(in)            :: point
    real(rkind)         , dimension(this%nNode) :: shapeFunc
    shapeFunc(1) = .5 - .5*point%getx()
    shapeFunc(2) = .5 + .5*point%getx()
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Line3D2NodeDT), intent(inout)                   :: this
    class(PointDT)      , intent(in)                      :: point
    real(rkind)         , dimension(this%dim, this%nNode) :: dShapeFunc
    dShapeFunc(1,1) = -.5
    dShapeFunc(1,2) = .5
    dShapeFunc(2,1) = 0._rkind
    dShapeFunc(2,2) = 0._rkind
  end function dShapeFunc
  
  function jacobianAllNodes(this, pointToValue, node)
    implicit none
    class(Line3D2NodeDT), intent(inout)                     :: this
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(this%nNode), intent(in) :: node
    real(rkind)         , dimension(this%dim, this%dim)     :: jacobianAllNodes
    integer(ikind)                                          :: i
    real(rkind)         , dimension(2, NNODE)               :: dsf
    jacobianAllNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, this%nNode
       jacobianAllNodes(1,1) = jacobianAllNodes(1,1) + dsf(1,i)*node(i)%getx()
       jacobianAllNodes(1,2) = jacobianAllNodes(1,2) + dsf(1,i)*node(i)%gety()
       jacobianAllNodes(1,3) = jacobianAllNodes(1,3) + dsf(1,i)*node(i)%getz()
    end do
  end function jacobianAllNodes

  function jacobianSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Line3D2NodeDT), intent(inout)                     :: this
    integer(ikind)      , dimension(:)         , intent(in) :: indexList
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(:)         , intent(in) :: node
    real(rkind)         , dimension(this%dim, this%dim)      :: jacobianSomeNodes
    integer(ikind)                                           :: i
    real(rkind)         , dimension(2, NNODE)                :: dsf
    jacobianSomeNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, size(node)
       jacobianSomeNodes(1,1) = jacobianSomeNodes(1,1) + dsf(1,i)*node(i)%getx()
       jacobianSomeNodes(1,2) = jacobianSomeNodes(1,2) + dsf(1,i)*node(i)%gety()
       jacobianSomeNodes(1,3) = jacobianSomeNodes(1,3) + dsf(1,i)*node(i)%getz()
    end do
  end function jacobianSomeNodes

  function jacobianAtGPoints(this, node)
    implicit none
    class(Line3D2NodeDT)                                , intent(inout) :: this
    class(NodePtrDT), dimension(this%nNode)             , intent(in)    :: node
    real(rkind), dimension(this%integrator%integTerms,this%dim,this%dim) :: jacobianAtGPoints
    integer(ikind)                                                      :: i, j
    real(rkind), dimension(this%integrator%integTerms,2,NNODE)           :: dsf
    jacobianAtGPoints = 0._rkind
    dsf = this%integrator%dShapeFunc
    do i = 1, this%integrator%integTerms
       do j = 1, this%nNode
          jacobianAtGPoints(i,1,1) = jacobianAtGPoints(i,1,1) + dsf(i,1,j)*node(j)%getx()
          jacobianAtGPoints(i,1,2) = jacobianAtGPoints(i,1,2) + dsf(i,1,j)*node(j)%gety()
          jacobianAtGPoints(i,1,3) = jacobianAtGPoints(i,1,3) + dsf(i,1,j)*node(j)%getz()
       end do
    end do
  end function jacobianAtGPoints

  real(rkind) function jacobianDetFromCoordAllNodes(this, pointToValue, node)
    implicit none
    class(Line3D2NodeDT), intent(inout)                     :: this
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(this%nNode), intent(in) :: node
    real(rkind)         , dimension(3,3)                    :: jacobian
    jacobian = this%jacobian(pointToValue, node)
    jacobianDetFromCoordAllNodes = sqrt(jacobian(1,1)**2 + jacobian(1,2)**2 + jacobian(1,3)**2)
  end function jacobianDetFromCoordAllNodes

  real(rkind) function jacobianDetFromCoordSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Line3D2NodeDT), intent(inout)                     :: this
    integer(ikind)      , dimension(:), intent(in)          :: indexList
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(:), intent(in)          :: node
    real(rkind)         , dimension(3,3)                    :: jacobian
    jacobian = this%jacobian(indexList, pointToValue, node)
    jacobianDetFromCoordSomeNodes = sqrt(jacobian(1,1)**2 + jacobian(1,2)**2 + jacobian(1,3)**2)
  end function jacobianDetFromCoordSomeNodes

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Line3D2NodeDT)       , intent(inout) :: this
    real(rkind), dimension(:,:), intent(in)    :: jacobian
    jacobianDetFromJacobian = sqrt(jacobian(1,1)**2 + jacobian(1,2)**2 + jacobian(1,3)**2)
  end function jacobianDetFromJacobian

  function jacobianDetAtGPointsFromCoord(this, node)
    implicit none
    class(Line3D2NodeDT)                   , intent(inout)      :: this
    class(NodePtrDT), dimension(this%nNode), intent(in)         :: node
    real(rkind)     , dimension(this%integrator%integTerms)     :: jacobianDetAtGPointsFromCoord
    integer(ikind)                                              :: i
    real(rkind)     , dimension(this%integrator%integTerms,3,3) :: jacobian
    jacobian = this%jacobianAtGPoints(node)
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromCoord(i) = &
            sqrt(jacobian(i,1,1)**2 + jacobian(i,1,2)**2 + jacobian(i,1,3)**2)
    end do
  end function jacobianDetAtGPointsFromCoord

  function jacobianDetAtGPointsFromJacobian(this, jacobian)
    implicit none
    class(Line3D2NodeDT)               , intent(inout) :: this
    real(rkind), dimension(:,:,:)      , intent(in)    :: jacobian
    real(rkind), dimension(this%integrator%integTerms) :: jacobianDetAtGPointsFromJacobian
    integer(ikind)                                     :: i
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromJacobian(i) = &
            sqrt(jacobian(i,1,1)**2 + jacobian(i,1,2)**2 + jacobian(i,1,3)**2)
    end do
  end function jacobianDetAtGPointsFromJacobian

  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Line3D2NodeDT), intent(inout) :: this
    integer(ikind)                      :: i
    integer(ikind)                      :: integTerms
    real(rkind)                         :: x
    integTerms = this%integrator%integTerms
    allocate(this%integrator%shapeFunc(integTerms, NNODE))
    allocate(this%integrator%dShapeFunc(integTerms, 2, NNODE))
    do i = 1, integTerms
       x = this%integrator%gPoint(i,1)
       this%integrator%shapeFunc(i, 1:NNODE) = this%shapeFunc(point(x))
       this%integrator%dShapeFunc(i, 1:2, 1:NNODE) = this%dShapeFunc(point(x))
    end do
  end subroutine valueShapeFuncAtGPoints

end module Line3D2NodeM
