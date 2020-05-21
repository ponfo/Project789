module Line1D3NodeM
  use UtilitiesM
  use DebuggerM

  use PointM
  use NodePtrM

  use IntegratorM

  use GeometryM

  implicit none

  private
  public :: Line1D3NodeDT, line1D3Node

  type, extends(GeometryDT) :: Line1D3NodeDT
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
  end type Line1D3NodeDT

  interface line1D3Node
     procedure :: constructor
  end interface line1D3Node

  integer(ikind), parameter :: NNODE = 3

contains

  type(Line1D3NodeDT) function constructor(gaussOrder)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    call constructor%init(gaussOrder)
  end function constructor

  subroutine init(this, gaussOrder)
    implicit none
    class(Line1D3NodeDT), intent(inout) :: this
    integer(ikind)      , intent(in)    :: gaussOrder
    this%nNode = NNODE
    this%dim = 1
    this%integrator = integrator(gaussOrder, 'line')
    call this%valueShapeFuncAtGPoints()
  end subroutine init

  real(rkind) function getLenght(this, node)
    implicit none
    class(Line1D3NodeDT)                       , intent(inout) :: this
    class(NodePtrDT)    , dimension(this%nNode), intent(in)    :: node
    getLenght = sqrt((node(1)%getx()-node(2)%getx())**2)
  end function getLenght

  function shapeFunc(this, point)
    implicit none
    class(Line1D3NodeDT), intent(inout)         :: this
    class(PointDT)      , intent(in)            :: point
    real(rkind)         , dimension(this%nNode) :: shapeFunc
    shapeFunc(1) = 0.5*point%getx()*(point%getx()-1)
    shapeFunc(3) = (1+point%getx())*(1-point%getx())
    shapeFunc(2) = 0.5*point%getx()*(1+point%getx())
  end function shapeFunc

  function dShapeFunc(this, point)
    implicit none
    class(Line1D3NodeDT), intent(inout)                   :: this
    class(PointDT)      , intent(in)                      :: point
    real(rkind)         , dimension(this%dim, this%nNode) :: dShapeFunc
    dShapeFunc(1,1) = point%getx()-0.5
    dShapeFunc(1,3) = -2*point%getx()
    dShapeFunc(1,2) = point%getx()+0.5
  end function dShapeFunc
  
  function jacobianAllNodes(this, pointToValue, node)
    implicit none
    class(Line1D3NodeDT), intent(inout)                     :: this
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(this%nNode), intent(in) :: node
    real(rkind)         , dimension(this%dim, this%dim)     :: jacobianAllNodes
    integer(ikind)                                          :: i
    real(rkind)         , dimension(1, NNODE)               :: dsf
    jacobianAllNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, this%nNode
       jacobianAllNodes(1,1) = jacobianAllNodes(1,1) + dsf(1,i)*node(i)%getx()
    end do
  end function jacobianAllNodes

  function jacobianSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Line1D3NodeDT), intent(inout)                     :: this
    integer(ikind)      , dimension(:)         , intent(in) :: indexList
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(:)         , intent(in) :: node
    real(rkind)         , dimension(this%dim, this%dim)      :: jacobianSomeNodes
    integer(ikind)                                           :: i
    real(rkind)         , dimension(1, NNODE)                :: dsf
    jacobianSomeNodes = 0._rkind
    dsf = this%dShapeFunc(pointToValue)
    do i = 1, size(node)
       jacobianSomeNodes(1,1) = jacobianSomeNodes(1,1) + dsf(1,i)*node(i)%getx()
    end do
  end function jacobianSomeNodes

  function jacobianAtGPoints(this, node)
    implicit none
    class(Line1D3NodeDT)                                , intent(inout) :: this
    class(NodePtrDT), dimension(this%nNode)             , intent(in)    :: node
    real(rkind), dimension(this%integrator%integTerms,this%dim,this%dim) :: jacobianAtGPoints
    integer(ikind)                                                      :: i, j
    real(rkind), dimension(this%integrator%integTerms,1,NNODE)           :: dsf
    jacobianAtGPoints = 0._rkind
    dsf = this%integrator%dShapeFunc
    do i = 1, this%integrator%integTerms
       do j = 1, this%nNode
          jacobianAtGPoints(i,1,1) = jacobianAtGPoints(i,1,1) + dsf(i,1,j)*node(j)%getx()
       end do
    end do
  end function jacobianAtGPoints

  real(rkind) function jacobianDetFromCoordAllNodes(this, pointToValue, node)
    implicit none
    class(Line1D3NodeDT), intent(inout)                     :: this
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(this%nNode), intent(in) :: node
    real(rkind)         , dimension(1,1)                    :: jacobian
    jacobian = this%jacobian(pointToValue, node)
    jacobianDetFromCoordAllNodes = jacobian(1,1)
  end function jacobianDetFromCoordAllNodes

  real(rkind) function jacobianDetFromCoordSomeNodes(this, indexList, pointToValue, node)
    implicit none
    class(Line1D3NodeDT), intent(inout)                     :: this
    integer(ikind)      , dimension(:), intent(in)          :: indexList
    class(PointDT)      , intent(in)                        :: pointToValue
    class(NodePtrDT)    , dimension(:), intent(in)          :: node
    real(rkind)         , dimension(1,1)                    :: jacobian
    jacobian = this%jacobian(indexList, pointToValue, node)
    jacobianDetFromCoordSomeNodes = jacobian(1,1)
  end function jacobianDetFromCoordSomeNodes

  real(rkind) function jacobianDetFromJacobian(this, jacobian)
    implicit none
    class(Line1D3NodeDT)       , intent(inout) :: this
    real(rkind), dimension(:,:), intent(in)    :: jacobian
    jacobianDetFromJacobian = jacobian(1,1)
  end function jacobianDetFromJacobian

  function jacobianDetAtGPointsFromCoord(this, node)
    implicit none
    class(Line1D3NodeDT)                   , intent(inout)      :: this
    class(NodePtrDT), dimension(this%nNode), intent(in)         :: node
    real(rkind)     , dimension(this%integrator%integTerms)     :: jacobianDetAtGPointsFromCoord
    integer(ikind)                                              :: i
    real(rkind)     , dimension(this%integrator%integTerms,1,1) :: jacobian
    jacobian = this%jacobianAtGPoints(node)
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromCoord(i) = jacobian(i,1,1)
    end do
  end function jacobianDetAtGPointsFromCoord

  function jacobianDetAtGPointsFromJacobian(this, jacobian)
    implicit none
    class(Line1D3NodeDT)               , intent(inout) :: this
    real(rkind), dimension(:,:,:)      , intent(in)    :: jacobian
    real(rkind), dimension(this%integrator%integTerms) :: jacobianDetAtGPointsFromJacobian
    integer(ikind)                                     :: i
    do i = 1, this%integrator%integTerms
       jacobianDetAtGPointsFromJacobian(i) = jacobian(i,1,1)
    end do
  end function jacobianDetAtGPointsFromJacobian

  subroutine valueShapeFuncAtGPoints(this)
    implicit none
    class(Line1D3NodeDT), intent(inout) :: this
    integer(ikind)                      :: i
    integer(ikind)                      :: integTerms
    real(rkind)                         :: x
    integTerms = this%integrator%integTerms
    allocate(this%integrator%shapeFunc(integTerms, NNODE))
    allocate(this%integrator%dShapeFunc(integTerms, 1, NNODE))
    do i = 1, integTerms
       x = this%integrator%gPoint(i,1)
       this%integrator%shapeFunc(i, 1:NNODE) = this%shapeFunc(point(x))
       this%integrator%dShapeFunc(i, 1:1, 1:NNODE) = this%dShapeFunc(point(x))
    end do
  end subroutine valueShapeFuncAtGPoints

end module Line1D3NodeM
