module GeometryM
  use UtilitiesM
  use DebuggerM

  use PointM
  use NodeM
  use NodePtrM

  use IntegratorM
  
  implicit none

  private
  public :: GeometryDT

  type, abstract :: GeometryDT
     integer(ikind)                 :: nNode
     integer(ikind)                 :: dim
     type(IntegratorDT)             :: integrator
     class(GeometryDT), allocatable :: boundaryGeometry
   contains
     procedure                                                   :: getnNode
     procedure                                                   :: getDim
     procedure                                                   :: getIntegrator

     procedure(getLenghtInterf)                       , deferred :: getLenght
     procedure(shapeFuncInterf)                       , deferred :: shapeFunc
     procedure(dShapeFuncInterf)                      , deferred :: dShapeFunc
     procedure(jacobianAllNodesInterf)                , deferred :: jacobianAllNodes
     procedure(jacobianSomeNodesInterf)               , deferred :: jacobianSomeNodes
     procedure(jacobianAtGPointsInterf)               , deferred :: jacobianAtGPoints
     procedure(jacobianDetFromCoordAllNodesInterf)    , deferred :: jacobianDetFromCoordAllNodes
     procedure(jacobianDetFromCoordSomeNodesInterf)   , deferred :: jacobianDetFromCoordSomeNodes
     procedure(jacobianDetFromJacobianInterf)         , deferred :: jacobianDetFromJacobian
     procedure(jacobianDetAtGPointsFromCoordInterf)   , deferred :: jacobianDetAtGPointsFromCoord
     procedure(jacobianDetAtGPointsFromJacobianInterf), deferred :: jacobianDetAtGPointsFromJacobian
     
     generic :: jacobian    => jacobianAllNodes     &
                             , jacobianSomeNodes
     generic :: jacobianDet => jacobianDetFromCoordAllNodes        &
                             , jacobianDetFromCoordSomeNodes       &
                             , jacobianDetFromJacobian
     generic :: jacobianDetAtGPoints => jacobianDetAtGPointsFromCoord     &
                                      , jacobianDetAtGPointsFromJacobian
  end type GeometryDT

  abstract interface
     real(rkind) function getLenghtInterf(this, node)
       use UtilitiesM
       import GeometryDT
       import NodePtrDT
       implicit none
       class(GeometryDT)                       , intent(inout) :: this
       class(NodePtrDT) , dimension(this%nNode), intent(in)    :: node
     end function getLenghtInterf
  end interface

  abstract interface
     function shapeFuncInterf(this, point)
       use UtilitiesM
       import GeometryDT
       import PointDT
       implicit none
       class(GeometryDT), intent(inout)         :: this
       class(PointDT)   , intent(in)            :: point
       real(rkind)      , dimension(this%nNode) :: shapeFuncInterf
     end function shapeFuncInterf
  end interface

  abstract interface
     function dShapeFuncInterf(this, point)
       use UtilitiesM
       import GeometryDT
       import PointDT
       implicit none
       class(GeometryDT), intent(inout)                   :: this
       class(PointDT)   , intent(in)                      :: point
       real(rkind)      , dimension(this%dim, this%nNode) :: dShapeFuncInterf
     end function dShapeFuncInterf
  end interface

  abstract interface
     function jacobianAllNodesInterf(this, pointToValue, node)
       use UtilitiesM
       import GeometryDT
       import PointDT
       import NodePtrDT
       implicit none
       class(GeometryDT)                       , intent(inout) :: this
       class(PointDT)                          , intent(in)    :: pointToValue
       class(NodePtrDT) , dimension(this%nNode), intent(in)    :: node
       real(rkind)      , dimension(this%dim, this%dim)        :: jacobianAllNodesInterf
     end function jacobianAllNodesInterf
  end interface

  abstract interface
     function jacobianSomeNodesInterf(this, indexList, pointToValue, node)
       use UtilitiesM
       import GeometryDT
       import PointDT
       import NodePtrDT
       implicit none
       class(GeometryDT)                       , intent(inout) :: this
       integer(ikind)   , dimension(:)         , intent(in)    :: indexList
       class(PointDT)                          , intent(in)    :: pointToValue
       class(NodePtrDT) , dimension(:)         , intent(in)    :: node
       real(rkind)      , dimension(this%dim, this%dim)        :: jacobianSomeNodesInterf
     end function jacobianSomeNodesInterf
  end interface

  abstract interface
     function jacobianAtGPointsInterf(this, node)
       use UtilitiesM
       import GeometryDT
       import NodePtrDT
       implicit none
       class(GeometryDT)                                   , intent(inout) :: this
       class(NodePtrDT), dimension(this%nNode)             , intent(in)    :: node
       real(rkind), dimension(this%integrator%integTerms,this%dim,this%dim) :: jacobianAtGPointsInterf
     end function jacobianAtGPointsInterf
  end interface
       
  abstract interface
     function jacobianDetFromCoordAllNodesInterf(this, pointToValue, node)
       use UtilitiesM
       import GeometryDT
       import PointDT
       import NodePtrDT
       implicit none
       class(GeometryDT)                       , intent(inout) :: this
       class(PointDT)                          , intent(in)    :: pointToValue
       class(NodePtrDT) , dimension(this%nNode), intent(in)    :: node
       real(rkind)                                             :: jacobianDetFromCoordAllNodesInterf
     end function jacobianDetFromCoordAllNodesInterf
  end interface

  abstract interface
     function jacobianDetFromCoordSomeNodesInterf(this, indexList, pointToValue, node)
       use UtilitiesM
       import GeometryDT
       import PointDT
       import NodePtrDT
       implicit none
       class(GeometryDT)              , intent(inout) :: this
       integer(ikind)   , dimension(:), intent(in)    :: indexList
       class(PointDT)                 , intent(in)    :: pointToValue
       class(NodePtrDT) , dimension(:), intent(in)    :: node
       real(rkind)                                    :: jacobianDetFromCoordSomeNodesInterf
     end function jacobianDetFromCoordSomeNodesInterf
  end interface

  abstract interface
     function jacobianDetFromJacobianInterf(this, jacobian)
       use UtilitiesM
       import GeometryDT
       implicit none
       class(GeometryDT)                , intent(inout) :: this
       real(rkind)      , dimension(:,:), intent(in)    :: jacobian
       real(rkind)                                      :: jacobianDetFromJacobianInterf
     end function jacobianDetFromJacobianInterf
  end interface

  abstract interface
     function jacobianDetAtGPointsFromCoordInterf(this, node)
       use UtilitiesM
       import GeometryDT
       import NodePtrDT
       implicit none
       class(GeometryDT)                      , intent(inout) :: this
       class(NodePtrDT), dimension(this%nNode), intent(in)    :: node
       real(rkind) , dimension(this%integrator%integTerms)    :: jacobianDetAtGPointsFromCoordInterf
     end function jacobianDetAtGPointsFromCoordInterf
  end interface

  abstract interface
     function jacobianDetAtGPointsFromJacobianInterf(this, jacobian)
       use UtilitiesM
       import GeometryDT
       implicit none
       class(GeometryDT)                 , intent(inout) :: this
       real(rkind), dimension(:,:,:)     , intent(in)    :: jacobian
       real(rkind), dimension(this%integrator%integTerms) :: jacobianDetAtGPointsFromJacobianInterf
     end function jacobianDetAtGPointsFromJacobianInterf
  end interface

contains

  integer(ikind) pure function getnNode(this)
    implicit none
    class(GeometryDT), intent(in) :: this
    getnNode = this%nNode
  end function getnNode

  integer(ikind) pure function getDim(this)
    implicit none
    class(GeometryDT), intent(in) :: this
    getDim = this%dim
  end function getDim

  type(IntegratorDT) pure function getIntegrator(this)
    implicit none
    class(GeometryDT), intent(in) :: this
    getIntegrator = this%integrator
  end function getIntegrator

end module GeometryM
