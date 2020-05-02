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
     integer(ikind)     :: nNode
     integer(ikind)     :: dim
     type(IntegratorDT) :: integrator
   contains
     procedure(shapeFuncInterf)                    , deferred :: shapeFunc
     procedure(dShapeFuncInterf)                   , deferred :: dShapeFunc
     procedure(jacobianAllNodesInterf)             , deferred :: jacobianAllNodes
     procedure(jacobianSomeNodesInterf)            , deferred :: jacobianSomeNodes
     procedure(jacobianDetFromCoordAllNodesInterf) , deferred :: jacobianDetFromCoordAllNodes
     procedure(jacobianDetFromCoordSomeNodesInterf), deferred :: jacobianDetFromCoordSomeNodes
     procedure(jacobianDetFromJacobianInterf)      , deferred :: jacobianDetFromJacobian
     generic                                                  :: jacobian    => jacobianAllNodes     &
                                                                              , jacobianSomeNodes
     generic                                                  :: jacobianDet =>                      &
                                                                   jacobianDetFromCoordAllNodes      &
                                                                 , jacobianDetFromCoordSomeNodes     &
                                                                 , jacobianDetFromJacobian
  end type GeometryDT

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
       
 
end module GeometryM
