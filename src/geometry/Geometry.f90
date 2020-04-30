module GeometryM
  use UtilitiesM
  use DebuggerM

  use PointM

  use IntegratorM
  
  implicit none

  private
  public :: GeometryDT

  type, abstract :: GeometryDT
     integer(ikind)     :: nNode
     integer(ikind)     :: dim
     type(IntegratorDT) :: integrator
   contains
     procedure(shapeFuncInterf)              , deferred :: shapeFunc
     procedure(dShapeFuncInterf)             , deferred :: dShapeFunc
     procedure(jacobianInterf)               , deferred :: jacobian
     procedure(jacobianDetFromCoordInterf)   , deferred :: jacobianDetFromCoord
     procedure(jacobianDetFromJacobianInterf), deferred :: jacobianDetFromJacobian
     generic                                            :: jacobianDet => jacobianDetFromCoord  &
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
     function jacobianInterf(this, pointToValue, nodalPoints)
       use UtilitiesM
       import GeometryDT
       import PointDT
       implicit none
       class(GeometryDT)                       , intent(inout) :: this
       class(PointDT)                          , intent(in)    :: pointToValue
       class(PointDT)   , dimension(this%nNode), intent(in)    :: nodalPoints
       real(rkind)      , dimension(this%dim, this%dim)        :: jacobianInterf
     end function jacobianInterf
  end interface

  abstract interface
     function jacobianDetFromCoordInterf(this, pointToValue, nodalPoints)
       use UtilitiesM
       import GeometryDT
       import PointDT
       implicit none
       class(GeometryDT)                       , intent(inout) :: this
       class(PointDT)                          , intent(in)    :: pointToValue
       class(PointDT)   , dimension(this%nNode), intent(in)    :: nodalPoints
       real(rkind)                                             :: jacobianDetFromCoordInterf
     end function jacobianDetFromCoordInterf
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
