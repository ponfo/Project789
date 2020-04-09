module GeometryM
  use UtilitiesM
  use DebuggerM

  use IntegratorM
  
  implicit none

  private
  public :: GeometryDT

  type, abstract :: GeometryDT
     private
     integer(ikind)     :: nNode
     type(IntegratorDT) :: integrator
  end type GeometryDT

end module GeometryM
