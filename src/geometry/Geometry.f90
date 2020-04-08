module GeometryM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use IntegratorM
  use IntegratorPtrM
  
  implicit none

  private
  public :: GeometryDT

  type, abstract :: GeometryDT
     private
     integer(ikind)                             :: nNode
     type(NodePtrDT), dimension(:), allocatable :: node
     type(IntegratorPtrDT)                      :: integrator
  end type GeometryDT

end module GeometryM
