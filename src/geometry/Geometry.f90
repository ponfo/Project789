module GeometryM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use IntegratorM
  
  implicit none

  private
  public :: GeometryDT

  type, abstract :: GeometryDT
     private
     integer(ikind)                             :: nNode
     type(NodePtrDT), dimension(:), allocatable :: node
     type(IntegratorDT)                         :: integrator
   contains
     procedure
