module NodePtrM
  use NodeM

  implicit none

  private
  public :: NodePtrDT

  type :: NodePtrDT
     class(NodeDT), pointer :: ptr
  end type NodePtrDT

end module NodePtrM
