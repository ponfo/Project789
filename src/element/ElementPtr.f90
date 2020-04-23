module ElementPtrM
  use ElementM

  implicit none

  private
  public :: ElementPtrDT

  type :: ElementPtrDT
     class(ElementDT), pointer :: ptr
  end type ElementPtrDT

end module ElementPtrM
