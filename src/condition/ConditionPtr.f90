module ConditionPtrM
  use ConditionM

  implicit none

  private
  public :: ConditionPtrDT

  type :: ConditionPtrDT
     class(ConditionDT), pointer :: ptr
  end type ConditionPtrDT

end module ConditionPtrM
