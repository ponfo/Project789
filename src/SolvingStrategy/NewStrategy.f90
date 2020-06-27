module NewStrategyM

  use ProcessM

  implicit none

  private
  public :: NewStrategyDT

  type, abstract :: NewStrategyDT
   contains
     procedure(NewStrategy_Procedure), nopass, deferred :: useNewStrategy
  end type NewStrategyDT

  abstract interface
     subroutine NewStrategy_Procedure(this)
       import ProcessDT
       class(ProcessDT), intent(inout) :: this
     end subroutine NewStrategy_Procedure
  end interface
  
end module NewStrategyM
