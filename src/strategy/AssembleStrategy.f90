module AssembleStrategyM
  use UtilitiesM
  use DebuggerM

  implicit none

  private
  public :: AssembleStrategyDT

  type, abstract :: AssembleStrategyDT
   contains
     procedure(assembleProcedureInterf), deferred :: assembleProcedure
  end type AssembleStrategyDT

  abstract interface
     subroutine assembleProcedureInterf(this)
       import :: AssembleStrategyDT
       implicit none
       class(AssembleStrategyDT), intent(inout) :: this
     end subroutine assembleProcedureInterf
  end interface

end module AssembleStrategyM
