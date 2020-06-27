module NewSolvingStrategyM

  use NewStrategyM
  use ProcessM
  
  implicit none

  private
  public :: NewSolvingStrategyDT
  
  type, abstract, extends(ProcessDT) :: NewSolvingStrategyDT
     class(NewStrategyDT), allocatable :: strategy
   contains
     procedure :: setStrategy
     procedure :: changeStrategy
     procedure :: useStrategy
  end type NewSolvingStrategyDT

contains

  subroutine setStrategy(this, strategy)
    implicit none
    class(NewSolvingStrategyDT), intent(inout) :: this
    class(NewStrategyDT)       , intent(in)    :: strategy
    if (allocated(this%strategy)) deallocate(this%strategy)
    allocate(this%strategy, source = strategy)
  end subroutine setStrategy

  subroutine useStrategy(NewSolvingStrategy)
    implicit none
    class(NewSolvingStrategyDT) :: NewSolvingStrategy
    if (allocated(NewSolvingStrategy%strategy)) then
       call NewSolvingStrategy%strategy%useNewStrategy(NewSolvingStrategy)
    else
       stop 'strategy: no strategy procedure available.'
    end if
  end subroutine useStrategy

  subroutine changeStrategy(this, newStrategy)
    implicit none
    class(NewSolvingStrategyDT), intent(inout) :: this
    class(NewStrategyDT)       , intent(in)    :: newStrategy
    deallocate(this%strategy)
    allocate(this%strategy, source = newStrategy)
  end subroutine changeStrategy
  
end module NewSolvingStrategyM
