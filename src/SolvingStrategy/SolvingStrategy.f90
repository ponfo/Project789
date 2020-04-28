module SolvingStrategyM
  
  use ProcessM

  use BuilderAndSolverM
  use SchemeM

  implicit none

  private
  public :: NewSolvingStrategyDT, SolvingStrategyDT, SetSolvingStrategy

  type, abstract, extends(ProcessDT) :: NewSolvingStrategyDT
     type(BuilderAndSolverDT)        :: builderAndSolver
     type(SchemeDT)                  :: scheme
   contains
     procedure(Strategy_Procedure), deferred :: useStrategy
  end type NewSolvingStrategyDT

  abstract interface
     subroutine Strategy_Procedure(this)
       import NewSolvingStrategyDT
       class(NewSolvingStrategyDT), intent(inout) :: this
     end subroutine Strategy_Procedure
  end interface

  interface SetSolvingStrategy
     procedure :: constructor
  end interface SetSolvingStrategy
  
  type SolvingStrategyDT
     class(NewSolvingStrategyDT)   , allocatable :: strategy
!!$     class(BuilderAndSolverDT)     , allocatable :: builderAndSolver
!!$     class(SchemeDT)               , allocatable :: scheme
   contains
     procedure :: init
     procedure :: change
     procedure :: sequence
!!$     procedure :: solve    
!!$     procedure :: calculateOutputData 
!!$     procedure :: Isconverged  !convergence check
  end type SolvingStrategyDT

contains

  type(SolvingStrategyDT) function constructor(strategy)
    implicit none
    class(NewSolvingStrategyDT), intent(in) :: strategy
    call constructor%init(strategy)
  end function constructor

  subroutine init(this, strategy)
    implicit none
    class(SolvingStrategyDT)   , intent(inout) :: this
    class(NewSolvingStrategyDT), intent(in)    :: strategy
    allocate(this%strategy, source = strategy)
  end subroutine init
  
  subroutine change(this, newStrategy)
    implicit none
    class(SolvingStrategyDT)   , intent(inout) :: this
    class(NewSolvingStrategyDT), intent(in)    :: newStrategy
    deallocate(this%strategy)
    allocate(this%strategy, source = newStrategy)
  end subroutine change

  subroutine sequence(this)
    implicit none
    class(SolvingStrategyDT), intent(inout) :: this
    call this%strategy%useStrategy()
  end subroutine sequence
  
!!$  subroutine solve(this)
!!$    implicit none
!!$    class(SolvingStrategyDT), intent(inout) :: this
!!$    call this%builderAndSolver%buildAndSolve()
!!$  end subroutine solve
!!$
!!$  subroutine calculateOutputData(this)
!!$    implicit none
!!$    class(SolvingStrategyDT), intent(inout) :: this
!!$    call this%scheme%outputdata()
!!$  end subroutine calculateOutputData
  
end module SolvingStrategyM
