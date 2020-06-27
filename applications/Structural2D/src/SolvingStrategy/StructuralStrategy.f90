module StructuralStrategyM
  
  use UtilitiesM

  use StructuralmodelM
  use ProcessM
  
  use NewStrategyM
  use SolvingStrategyM

  use StructuralSchemeM
  use StructuralBuilderAndSolverM
  use SchemeM
  use BuilderAndSolverM

  implicit none

  private
  public :: StructuralStrategyDT

  type, extends(NewStrategyDT) :: StructuralStrategyDT
   contains
     procedure, nopass :: useNewStrategy => StructuralStrategy  
  end type StructuralStrategyDT

contains

  subroutine StructuralStrategy(this)
    implicit none
    class(ProcessDT)        , intent(inout) :: this
    type(StructuralSchemeDT)                :: directScheme
    type(StructuralBuilderAndSolverDT)      :: directBAndS
    select type(this)
    class is(SolvingStrategyDT)
    this%scheme           = SetScheme(directScheme)
    this%builderAndSolver = SetBuilderAndSolver(directBAndS)
    call directBAndS%buildAndSolve(this%structuralModel)
    call directScheme%calculatePost(this%structuralModel)
    class default
       stop 'strategy: unsupported class.'
    end select
  end subroutine StructuralStrategy
  
end module StructuralStrategyM
