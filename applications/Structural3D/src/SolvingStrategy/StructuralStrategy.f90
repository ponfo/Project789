module StructuralStrategyM
  use UtilitiesM

  use StructuralmodelM
  
  use SolvingStrategyM

  use StructuralSchemeM
  use StructuralBuilderAndSolverM
  use SchemeM
  use BuilderAndSolverM

  implicit none

  private
  public :: StructuralStrategyDT

  type, extends(NewSolvingStrategyDT) :: StructuralStrategyDT
   contains
     procedure :: useProcess  => process
     procedure :: buildStrategyAndSolve 
  end type StructuralStrategyDT

contains

  subroutine process(this)
    implicit none
    class(StructuralStrategyDT), intent(inout) :: this
  end subroutine process

  subroutine buildStrategyAndSolve(this, model)
    implicit none
    class(StructuralStrategyDT), intent(inout) :: this
    class(StructuralmodelDT)   , intent(inout) :: model
    type(StructuralSchemeDT)                   :: directScheme
    type(StructuralBuilderAndSolverDT)         :: directBAndS
    integer(ikind) :: i
    allocate(this%scheme, source = SetScheme(directScheme))
    allocate(this%builderAndSolver, source = SetBuilderAndSolver(directBAndS))
    call directBAndS%buildAndSolve(model)
    call directScheme%calculatePost(model)
  end subroutine buildStrategyAndSolve
  
end module StructuralStrategyM
