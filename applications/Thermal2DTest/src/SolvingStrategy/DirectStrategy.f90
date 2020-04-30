module ThermalStrategyM

  use ThermalmodelM
  
  use SolvingStrategyM

  !use DirectSchemeM
  use DirectBuilderAndSolverM
  use SchemeM
  use BuilderAndSolverM

  implicit none

  private
  public :: ThermalStrategyDT

  type, extends(NewSolvingStrategyDT) :: ThermalStrategyDT
   contains
     procedure :: useProcess  => process
     procedure :: buildStrategyAndSolve 
  end type ThermalStrategyDT

contains

  subroutine process(this)
    implicit none
    class(ThermalStrategyDT), intent(inout) :: this
  end subroutine process

  subroutine buildStrategyAndSolve(this, model)
    implicit none
    class(ThermalStrategyDT), intent(inout) :: this
    class(ThermalmodelDT)   , intent(inout) :: model
    !type(DirectSchemeDT)                    :: directScheme
    type(DirectBuilderAndSolverDT)          :: directBAndS
    !allocate(this%scheme, source = SetScheme(directScheme))
    allocate(this%builderAndSolver, source = SetBuilderAndSolver(directBAndS))
    call directBAndS%buildAndSolve(model)
    !call DirectScheme%calculateFlux(model)
  end subroutine buildStrategyAndSolve
  
end module ThermalStrategyM
