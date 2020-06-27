module ThermalStrategyM
  
  use NewStrategyM
  use SolvingStrategyM
  use ProcessM

  use ThermalSchemeM
  use ThermalBuilderAndSolverM
  use SchemeM
  use BuilderAndSolverM

  implicit none

  private
  public :: ThermalStrategyDT

  type, extends(NewStrategyDT) :: ThermalStrategyDT
   contains
     procedure, nopass :: useNewStrategy => ThermalStrategy 
  end type ThermalStrategyDT

contains

  subroutine ThermalStrategy(this)
    implicit none
    class(ProcessDT), intent(inout) :: this
    type(ThermalSchemeDT)           :: directScheme
    type(ThermalBuilderAndSolverDT) :: directBAndS
    select type(this)
    class is(SolvingStrategyDT)
    this%scheme           = SetScheme(directScheme)
    this%builderAndSolver = SetBuilderAndSolver(directBAndS)
    call directBAndS%buildAndSolve(this%thermalModel)
    call DirectScheme%calculateFlux(this%thermalModel)
    class default
       stop 'strategy: unsupported class.'
    end select
  end subroutine ThermalStrategy
  
end module ThermalStrategyM
