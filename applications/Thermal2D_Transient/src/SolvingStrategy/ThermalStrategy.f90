module ThermalStrategyM
  
  use UtilitiesM
  use SparseKit

  use ThermalmodelM

  use SolvingStrategyM

  use Calculate_dtM
  use PrintM

  use ThermalSchemeM
  use ThermalBuilderAndSolverM
  use SchemeM
  use Poisson2DM
  use RK4M
  use BuilderAndSolverM

  implicit none

  private
  public :: ThermalStrategyDT

  type, extends(NewSolvingStrategyDT) :: ThermalStrategyDT
   contains
     procedure :: buildStrategyAndSolve 
  end type ThermalStrategyDT

contains

  subroutine buildStrategyAndSolve(this, model)
    implicit none
    class(ThermalStrategyDT), intent(inout) :: this
    class(ThermalmodelDT)   , intent(inout) :: model
    type(Sparse)                            :: inverseMatrix
    type(ThermalSchemeDT)                   :: directScheme
    type(ThermalBuilderAndSolverDT)         :: directBAndS
    type(Poisson2DDT)                       :: poisson2D
    type(RK4DT)                             :: rk4
    type(Calculate_dtDT)                    :: calculate_dt
    type(PrintDT)                           :: writeOutput
    real(rkind), dimension(:), allocatable  :: rhs
    real(rkind)                             :: t
    real(rkind)                             :: dt
    real(rkind)                             :: error
    integer(ikind)                          :: step1
    integer(ikind)                          :: step2
    write(*,*) '*** Transient Strategy ***'
    allocate(this%scheme, source = SetScheme(directScheme))
    allocate(this%builderAndSolver, source = SetBuilderAndSolver(directBAndS))
    model%dof    = 0._rkind
    step1       = 0
    step2       = 0
    t           = 0._rkind
    error       = 1._rkind
    call directBAndS%buildAndSolve(model)
    !call DirectScheme%calculateFlux(model)
    allocate(this%process , source = calculate_dt)
    dt    = calculate_dt%calculate(model)*500
    inverseMatrix = inverse(model%mass)
    deallocate(this%process)
    allocate(this%process, source = WriteOutput)
    call WriteOutput%initPrint()
    do while(error .ge. 5e-8)
    poisson2D = SetPoisson2D(model%dof, model%lhs&
         , model%rhs, inverseMatrix, rk4         )
       if (step1 == step2) then
          call writeOutput%print(model%dof, model%heatFlux, step1)
          step2 = step2 + 500
       write(*,*) 't = ', t, 'error = ', error
       end if      
       call poisson2D%integrate(dt)
       model%dof = poisson2D%getState()
       t         = t + dt
       call applyDirichlet(model)
       !call directScheme%calculateFlux(model)
       rhs  = (model%rhs-model%lhs*model%dof)
       error = sqrt(dot_product(rhs, rhs))
       step1 = step1 + 1
    end do
    call directScheme%calculateFlux(model)
    write(*,*) 't final = ', t, 'error = ', error 
    write(*,*) '*** Finished Integration ***'
  end subroutine buildStrategyAndSolve
  
end module ThermalStrategyM
