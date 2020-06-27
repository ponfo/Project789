module ThermalStrategyM

  use UtilitiesM
  use SparseKit

  use NewStrategyM
  use SolvingStrategyM

  use ProcessM
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

  type, extends(NewStrategyDT) :: ThermalStrategyDT
   contains
     procedure, nopass :: useNewStrategy => ThermalStrategy
  end type ThermalStrategyDT

contains

  subroutine ThermalStrategy(this)
    implicit none
    class(ProcessDT), intent(inout)         :: this
    type(Sparse)                            :: inverseMatrix
    type(ThermalSchemeDT)                   :: directScheme
    type(ThermalBuilderAndSolverDT)         :: directBAndS
    type(Poisson2DDT)                       :: poisson2D
    type(RK4DT)                             :: rk4
    type(PrintDT)                           :: writeOutput
    real(rkind), dimension(:), allocatable  :: rhs
    real(rkind)                             :: t
    real(rkind)                             :: dt
    real(rkind)                             :: error
    real(rkind)                             :: errorTol
    integer(ikind)                          :: step1
    integer(ikind)                          :: step2
    integer(ikind)                          :: printStep
    logical                                 :: multi_step = .false.
    select type(this)
    class is(SolvingStrategyDT)
       write(*,*) '*** Transient Strategy ***'
       this%scheme           = SetScheme(directScheme)
       this%builderAndSolver = SetBuilderAndSolver(directBAndS)
       step1       = 0
       step2       = 0
       t           = this%application%model%processInfo%getT0()
       errorTol    = this%application%model%processInfo%getErrorTol()
       error       = errorTol+1
       printStep   = this%application%model%processInfo%getPrintStep()
       call directBAndS%buildAndSolve(this%application%model)
       call DirectScheme%calculateFlux(this%application%model)
       dt    = this%application%model%processInfo%getDt()*500
       inverseMatrix = inverse(this%application%model%mass)
       allocate(this%process, source = WriteOutput)
       call WriteOutput%initPrint()
       do while(error .ge. errorTol)
          call this%application%model%processInfo%setStep(step1)
          poisson2D = SetPoisson2D(this%application%model%dof, this%application%model%lhs&
               , this%application%model%rhs, inverseMatrix, rk4         )
          if (step1 == step2) then
             call writeOutput%print(this%application%model%dof, this%application%model%heatFlux, step1)
             step2 = step2 + printStep
             write(*,*) 't = ', t, 'error = ', error
          end if
          call poisson2D%integrate(dt, multi_step)
          this%application%model%dof = poisson2D%getState()
          t         = t + dt
          call applyDirichlet(this%application%model)
          !call directScheme%calculateFlux(this%application%model)
          rhs  = (this%application%model%rhs-this%application%model%lhs*this%application%model%dof)
          error = sqrt(dot_product(rhs, rhs))
          step1 = step1 + 1
       end do
       call directScheme%calculateFlux(this%application%model)
       write(*,*) 't final = ', t, 'error = ', error 
       write(*,*) '*** Finished Integration ***'
    class default
       stop 'strategy: unsupported class.'
    end select
  end subroutine ThermalStrategy

end module ThermalStrategyM
