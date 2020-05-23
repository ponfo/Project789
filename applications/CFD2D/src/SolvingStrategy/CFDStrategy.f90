module CFDStrategyM
  
  use UtilitiesM
  use SparseKit

  use CFDmodelM
  
  use SolvingStrategyM

  use PrintM

  use CFDSchemeM
  use CFDBuilderAndSolverM
  use SchemeM
  use NavierStokes2DM
  use RK4M
  use BuilderAndSolverM

  implicit none

  private
  public :: CFDStrategyDT

  type, extends(NewSolvingStrategyDT) :: CFDStrategyDT
   contains
     procedure :: buildStrategyAndSolve 
  end type CFDStrategyDT

contains

  subroutine buildStrategyAndSolve(this, model)
    implicit none
    class(CFDStrategyDT), intent(inout)    :: this
    class(CFDmodelDT)   , intent(inout)    :: model
    type(CFDSchemeDT)                      :: scheme
    type(CFDBuilderAndSolverDT)            :: builAndSolve
    type(Sparse)                           :: inverseMatrix
    type(NavierStokes2DDT)                 :: navierStokes2D
    type(RK4DT)                            :: rk4
    type(PrintDT)                          :: writeOutput
    real(rkind), dimension(:), allocatable :: rhs
    real(rkind)                            :: t
    real(rkind)                            :: dt
    real(rkind)                            :: error
    real(rkind)                            :: errorTol
    integer(ikind)                         :: step1
    integer(ikind)                         :: step2
    integer(ikind)                         :: printStep
    write(*,*) '*** Transient Strategy ***'
    allocate(this%scheme, source = SetScheme(scheme))
    allocate(this%builderAndSolver, source = SetBuilderAndSolver(builAndSolve))
    step1       = 0
    step2       = 0
    t           = model%processInfo%getT0()
    errorTol    = model%processInfo%getErrorTol()
    error       = errorTol+1
    printStep   = model%processInfo%getPrintStep()
    call model%processInfo%setStep(step1)
    call builAndSolve%buildAndSolve(model)
    call scheme%calculateOutputs(model)
    !::::::::::::::::::::::::::::::::::::::::::::::
    dt    = model%processInfo%getDt()
    inverseMatrix = inverse(model%mass)
    allocate(this%process, source = WriteOutput)
    call WriteOutput%initPrint()
    !::::::::::::::::::::::::::::::::::::::::::::::
    do while(error .ge. errorTol)
       call model%processInfo%setStep(step1)
       navierStokes2D = SetNavierStokes2D(model%dof, model%lhs&
            , model%rhs, inverseMatrix, rk4                   )
       if (step1 == step2) then
          call writeOutput%print(step1, model%results%density          &
               , model%results%internalEnergy, model%results%mach       &
               , model%results%pressure      , model%results%temperature&
               , model%results%velocity                                )
          step2 = step2 + printStep
          write(*,*) 't = ', t, 'error = ', error
       end if
       call navierStokes2D%integrate(dt)
       model%dof = navierStokes2D%getState()
       t         = t + dt
       call builAndSolve%buildAndSolve(model)
       call scheme%calculateOutputs(model)
       rhs  = (model%rhs-model%lhs*model%dof)
       error = sqrt(dot_product(rhs, rhs))
       step1 = step1 + 1
    end do
    !:::::::::::::::::::::::::::::::::::::::::::::
    call scheme%calculateOutputs(model)
    write(*,*) 't final = ', t, 'error = ', error 
    write(*,*) '*** Finished Integration ***'
  end subroutine buildStrategyAndSolve
  
end module CFDStrategyM
