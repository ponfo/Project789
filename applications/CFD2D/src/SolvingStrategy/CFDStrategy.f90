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
    print*, 'Build and Solve'
    call builAndSolve%buildAndSolve(model)
    print*, 'Scheme'
    call scheme%calculateOutputs(model)
    !::::::::::::::::::::::::::::::::::::::::::::::
    dt    = model%processInfo%getDt()
    print*, 'inverse'
    !call model%mass%printNonZeros()
    inverseMatrix = inverseLumped(model%mass)
    allocate(this%process, source = WriteOutput)
    print*, 'Init output'
    call WriteOutput%initPrint()
    print*, 'Init iterations'
    call model%processInfo%setStep(step1)
    !::::::::::::::::::::::::::::::::::::::::::::::
    do while(error .ge. errorTol)
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
       step1 = step1 + 1
       call model%processInfo%setStep(step1)
       call builAndSolve%buildAndSolve(model)
       call scheme%calculateOutputs(model)
       rhs  = (model%rhs-model%lhs*model%dof)
       error = sqrt(dot_product(rhs, rhs))
    end do
    !:::::::::::::::::::::::::::::::::::::::::::::
    call scheme%calculateOutputs(model)
    write(*,*) 't final = ', t, 'error = ', error 
    write(*,*) '*** Finished Integration ***'
  end subroutine buildStrategyAndSolve

  type(Sparse) function inverseLumped(matrix)
    implicit none
    class(Sparse), intent(inout) :: matrix
    logical                      :: sortRows = .false.
    integer(ikind) :: i
    inverseLumped = Sparse(nnz = matrix%getn(), rows = matrix%getn())
    do i = 1, matrix%getn()
       call inverseLumped%append(1._rkind/matrix%get(i,i),i,i)
    end do
    call inverseLumped%makeCRS(sortRows)
  end function inverseLumped
    
end module CFDStrategyM
