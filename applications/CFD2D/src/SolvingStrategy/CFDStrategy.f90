
module CFDStrategyM

  use UtilitiesM
  use SparseKit
  use DebuggerM

  use CFDmodelM

  use SolvingStrategyM

  use PrintM

  use CFDSchemeM
  use CFDBuilderAndSolverM
  use SchemeM
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
    class(CFDStrategyDT), intent(inout)      :: this
    class(CFDmodelDT)   , intent(inout)      :: model
    type(CFDSchemeDT)                        :: scheme
    type(CFDBuilderAndSolverDT)              :: builAndSolve
    type(PrintDT)                            :: writeOutput
    real(rkind), dimension(:,:), allocatable :: oldDof, newDof
    real(rkind), dimension(:,:), allocatable :: rhs1, rhs2, rhs3
    real(rkind)                              :: dtMin, dtMin1, t
    real(rkind)                              :: factor, error, porc
    real(rkind)                              :: errorTol, error1(4), error2(4)
    integer(ikind)                           :: maxIter, iNode, nNode, i
    integer(ikind)                           :: step1, step2, printStep
    integer(ikind)                           :: flagg, stab, RK
    call debugLog('  *** Transient Strategy ***')
    print'(A)', '*** Transient Strategy ***'
    nNode = model%getnNode()
    allocate(model%processInfo%mat(4,nNode))
    allocate(model%processInfo%vect(nNode))
    allocate(this%scheme, source = SetScheme(scheme))
    allocate(this%builderAndSolver, source = SetBuilderAndSolver(builAndSolve))
    allocate(this%process, source = WriteOutput)
    allocate(model%processInfo%mat2(3,model%getnElement()))
    allocate(model%processInfo%vect2(model%getnElement()))
    allocate(newDof(4,nNode),oldDof(4,nNode),rhs1(4,nNode),rhs2(4,nNode),rhs3(4,nNode))
    call model%processInfo%initProcess(3)
    call model%processInfo%setProcess(2,1)
    printStep = model%processInfo%getPrintStep()
    errorTol  = model%processInfo%getErrorTol()
    maxIter   = model%processInfo%getMaxIter()
    error     = errorTol+1
    error1    = errorTol+1
    error2    = 1.d0
    step1     = 0
    step2     = 1
    t         = 0.d0
    newDof    = 0.d0
    flagg     = 1
    stab      = 1
    RK        = 4
    call initDof(model)
    call model%processInfo%setProcess(1,1)
    call calculateMass(model)
    call model%processInfo%setProcess(1,0)
    call WriteOutput%initPrint()
    call applyDirichlet(model)
    do while (step1 .lt. maxIter .and. error .gt. errorTol)
       step1 = step1 + 1
       call calculateDt(model)
       dtMin = model%processInfo%getDt()
       if (flagg == 1) then
          dtmin1 = dtmin
          flagg  = 2
       end if
       porc = abs((dtMin-dtMin1)/dtMin)
       if (100.d0*porc .le. 1.d0) then
          dtMin = dtMin1
       else
          dtMin1 = dtMin
          flagg  = 2
       end if
       t      = t + model%processInfo%getDt()
       oldDof = model%dof
       if (flagg .le. 4) then
          do i = 1, RK
             factor = (1.d0/(RK+1.d0-i))
             if (i == 1) then
                call model%processInfo%setProcess(3,1)
             else
                call model%processInfo%setProcess(3,0)
             end if
             model%rhs = 0.d0
             call builAndSolve%buildAndSolve(model)
             !$OMP PARALLEL DO PRIVATE(iNode)
             do iNode = 1, nNode
                newDof(:,iNode) = model%dof(:,iNode) + factor*model%rhs(:,iNode)*dtMin
             end do
             !$OMP END PARALLEL DO
             model%dof = newDof
             call applyDirichlet(model)
             rhs3 = rhs2
             rhs2 = rhs1
             rhs1 = model%rhs
          end do
       else
          if (stab == 4) stab = 1
          if (stab == 2) then
             call model%processInfo%setProcess(3,1)
          else
             call model%processInfo%setProcess(3,0)
          end if
          stab = stab + 1
          model%rhs = 0.d0
          call builAndSolve%buildAndSolve(model)
          factor = 24.d0
          !$OMP PARALLEL DO PRIVATE(iNode)
          do iNode = 1, nNode
             newDof(:,iNode) = model%dof(:,iNode)+ (55.d0*model%rhs(:,iNode)&
                  -59.d0*rhs1(:,iNode)+37.d0*rhs2(:,iNode)-9.d0*rhs3(:,iNode))*dtMin/factor
          end do
          !$OMP END PARALLEL DO
          rhs3 = rhs2
          rhs2 = rhs1
          rhs1 = model%rhs
          model%dof = newDof
          call applyDirichlet(model)
       end if
       if (step1 == step2 .or. step1 == maxIter) then
          error1 = 0.d0
          error2 = 0.d0
          !$OMP PARALLEL DO REDUCTION(+:error1, error2) 
          do iNode = 1, nNode
             error1(:)  = error1(:) + (model%dof(:,iNode)-oldDof(:,iNode))**2
             error2(:)  = error2(:) + oldDof(:,iNode)**2
          end do
          !$OMP END PARALLEL DO
          error = maxval(sqrt(error1/error2))
          if (error.gt.1.d2) then
             call debugLog('CONVERGENCE ERROR')
             print'(A)', 'CONVERGENCE ERROR'
             stop
          end if
          call scheme%calculateOutputs(model)
          call writeOutput%print(step1, model%results%density           &
               , model%results%internalEnergy, model%results%mach       &
               , model%results%pressure      , model%results%temperature&
               , model%results%velocity                                )
          call debugLog('::::::::::::::::::::::::::::::::::::::::')
          call debugLog('Step     : ', step1)
          call debugLog('Error Ec. de Continuidad  = ', sqrt(error1(1)/error2(1)))
          call debugLog('Error Ec. de Momento x    = ', sqrt(error1(2)/error2(2)))
          call debugLog('Error Ec. de Momento y    = ', sqrt(error1(3)/error2(3)))
          call debugLog('Error Ec. de Energia      = ', sqrt(error1(4)/error2(4)))
          call debugLog('t        : ', t)
          call debugLog('dt       : ', dtMin)
          call debugLog('Mach Max = ', maxval(model%results%mach))
          call debugLog('::::::::::::::::::::::::::::::::::::::::')
          print'(A40)'      , '::::::::::::::::::::::::::::::::::::::::' 
          print'(A11,I5   )', 'Step     : ', step1
          print'(A29,E10.3)', 'Error Ec. de Continuidad  = ', sqrt(error1(1)/error2(1))
          print'(A29,E10.3)', 'Error Ec. de Momento x    = ', sqrt(error1(2)/error2(2))
          print'(A29,E10.3)', 'Error Ec. de Momento y    = ', sqrt(error1(3)/error2(3))
          print'(A29,E10.3)', 'Error Ec. de Energia      = ', sqrt(error1(4)/error2(4))
          print'(A11,E10.3)', 't        : ', t
          print'(A11,E10.3)', 'dt       : ', dtMin
          print'(A11,E10.3)', 'Mach Max = ', maxval(model%results%mach)
          print'(A40)'      , '::::::::::::::::::::::::::::::::::::::::' 
          step2 = step2 + printStep
       end if
       flagg = flagg + 1
    end do
    call debugLog('*** Finished Integration ***')
    print'(A)', '*** Finished Integration ***'
  end subroutine buildStrategyAndSolve

  subroutine initDof(model)
    implicit none
    class(CFDModelDT), intent(inout) :: model
    real(rkind) :: Vx, Vy, T, P, rho, mach, Cv
    integer(ikind) :: i
    Vx   = model%processInfo%getConstants(7)
    Vy   = model%processInfo%getConstants(8)
    T    = model%processInfo%getConstants(9)
    P    = model%processInfo%getConstants(10)
    rho  = model%processInfo%getConstants(11)
    mach = model%processInfo%getConstants(12)
    Cv   = model%processInfo%getConstants(4)
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, model%getnNode()
       model%dof(1,i) = rho
       model%dof(2,i) = rho*Vx
       model%dof(3,i) = rho*Vy
       model%dof(4,i) = rho*(Cv*T+0.5d0*(Vx**2+Vy**2)) 
    end do
    !$OMP END PARALLEL DO
  end subroutine initDof

end module CFDStrategyM
