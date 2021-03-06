module CFDStrategyM

  use UtilitiesM
  use SparseKit
  use DebuggerM

  use ProcessM
  
  use NewStrategyM
  use SolvingStrategyM

  use PrintM

  use CFDSchemeM
  use CFDBuilderAndSolverM
  use SchemeM
  use BuilderAndSolverM

  use NavierStokes2DM

  use AdamsB4M
  use ExplicitEulerM

  implicit none

  private
  public :: CFDStrategyDT

  type, extends(NewStrategyDT) :: CFDStrategyDT
   contains
     procedure, nopass :: useNewStrategy => CFDStrategy 
  end type CFDStrategyDT

contains

  subroutine CFDStrategy(this)
    implicit none
    class(ProcessDT), intent(inout)          :: this
    type(CFDSchemeDT)                        :: scheme
    type(CFDBuilderAndSolverDT)              :: builAndSolve
    type(PrintDT)                            :: writeOutput
    type(NavierStokes2DDT)                   :: NavierStokes2D
    type(ExplicitEulerDT)                    :: ExplicitEuler
    type(AdamsB4DT)                          :: AdamsBash4 
    real(rkind), dimension(:,:), allocatable :: oldDof 
    real(rkind), dimension(:)  , allocatable :: auxDof, auxRhs
    real(rkind)                              :: dtMin, dtMin1, t
    real(rkind)                              :: factor, error, porc
    real(rkind)                              :: errorTol, error1(4), error2(4)
    integer(ikind)                           :: maxIter, iNode, nNode, i
    integer(ikind)                           :: step1, step2, printStep
    integer(ikind)                           :: flagg, stab, RK
    logical                                  :: multi_step = .true.
    select type(this)
    class is(SolvingStrategyDT)
       call debugLog('  *** Transient Strategy ***')
       print'(A)', '*** Transient Strategy ***'
       nNode = this%application%model%getnNode()
       allocate(this%process, source = WriteOutput)
       allocate(this%application%model%processInfo%mat(4,nNode))
       allocate(this%application%model%processInfo%vect(nNode))
       allocate(oldDof(4,nNode), auxRhs(nNode*4), auxDof(nNode*4))
       call this%application%model%processInfo%initProcess(3)
       call this%application%model%processInfo%setProcess(2,1) 
       this%scheme           = SetScheme(scheme)
       this%builderAndSolver = SetBuilderAndSolver(builAndSolve)
       printStep = this%application%model%processInfo%getPrintStep()
       errorTol  = this%application%model%processInfo%getErrorTol()
       maxIter   = this%application%model%processInfo%getMaxIter()
       error     = errorTol+1
       error1    = errorTol+1
       error2    = 1.d0
       step1     = 0
       step2     = printStep
       t         = 0.d0
       auxDof    = 0.d0
       auxRhs    = 0.d0
       flagg     = 1
       stab      = 1
       RK        = 4
       call this%application%model%processInfo%setProcess(1,1)
       call calculateMass(this%application%model)
       call this%application%model%processInfo%setProcess(1,0)
       call WriteOutput%initPrint()
       call applyDirichlet(this%application%model)
       do while (step1 .lt. maxIter .and. error .gt. errorTol)
          step1 = step1 + 1
          call calculateDt(this%application%model)
          dtMin = this%application%model%processInfo%getDt()
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
          t      = t + dtmin
          oldDof = this%application%model%dof
          call this%application%model%processInfo%setStep(step1)
          if (flagg .le. 5) then
             do i = 1, RK
                factor = (1.d0/(RK+1.d0-i))
                if (i == 1) then
                   call this%application%model%processInfo%setProcess(3,1)
                else
                   call this%application%model%processInfo%setProcess(3,0)
                end if
                this%application%model%rhs = 0.d0
                call builAndSolve%buildAndSolve(this%application%model)
                do iNode = 1, nNode
                   auxDof(iNode*4-3:iNode*4) = this%application%model%dof(:,iNode)
                   auxRhs(iNode*4-3:iNode*4) = this%application%model%rhs(:,iNode)
                end do
                navierStokes2D = SetNavierStokes2D(auxDof, auxRhs, ExplicitEuler, step1)
                call navierStokes2D%integrate(factor*dtMin, multi_step)
                auxDof = navierStokes2D%getState()
                do iNode = 1, nNode
                   this%application%model%dof(:,iNode) = auxDof(iNode*4-3:iNode*4)
                end do
                call applyDirichlet(this%application%model)
             end do
          else
             if (stab == 4) stab = 1
             if (stab == 2) then
                call this%application%model%processInfo%setProcess(3,1)
             else
                call this%application%model%processInfo%setProcess(3,0)
             end if
             stab = stab + 1
             this%application%model%rhs = 0.d0
             call builAndSolve%buildAndSolve(this%application%model)
             do iNode = 1, nNode
                auxDof(iNode*4-3:iNode*4) = this%application%model%dof(:,iNode)
                auxRhs(iNode*4-3:iNode*4) = this%application%model%rhs(:,iNode)
             end do
             navierStokes2D = SetNavierStokes2D(auxDof, auxRhs, AdamsBash4, step1)
             call navierStokes2D%integrate(dtMin, multi_step)
             auxDof = navierStokes2D%getState()
             do iNode = 1, nNode
                this%application%model%dof(:,iNode) = auxDof(iNode*4-3:iNode*4)
             end do
             call applyDirichlet(this%application%model)
          end if
          if (step1 == step2 .or. step1 == maxIter) then
             error1 = 0.d0
             error2 = 0.d0 
             do iNode = 1, nNode
                error1(:)  = error1(:) + (this%application%model%dof(:,iNode)-oldDof(:,iNode))**2
                error2(:)  = error2(:) + oldDof(:,iNode)**2
             end do
             error = maxval(sqrt(error1/error2))
             if (error.gt.1.d2) then
                call debugLog('CONVERGENCE ERROR')
                print'(A)', 'CONVERGENCE ERROR'
                stop
             end if
             call scheme%calculateOutputs(this%application%model)
             call writeOutput%print(step1, this%application%model%results%density           &
                  , this%application%model%results%internalEnergy, this%application%model%results%mach       &
                  , this%application%model%results%pressure      , this%application%model%results%temperature&
                  , this%application%model%results%velocity                                )
             call debugLog('::::::::::::::::::::::::::::::::::::::::')
             call debugLog('Step     : ', step1)
             call debugLog('Error Ec. de Continuidad  = ', sqrt(error1(1)/error2(1)))
             call debugLog('Error Ec. de Momento x    = ', sqrt(error1(2)/error2(2)))
             call debugLog('Error Ec. de Momento y    = ', sqrt(error1(3)/error2(3)))
             call debugLog('Error Ec. de Energia      = ', sqrt(error1(4)/error2(4)))
             call debugLog('t        : ', t)
             call debugLog('dt       : ', dtMin)
             call debugLog('Mach Max = ', maxval(this%application%model%results%mach))
             call debugLog('::::::::::::::::::::::::::::::::::::::::')
             print'(A40)'      , '::::::::::::::::::::::::::::::::::::::::' 
             print'(A11,I10   )', 'Step     : ', step1
             print'(A29,E10.3)', 'Error Ec. de Continuidad  = ', sqrt(error1(1)/error2(1))
             print'(A29,E10.3)', 'Error Ec. de Momento x    = ', sqrt(error1(2)/error2(2))
             print'(A29,E10.3)', 'Error Ec. de Momento y    = ', sqrt(error1(3)/error2(3))
             print'(A29,E10.3)', 'Error Ec. de Energia      = ', sqrt(error1(4)/error2(4))
             print'(A11,E10.3)', 't        : ', t
             print'(A11,E10.3)', 'dt       : ', dtMin
             print'(A11,E10.3)', 'Mach Max = ', maxval(this%application%model%results%mach)
             print'(A40)'      , '::::::::::::::::::::::::::::::::::::::::' 
             step2 = step2 + printStep
          end if
          flagg = flagg + 1
       end do
       call debugLog('*** Finished Integration ***')
       print'(A)', '*** Finished Integration ***'
    class default
       stop 'strategy: unsupported class.'
    end select
  end subroutine CFDStrategy

end module CFDStrategyM
