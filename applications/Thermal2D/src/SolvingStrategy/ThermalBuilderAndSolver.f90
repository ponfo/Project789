module ThermalBuilderAndSolverM

  use UtilitiesM
  use DebuggerM

  use SparseKit

  use NodePtrM
  use ElementPtrM
  use ConditionPtrM

  use LeftHandSideM

  use ThermalModelM

  use BuilderAndSolverM

  use LinearSolverM
  use mklPardisoM

  implicit none

  private
  public :: ThermalBuilderAndSolverDT

  type, extends(NewBuilderAndSolverDT) :: ThermalBuilderAndSolverDT
   contains
     procedure :: buildAndSolve
  end type ThermalBuilderAndSolverDT

contains

  subroutine buildAndSolve(this, model)
    implicit none
    class(ThermalBuilderAndSolverDT), intent(inout) :: this
    class(ThermalModelDT)          , intent(inout) :: model
    integer(ikind) :: i
    write(*,*) '*** Direct Builder And Solver ***'
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call assembleSystem(model)
    call applyBC(model)
    write(*,*) '*** Init Linear Solver ***'
    call solve(model)
  end subroutine buildAndSolve

  subroutine assembleSystem(model)
    implicit none
    class(ThermalModelDT)    , intent(inout) :: model
    integer(ikind)                           :: i, j, iElem, nElem, nNode, row, col
    type(LeftHandSideDT)                     :: localLHS
    real(rkind), dimension(:)  , allocatable :: localRHS
    type(ElementPtrDT)                       :: element
    nElem = model%getnElement()
    model%rhs = 0._rkind
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(model%processInfo, localLHS, localRHS)
       do i = 1, nNode
          row = element%getNodeID(i)
          do j = 1, nNode
             col = element%getNodeID(j)
             call model%LHS%append(val = localLHS%stiffness(i,j)  &
                  , row = row                                    &
                  , col = col                                    )
          end do
          model%rhs(row) = model%rhs(row) + localRHS(i)
       end do
       call localLHS%free()
       deallocate(localRHS)
    end do
    call model%lhs%makeCRS()
  end subroutine assembleSystem

  subroutine applyBC(model)
    implicit none
    class(ThermalModelDT), intent(inout) :: model
    call applyNewmann(model)
    call applyDirichlet(model)
  end subroutine applyBC

  subroutine applyNewmann(model)
    implicit none
    class(ThermalModelDT)      , intent(inout) :: model
    integer(ikind)                             :: i, j, iCond, nCond, nNode, row, col
    type(LeftHandSideDT)                       :: localLHS
    real(rkind), dimension(:)  , allocatable   :: localRHS
    type(ConditionPtrDT)                       :: condition
    nCond = model%getnCondition()
    do iCond = 1, nCond
       condition = model%getCondition(iCond)
       nNode = condition%getnNode()
       call condition%calculateLocalSystem(model%processInfo, localLHS, localRHS)
       if(condition%getAffectsLHS() == .true.) then
          do i = 1, nNode
             row = condition%getNodeID(i)
             do j = 1, nNode
                col = condition%getNodeID(j)
                call model%LHS%appendPostCRS(val = localLHS%stiffness(i,j)  &
                     , row = row                                           &
                     , col = col                                           )
             end do
          end do
       call localLHS%free()
       end if
       if(condition%getAffectsRHS() == .true.) then
          do i = 1, nNode
             row = condition%getNodeID(i)
             model%rhs(row) = model%rhs(row) + localRHS(i)
          end do
          deallocate(localRHS)
       end if
    end do
  end subroutine applyNewmann

  subroutine applyDirichlet(model)
    implicit none
    class(ThermalModelDT), intent(inout) :: model
    integer(ikind)                       :: i, j, nNode, nodeID
    type(NodePtrDT)                      :: node
    nNode = model%getnNode()
    do i = 1, nNode
       node = model%getNode(i)
       do j = 1, node%getnDof()
          if(node%ptr%dof(j)%getName() == 'TEMPERATURE') then
             if(node%ptr%dof(j)%isFixed) then
                nodeID = node%ptr%getID()
                call model%lhs%setDirichlet(nodeID)
                model%rhs(nodeID) = node%ptr%dof(j)%fixedVal
             end if
          end if
       end do
    end do
  end subroutine applyDirichlet

  subroutine solve(model)
    implicit none
    class(ThermalmodelDT), intent(inout) :: model
    class(LinearSolverDT), allocatable   :: UserSolver
    type(MKLpardisoDT)                   :: MKLPardiso
    real(rkind)                          :: start
    real(rkind)                          :: finish
    real(rkind)                          :: pt(64)
    integer(ikind)                       :: maxfct 
    integer(ikind)                       :: mnum 
    integer(ikind)                       :: mtype 
    integer(ikind)                       :: phase
    integer(ikind)                       :: idum(1)
    integer(ikind)                       :: nrhs 
    integer(ikind)                       :: iparm(64)
    integer(ikind)                       :: msglvl
    integer(ikind)                       :: error
    integer(ikind), dimension(137)       :: data
    integer(ikind)                       :: i
    call debuglog('solving linear system')
    call cpu_time(start)
    allocate(UserSolver, source = SetLinearSolver(MKLPardiso))
    pt            = 0
    maxfct        = 1
    mnum          = 1
    mtype         = 1      ! real and structurally symmetric 
    phase         = 13     ! analisys, numerical factorization, solve,
    ! iterative refinement
    idum          = 0
    nrhs          = 1
    iparm         = 0
    iparm(1)      = 1      ! user defined iparms
    iparm(2)      = 2      ! 3 The parallel (OpenMP) version of the
    !nested dissection algorithm.
    ! 2 The nested dissection algorithm from
    !the METIS package.
    ! 0 The minimum degree algorithm.
    iparm(4)      = 61     ! LU-preconditioned CGS iteration with a
    ! stopping criterion of
    ! 1.0E-6 for nonsymmetric matrices.
    iparm(24)     = 1      ! two-level factorization algorithm.
    iparm(60)     = 1      ! in-core mode or out-core mode
    msglvl        = 0      ! non-print statistical information.
    error         = 0     
    data(1:64)    = pt    
    data(65)      = maxfct
    data(66)      = mnum
    data(67)      = mtype
    data(68)      = phase
    data(69)      = idum(1)  
    data(70)      = nrhs
    do i = 1, 64
       data(70+i) = iparm(i)
    end do
    data(136)     = msglvl
    data(137)     = error

    call UserSolver%solve(model%rhs, model%lhs, model%dof, data)

    call cpu_time(finish)
    print'(a,e14.7)', 'solver time => ', (finish-start)
    call debuglog('done solving')
  end subroutine solve

end module ThermalBuilderAndSolverM
