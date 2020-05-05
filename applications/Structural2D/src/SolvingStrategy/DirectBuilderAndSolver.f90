module DirectBuilderAndSolverM

  use UtilitiesM
  use DebuggerM

  use SparseKit

  use NodePtrM
  use ElementPtrM
  use ConditionPtrM

  use StructuralModelM

  use BuilderAndSolverM

  use LinearSolverM
  use mklPardisoM

  implicit none

  private
  public :: DirectBuilderAndSolverDT

  type, extends(NewBuilderAndSolverDT) :: DirectBuilderAndSolverDT
   contains
     procedure :: buildAndSolve
  end type DirectBuilderAndSolverDT

contains

  subroutine buildAndSolve(this, model)
    implicit none
    class(DirectBuilderAndSolverDT), intent(inout) :: this
    class(StructuralModelDT)       , intent(inout) :: model
    integer(ikind) :: i
    write(*,*) '*** Direct Builder And Solver ***'
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call assembleSystem(model)
    call applyBC(model)
    print*, 'SYSTEM'
    print*, 'LHS'
    call model%lhs%printNonZeros()
    print*, 'RHS'
    do i = 1, size(model%rhs)
       print'(A,I0,A,E16.8)', 'rhs(', i, ') = ', model%rhs(i)
    end do
    write(*,*) '*** Init Linear Solver ***'
    call solve(model)
  end subroutine buildAndSolve

  subroutine assembleSystem(model)
    implicit none
    class(StructuralModelDT)    , intent(inout) :: model
    integer(ikind)                              :: iNode, jNode, iDof, jDof, iNodeID, jNodeID
    integer(ikind)                              :: iElem, nElem, nNode, nDof
    real(rkind), dimension(:,:), allocatable    :: localLHS
    real(rkind), dimension(:)  , allocatable    :: localRHS
    type(ElementPtrDT)                          :: element
    nElem = model%getnElement()
    nDof = 2
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(localLHS, localRHS)
       do iNode = 1, nNode
          iNodeID = element%getNodeID(iNode)
          do jNode = 1, nNode
             jNodeID = element%getNodeID(jNode)
             do iDof = 1, nDof
                do jDof = 1, nDof
                   call model%LHS%append(                                                &
                          val = localLHS(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof))  &
                        , row = iNodeID*nDof-(nDof-iDof)                                 &
                        , col = jNodeID*nDof-(nDof-jDof)                                 )
                end do
             end do
          end do
          model%rhs(iNodeID*nDof-1) = localRHS(iNode*nDof-1)
          model%rhs(iNodeID*nDof)   = localRHS(iNode*nDof)
       end do
       deallocate(localLHS)
       deallocate(localRHS)
    end do
    call model%lhs%makeCRS()
  end subroutine assembleSystem

  subroutine applyBC(model)
    implicit none
    class(StructuralModelDT), intent(inout) :: model
    call applyDirichlet(model)
  end subroutine applyBC

  subroutine applyDirichlet(model)
    implicit none
    class(StructuralModelDT), intent(inout) :: model
    integer(ikind)                          :: i, nNode, nodeID, nDof
    type(NodePtrDT)                         :: node
    nNode = model%getnNode()
    nDof = 2
    do i = 1, nNode
       node = model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-1)
          model%rhs(nodeID*nDof-1) = node%ptr%dof(1)%fixedVal
       end if
       if(node%ptr%dof(2)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof)
          model%rhs(nodeID*nDof) = node%ptr%dof(2)%fixedVal
       end if
    end do
  end subroutine applyDirichlet

  subroutine solve(model)
    implicit none
    class(StructuralmodelDT), intent(inout) :: model
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

end module DirectBuilderAndSolverM
