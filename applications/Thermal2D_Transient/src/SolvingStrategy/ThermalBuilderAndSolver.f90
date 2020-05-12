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
  public :: ThermalBuilderAndSolverDT, applyDirichlet, applyNewmann

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
!!$    print*, 'SYSTEM'
!!$    print*, 'LHS'
!!$    call model%lhs%printNonZeros()
!!$    print*, 'Mass'
!!$    call model%mass%printNonZeros()
!!$    print*, 'RHS'
!!$    do i = 1, size(model%rhs)
!!$       print'(A,I0,A,E16.8)', 'rhs(', i, ') = ', model%rhs(i)
!!$    end do
  end subroutine buildAndSolve

  subroutine assembleSystem(model)
    implicit none
    class(ThermalModelDT)    , intent(inout) :: model
    integer(ikind)                           :: i, j, iElem, nElem, nNode, row, col
    type(LeftHandSideDT)                     :: localLHS
    real(rkind), dimension(:)  , allocatable :: localRHS
    real(rkind)                              :: adder
    type(ElementPtrDT)                       :: element
    nElem = model%getnElement()
    model%rhs = 0._rkind
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(localLHS, localRHS)
       do i = 1, nNode
          row = element%getNodeID(i)
          adder = 0._rkind
          do j = 1, nNode
             col = element%getNodeID(j)
             call model%LHS%append(val = localLHS%stiffness(i,j)  &
                  , row = row                                    &
                  , col = col                                    )
             adder = adder + localLHS%mass(i,j)
          end do
          model%rhs(row) = model%rhs(row) + localRHS(i)
          call model%mass%append(val = adder                     &
               , row = row                                      &
               , col = row                                      )
       end do
       call localLHS%free()
       deallocate(localRHS)
    end do
    call model%lhs%makeCRS()
    call model%mass%makeCRS()
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
       call condition%calculateLocalSystem(localLHS, localRHS)
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
    integer(ikind)                       :: i, nNode, nodeID
    type(NodePtrDT)                      :: node
    nNode = model%getnNode()
    do i = 1, nNode
       node = model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID)
          model%rhs(nodeID) = node%ptr%dof(1)%fixedVal
       end if
    end do
  end subroutine applyDirichlet

end module ThermalBuilderAndSolverM
