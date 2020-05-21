module CFDBuilderAndSolverM

  use UtilitiesM
  use DebuggerM

  use SparseKit

  use NodePtrM
  use ElementPtrM
  use ConditionPtrM

  use LeftHandSideM

  use CFDModelM

  use BuilderAndSolverM

  implicit none

  private
  public :: CFDBuilderAndSolverDT

  type, extends(NewBuilderAndSolverDT) :: CFDBuilderAndSolverDT
   contains
     procedure :: buildAndSolve
  end type CFDBuilderAndSolverDT

contains

  subroutine buildAndSolve(this, model)
    implicit none
    class(CFDBuilderAndSolverDT), intent(inout) :: this
    class(CFDModelDT)           , intent(inout) :: model
    integer(ikind) :: i
    call applyBC(model)
    call assembleSystem(model)
!!$    print*, 'SYSTEM'
!!$    print*, 'LHS'
!!$    call model%lhs%printNonZeros()
!!$    print*, 'MASS'
!!$    call model%mass%printNonZeros()
!!$    print*, 'RHS'
!!$    do i = 1, size(model%rhs)
!!$       print'(A,I0,A,E16.8)', 'rhs(', i, ') = ', model%rhs(i)
!!$    end do
  end subroutine buildAndSolve

  subroutine assembleSystem(model)
    implicit none
    class(CFDModelDT)        , intent(inout) :: model
    integer(ikind)                           :: iNode, jNode, iDof, jDof, iNodeID, jNodeID
    integer(ikind)                           :: iElem, nElem, nNode, nDof
    type(LeftHandSideDT)                     :: localLHS
    real(rkind), dimension(:)  , allocatable :: localRHS
    real(rkind)                              :: adder
    type(ElementPtrDT)                       :: element
    nElem = model%getnElement()
    nDof = 4
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(localLHS, localRHS)
       do iNode = 1, nNode
          iNodeID = element%getNodeID(iNode)
          do jNode = 1, nNode
             jNodeID = element%getNodeID(jNode)
             do iDof = 1, nDof
                adder = 0._rkind
                do jDof = 1, nDof
                        adder = adder                                                 &
                        + localLHS%mass(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof))  
                end do
                call model%mass%append(                                                &
                        val = adder                                                   &
                        , row = iNodeID*nDof-(nDof-iDof)                              &
                        , col = iNodeID*nDof-(nDof-iDof)                              )
             end do
          end do
          model%rhs(iNodeID*nDof-3) = model%rhs(iNodeID*nDof-3) + localRHS(iNode*nDof-3)
          model%rhs(iNodeID*nDof-2) = model%rhs(iNodeID*nDof-2) + localRHS(iNode*nDof-2)
          model%rhs(iNodeID*nDof-1) = model%rhs(iNodeID*nDof-1) + localRHS(iNode*nDof-1)
          model%rhs(iNodeID*nDof)   = model%rhs(iNodeID*nDof)   + localRHS(iNode*nDof)
       end do
       call localLHS%free()
       deallocate(localRHS)
    end do
    call model%lhs%makeCRS()
  end subroutine assembleSystem

  subroutine applyBC(model)
    implicit none
    class(CFDModelDT), intent(inout) :: model
    call applyNewmann(model)
    call applyDirichlet(model)
  end subroutine applyBC

  subroutine applyNewmann(model)
    implicit none
    class(CFDModelDT)                   , intent(inout) :: model
    integer(ikind)                                      :: iCond, nCond, iNode, nNode, iNodeID
    real(rkind)             , dimension(:), allocatable :: localRHS
    type(ConditionPtrDT)                                :: condition
    nCond = model%getnCondition()
    do iCond = 1, nCond
       condition = model%getCondition(iCond)
       nNode = condition%getnNode()
       call condition%calculateRHS(localRHS)
       do iNode = 1, nNode
          iNodeID = condition%getNodeID(iNode)
          model%rhs(iNodeID*4-3) = model%rhs(iNodeID*4-3) + localRHS(iNode*4-3)
          model%rhs(iNodeID*4-2) = model%rhs(iNodeID*4-2) + localRHS(iNode*4-2)
          model%rhs(iNodeID*4-1) = model%rhs(iNodeID*4-1) + localRHS(iNode*4-1)
          model%rhs(iNodeID*4  ) = model%rhs(iNodeID*4  ) + localRHS(iNode*4  )
       end do
       deallocate(localRHS)
    end do
  end subroutine applyNewmann

  subroutine applyDirichlet(model)
    implicit none
    class(CFDModelDT), intent(inout) :: model
    integer(ikind)                   :: i, nNode, nodeID, nDof
    type(NodePtrDT)                  :: node
    nNode = model%getnNode()
    nDof = 4
    do i = 1, nNode
       node = model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-3)
          model%rhs(nodeID*nDof-3) = node%ptr%dof(1)%fixedVal
       end if
       if(node%ptr%dof(2)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-2)
          model%rhs(nodeID*nDof-2) = node%ptr%dof(2)%fixedVal
       end if
       if(node%ptr%dof(3)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-1)
          model%rhs(nodeID*nDof-1) = node%ptr%dof(3)%fixedVal
       end if
       if(node%ptr%dof(4)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof)
          model%rhs(nodeID*nDof  ) = node%ptr%dof(4)%fixedVal
       end if
    end do
  end subroutine applyDirichlet
  
end module CFDBuilderAndSolverM
