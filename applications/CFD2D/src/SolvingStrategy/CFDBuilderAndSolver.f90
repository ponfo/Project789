module CFDBuilderAndSolverM

  use UtilitiesM
  use SparseKit
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
    call applyNewmann(model)
    call assembleSystem(model)
    call applyDirichlet(model)
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
    call model%lhs%free()
    model%lhs = Sparse(nnz = model%getnElement()*256, rows = model%getnNode()*4)
    nElem = model%getnElement()
    nDof = 4
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(model%processInfo,localLHS, localRHS)
       do iNode = 1, nNode
          iNodeID = element%getNodeID(iNode)
          do jNode = 1, nNode
             jNodeID = element%getNodeID(jNode)
             do iDof = 1, nDof
                adder = 0._rkind
                do jDof = 1, nDof
                   call model%LHS%append(                                                        &
                        val = localLHS%stiffness(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof)) &
                        , row = iNodeID*nDof-(nDof-iDof)                                        &
                        , col = jNodeID*nDof-(nDof-jDof)                                        )
                   adder = adder                                                                &
                        + localLHS%mass(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof))  
                end do
                if (model%processInfo%getStep() .eq. 0) then
                   call model%mass%append(                                                       &
                        val = adder                                                             &
                        , row = iNodeID*nDof-(nDof-iDof)                                        &
                        , col = iNodeID*nDof-(nDof-iDof)                                        )
                end if
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
    if (model%processInfo%getStep() .eq. 0) then
       call model%mass%makeCRS()
    end if
    call model%lhs%makeCRS()
  end subroutine assembleSystem

  subroutine applyNewmann(model)
    implicit none
    class(CFDModelDT)                   , intent(inout) :: model
    integer(ikind)                                      :: iCond, nCond, iNode, nNode, iNodeID
    integer(ikind)          , dimension(:), allocatable :: position 
    real(rkind)             , dimension(:), allocatable :: localRHS
    type(ConditionPtrDT)                                :: condition
    nCond = model%getnCondition()
    allocate(position(model%getnNode()))
    position = 0
    do iCond = 1, nCond
       condition = model%getCondition(iCond)
       nNode = condition%getnNode()
       call condition%calculateRHS(model%processInfo, localRHS)
       do iNode = 1, nNode
          iNodeID = condition%getNodeID(iNode)
          position(iNodeID) = position(iNodeID) + 1
          if (position(iNodeID) == 2) then
             model%dof(iNodeID*4-2) = (model%dof(iNodeID*4-2)&
                  +(model%dof(iNodeID*4-2)*localRHS(iNode*2-1)+model%dof(iNodeID*4-1)*localRHS(iNode*2))*localRHS(iNode*2-1))/2.d0
             model%dof(iNodeID*4-1) = (model%dof(iNodeID*4-1)&
                  +(model%dof(iNodeID*4-2)*localRHS(iNode*2-1)+model%dof(iNodeID*4-1)*localRHS(iNode*2))*localRHS(iNode*2))/2.d0
          else
             model%dof(iNodeID*4-2) = (model%dof(iNodeID*4-2)*localRHS(iNode*2-1)+model%dof(iNodeID*4-1)*localRHS(iNode*2))&
                  *localRHS(iNode*2-1)
             model%dof(iNodeID*4-1) = (model%dof(iNodeID*4-2)*localRHS(iNode*2-1)+model%dof(iNodeID*4-1)*localRHS(iNode*2))&
                  *localRHS(iNode*2)
          end if
       end do
       deallocate(localRHS)
    end do
    deallocate(position)
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
