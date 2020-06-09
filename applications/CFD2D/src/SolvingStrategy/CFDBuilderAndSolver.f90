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
  public :: CFDBuilderAndSolverDT, applyDirichlet, calculateDt, calculateMass

  type, extends(NewBuilderAndSolverDT) :: CFDBuilderAndSolverDT
   contains
     procedure :: buildAndSolve
  end type CFDBuilderAndSolverDT
  
contains

  subroutine buildAndSolve(this, model)
    implicit none
    class(CFDBuilderAndSolverDT), intent(inout) :: this
    class(CFDModelDT)           , intent(inout) :: model
    call calculateStab(model)
    call assembleSystem(model)
  end subroutine buildAndSolve

  subroutine assembleSystem(model)
    implicit none
    class(CFDModelDT)      , intent(inout)     :: model
    integer(ikind)                             :: iNode, i, iNodeID, j
    integer(ikind)                             :: iElem, nElem, nNode
    type(LeftHandSideDT)                       :: localLHS
    real(rkind), dimension(:)    , allocatable :: localRHS
    real(rkind), dimension(:,:,:), allocatable :: resultMat 
    type(ElementPtrDT)                         :: element
    nElem = model%getnElement()
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(model%processInfo,localLHS, localRHS)
       do iNode = 1, nNode
          iNodeID = element%getNodeID(iNode)
          model%rhs(1:4,iNodeID) = model%rhs(1:4,iNodeID) + localLHS%stiffness(1:4,iNode)
       end do
       call localLHS%free()
    end do
  end subroutine assembleSystem

  subroutine calculateStab(model)
    implicit none
    class(CFDModelDT)          , intent(inout) :: model
    integer(ikind)                             :: iNode, i, iNodeID, j
    integer(ikind)                             :: iElem, nElem
    real(rkind), dimension(:,:,:), allocatable :: resultMat 
    type(ElementPtrDT)                         :: element
    nElem = model%getnElement()
    if (model%processInfo%getProcess(3) == 1) then
       model%processInfo%mat = 0._rkind
       do iElem = 1, nElem
          element = model%getElement(iElem)
          call element%calculateResults(model%processInfo,resultMat)
       end do
    end if
  end subroutine calculateStab
  
  subroutine calculateMass(model)
    implicit none
    class(CFDModelDT)          , intent(inout) :: model
    integer(ikind)                             :: iElem, nElem
    real(rkind), dimension(:,:,:), allocatable :: resultMat 
    type(ElementPtrDT)                         :: element
    nElem = model%getnElement()
    model%processInfo%vect = 0._rkind
    do iElem = 1, nElem
       element = model%getElement(iElem)
       call element%calculateResults(model%processInfo,resultMat)
    end do
  end subroutine calculateMass
  
  subroutine calculateDt(model)
    implicit none
    class(CFDModelDT)          , intent(inout) :: model
    integer(ikind)                             :: iElem, nElem
    real(rkind), dimension(:,:,:), allocatable :: resultMat 
    type(ElementPtrDT)                         :: element
    nElem = model%getnElement()
    do iElem = 1, nElem
       element = model%getElement(iElem)
       call element%calculateResults(model%processInfo,resultMat)
    end do
  end subroutine calculateDt
  
  subroutine applyDirichlet(model)
    implicit none
    class(CFDModelDT)                   , intent(inout) :: model
    integer(ikind)                                      :: iCond, nCond, iNode, iNodeID
    real(rkind)             , dimension(:), allocatable :: localRHS
    real(rkind)                                         :: Vx, Vy, filter, nx, ny
    integer(ikind)                                      :: i, nNode, nodeID, nDof
    integer(ikind)          , dimension(:), allocatable :: position
    type(NodePtrDT)                                     :: node
    type(ConditionPtrDT)                                :: condition
    nCond = model%getnCondition()
    nNode = model%getnNode()
    nDof = 4
    do i = 1, nNode
       node = model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(1,nodeID) = node%ptr%dof(1)%fixedVal
       end if
       if(node%ptr%dof(2)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(2,nodeID) = node%ptr%dof(2)%fixedVal
       end if
       if(node%ptr%dof(3)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(3,nodeID) = node%ptr%dof(3)%fixedVal
       end if
       if(node%ptr%dof(4)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(4,nodeID) = node%ptr%dof(4)%fixedVal
       end if
    end do
    allocate(position(nNode))
    position = 0
    do iCond = 1, nCond
       condition = model%getCondition(iCond)
       nNode = condition%getnNode()
       call condition%calculateRHS(model%processInfo, localRHS)
       do iNode = 1, nNode
          iNodeID = condition%getNodeID(iNode)
          position(iNodeID) = position(iNodeID) + 1
          if (position(iNodeID) .le. 1) then
             i = i + 1
             Vx  = model%dof(2,iNodeID)
             Vy  = model%dof(3,iNodeID)
             nx  = localRHS(iNode*2-1)
             ny  = localRHS(iNode*2  )
             model%dof(2,iNodeID) = (Vx*nx + Vy*ny)*nx
             model%dof(3,iNodeID) = (Vx*nx + Vy*ny)*ny
          end if
       end do
       deallocate(localRHS)
    end do
    deallocate(position)
  end subroutine applyDirichlet
  
end module CFDBuilderAndSolverM
