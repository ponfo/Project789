module CFDBuilderAndSolverM

  use UtilitiesM
  use SparseKit
  use DebuggerM

  use SparseKit

  use NodePtrM
  use ElementPtrM
  use ConditionPtrM

  use CFDElementM 

  use LeftHandSideM

  use CFDModelM
  use CFDApplicationM

  use BuilderAndSolverM

  implicit none

  private
  public :: CFDBuilderAndSolverDT

  type, extends(NewBuilderAndSolverDT) :: CFDBuilderAndSolverDT
   contains
     procedure :: buildAndSolve
     procedure :: update
  end type CFDBuilderAndSolverDT
  
contains

  subroutine buildAndSolve(this, app)
    implicit none
    class(CFDBuilderAndSolverDT), intent(inout) :: this
    class(CFDApplicationDT)     , intent(inout) :: app
    call calculateStab(app)
    call assembleSystem(app)
  end subroutine buildAndSolve

  subroutine update(this, app)
    implicit none
    class(CFDBuilderAndSolverDT), intent(inout) :: this
    class(CFDApplicationDT)     , intent(inout) :: app
    call applyDirichlet(app)
  end subroutine update

  subroutine assembleSystem(app)
    implicit none
    class(CFDApplicationDT), intent(inout)     :: app
    integer(ikind)                             :: iNode, iNodeID
    integer(ikind)                             :: iElem, nElem, nNode
    type(LeftHandSideDT)                       :: localLHS
    real(rkind), dimension(:)    , allocatable :: localRHS
    type(ElementPtrDT)                         :: element
    app%model%rhs = 0._rkind
    nElem = app%model%getnElement()
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(element, nNode, iNode, iNodeID, localRHS, localLHS) &
    !$OMP SHARED(app, nElem)
    do iElem = 1, nElem
       element = app%model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(app%model%processInfo, localLHS, localRHS)
       do iNode = 1, nNode
          iNodeID = element%getNodeID(iNode)
          app%model%rhs(iNodeID*4-3) = app%model%rhs(iNodeID*4-3) &
               + localRHS(iNode*4-3)
          app%model%rhs(iNodeID*4-2) = app%model%rhs(iNodeID*4-2) &
               + localRHS(iNode*4-2)
          app%model%rhs(iNodeID*4-1) = app%model%rhs(iNodeID*4-1) &
               + localRHS(iNode*4-1)
          app%model%rhs(iNodeID*4) = app%model%rhs(iNodeID*4) &
               + localRHS(iNode*4)
       end do
       deallocate(localRHS)
    end do
    !$OMP END PARALLEL DO
  end subroutine assembleSystem
  
  subroutine calculateStab(app)
    implicit none
    class(CFDApplicationDT), intent(inout) :: app
    integer(ikind)                         :: iElem, nElem
    call zeroStabMat()
    nElem = app%model%getnElement()
    do iElem = 1, nElem
       call app%element(iElem)%calculateStabMat(app%model%processInfo)
    end do
  end subroutine calculateStab

  subroutine applyDirichlet(app)
    implicit none
    class(CFDApplicationDT)             , intent(inout) :: app
    integer(ikind)                                      :: iCond, nCond, iNode, iNodeID
    real(rkind)             , dimension(:), allocatable :: localRHS
    real(rkind)                                         :: Vx, Vy, filter, nx, ny
    integer(ikind)                                      :: i, nNode, nodeID, nDof
    integer(ikind)          , dimension(:), allocatable :: position
    type(NodePtrDT)                                     :: node
    type(ConditionPtrDT)                                :: condition
    nCond = app%model%getnCondition()
    nNode = app%model%getnNode()
    nDof = 4
    do i = 1, nNode
       node = app%model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          app%model%dof(nodeID*4-3) = node%ptr%dof(1)%fixedVal
       end if
       if(node%ptr%dof(2)%isFixed) then
          nodeID = node%ptr%getID()
          app%model%dof(nodeID*4-2) = node%ptr%dof(2)%fixedVal
       end if
       if(node%ptr%dof(3)%isFixed) then
          nodeID = node%ptr%getID()
          app%model%dof(nodeID*4-1) = node%ptr%dof(3)%fixedVal
       end if
       if(node%ptr%dof(4)%isFixed) then
          nodeID = node%ptr%getID()
          app%model%dof(nodeID*4) = node%ptr%dof(4)%fixedVal
       end if
    end do
    allocate(position(nNode))
    position = 0
    do iCond = 1, nCond
       condition = app%model%getCondition(iCond)
       nNode = condition%getnNode()
       call condition%calculateRHS(app%model%processInfo, localRHS)
       do iNode = 1, nNode
          iNodeID = condition%getNodeID(iNode)
          position(iNodeID) = position(iNodeID) + 1
          if (position(iNodeID) .le. 1) then
             i = i + 1
             Vx  = app%model%dof(iNodeID*4-2)
             Vy  = app%model%dof(iNodeID*4-1)
             nx  = localRHS(iNode*2-1)
             ny  = localRHS(iNode*2  )
             app%model%dof(iNodeID*4-2) = (Vx*nx + Vy*ny)*nx
             app%model%dof(iNodeID*4-1) = (Vx*nx + Vy*ny)*ny
          end if
       end do
       deallocate(localRHS)
    end do
    deallocate(position)
  end subroutine applyDirichlet

end module CFDBuilderAndSolverM
