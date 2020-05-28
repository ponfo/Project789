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

  real(rkind), dimension(:,:), allocatable, save :: sumTau
  
contains

  subroutine buildAndSolve(this, model)
    implicit none
    class(CFDBuilderAndSolverDT), intent(inout) :: this
    class(CFDModelDT)           , intent(inout) :: model
    integer(ikind) :: i
    if (model%processInfo%getStep()==0) call initDof(model)
    call applyNewmannAndDirichletDOF(model)
    call assembleSystem(model)
    call applyDirichletAndNeumannLHS(model)
  end subroutine buildAndSolve

  subroutine assembleSystem(model)
    implicit none
    class(CFDModelDT)          , intent(inout) :: model
    real(rkind), dimension(:,:,:), allocatable :: localResultMat
    real(rkind), dimension(:,:)  , allocatable :: A1
    real(rkind), dimension(:,:)  , allocatable :: A2
    real(rkind), dimension(:,:)  , allocatable :: MTau
    real(rkind), dimension(:,:)  , allocatable :: sumNdNA
    real(rkind), dimension(:,:)  , allocatable :: sumdN
    real(rkind), dimension(:,:)  , allocatable :: inv
    real(rkind), dimension(:,:)  , allocatable :: lhs
    integer(ikind)                             :: iNode, jNode, iDof, jDof, iNodeID, jNodeID
    integer(ikind)                             :: iElem, nElem, nNode, nDof, ii, jj, i, j
    type(LeftHandSideDT)                       :: localLHS
    real(rkind), dimension(:)    , allocatable :: localRHS, constants
    real(rkind)                                :: adder, v(2), Vx, Vy, gamma, E, nu, Tau(4,4) 
    type(ElementPtrDT)                         :: element
    allocate(constants(12+size(model%dof)))
    do i = 1, 12
       constants(i) = model%processInfo%getConstants(i)
    end do
    do i = 1, size(model%dof)
       constants(12+i) = model%dof(i)
    end do
    call model%processInfo%setConstants(12+size(model%dof),constants)
    call model%lhs%free()
    model%lhs = Sparse(nnz = model%getnElement()*256, rows = model%getnNode()*4)
    nElem = model%getnElement()
    nDof = 4
    if (model%processInfo%getStep() .eq. 0) then
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
                      adder = adder                                                             &
                           + localLHS%mass(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof))  
                   end do
                   call model%mass%append(                                                       &
                        val = adder                                                             &
                        , row = iNodeID*nDof-(nDof-iDof)                                        &
                        , col = iNodeID*nDof-(nDof-iDof)                                        )
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
       call model%mass%makeCRS()
    end if
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateLocalSystem(model%processInfo,localLHS, localRHS)
       gamma = model%processInfo%getConstants(6)
       call element%calculateResults(model%processInfo, localResultMat)
       tau = localResultMat(1:4,1:4,6)
       nu  = localResultMat(5,1,6)
       do iNode = 1, nNode
          allocate(A1(nNode*nDof,nNode*nDof))
          allocate(A2(nNode*nDof,nNode*nDof))
          allocate(MTau(nNode*nDof,nNode*nDof))
          if (model%processInfo%getStep() == 0) allocate(inv(nNode*nDof,nNode*nDof))
          inv = 0.d0
          iNodeID = element%getNodeID(iNode)
          if (model%dof(iNode*4-3) == 0.d0) then
             Vx = 0.d0
             Vy = 0.d0
             E  = 0.d0
          else
             Vx = model%dof(iNode*4-2)/model%dof(iNode*4-3)
             Vy = model%dof(iNode*4-1)/model%dof(iNode*4-3)
             E  = model%dof(iNode*4)/model%dof(iNode*4-3)
          end if
          v  = (/Vx,Vy/)
          do jNode = 1, nNode
             if (model%processInfo%getStep() == 0) then
                do iDof = 1, nDof
                   adder = 0._rkind
                   do jDof = 1, nDof
                      adder = adder                                                             &
                           + localLHS%mass(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof))  
                   end do
                   inv(iDof,iDof) = 1.d0/adder
                end do
             end if
             ii = nDof*iNode-3
             jj = nDof*jNode-3
             MTau(ii:ii+3,jj:jj+3) = Tau                    
             A1(ii:ii+3,jj:jj+3)   = MA1(v,gamma,E)
             A2(ii:ii+3,jj:jj+3)   = MA2(v,gamma,E)
          end do
          sumNdNA = matmul(localResultMat(:,:,1),A1) + matmul(localResultMat(:,:,2),A2)
          if (model%processInfo%getStep() == 0) then
             sumTau = matmul(MTau,(matmul(localResultMat(:,:,3),matmul(A1,A1))&
                  +matmul(localResultMat(:,:,4),matmul(A2,A2))&
                  +matmul(localResultMat(:,:,5),matmul(A1,A2))&
                  +matmul(localResultMat(:,:,5),matmul(A2,A1))))&
                  -matmul(inv, matmul(MTau,matmul(A1+A2,sumNdNA)))
          end if
          sumdN   = localResultMat(:,:,3) + localResultMat(:,:,4) + 2.d0 * localResultMat(:,:,5)
          lhs     = -sumNdNA - sumTau - nu * sumdN 
          do jNode = 1, nNode
             jNodeID = element%getNodeID(jNode)
             do iDof = 1, nDof
                do jDof = 1, nDof
                   call model%LHS%append(                                      &
                        val = lhs(iNode*nDof-(nDof-iDof),jNode*nDof-(nDof-jDof)) &
                        , row = iNodeID*nDof-(nDof-iDof)                      &
                        , col = jNodeID*nDof-(nDof-jDof)                      ) 
                end do
             end do
          end do
          deallocate(MTau, A1, A2)
          if (model%processInfo%getStep() == 0) deallocate(inv)
       end do
       call localLHS%free()
       deallocate(localRHS)
       deallocate(localResultMat)
    end do
    call model%lhs%makeCRS()
  end subroutine assembleSystem

  subroutine applyNewmannAndDirichletDOF(model)
    implicit none
    class(CFDModelDT)                   , intent(inout) :: model
    integer(ikind)                                      :: iCond, nCond, iNode, iNodeID
    integer(ikind)          , dimension(:), allocatable :: position 
    real(rkind)             , dimension(:), allocatable :: localRHS
    integer(ikind)                                      :: i, nNode, nodeID, nDof
    type(NodePtrDT)                                     :: node
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
    nNode = model%getnNode()
    nDof = 4
    do i = 1, nNode
       node = model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(nodeID*nDof-3) = node%ptr%dof(1)%fixedVal
       end if
       if(node%ptr%dof(2)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(nodeID*nDof-2) = node%ptr%dof(2)%fixedVal
       end if
       if(node%ptr%dof(3)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(nodeID*nDof-1) = node%ptr%dof(3)%fixedVal
       end if
       if(node%ptr%dof(4)%isFixed) then
          nodeID = node%ptr%getID()
          model%dof(nodeID*nDof  ) = node%ptr%dof(4)%fixedVal
       end if
    end do
  end subroutine applyNewmannAndDirichletDOF

  subroutine applyDirichletAndNeumannLHS(model)
    implicit none
    class(CFDModelDT), intent(inout)           :: model
    integer(ikind)                             :: i, j, nNode, nodeID, nDof, nCond, iCond
    integer(ikind)                             :: iNode, iNodeID , nElem, iElem
    type(ElementPtrDT)                         :: element
    type(NodePtrDT)                            :: node
    type(ConditionPtrDT)                       :: condition
    real(rkind)                                :: dt_min, alpha, deltaTU, dt, area   
    real(rkind)                                :: val1, val2, V, dt_elem, deltaTC    
    real(rkind)                                :: Vx, Vy, Vxmax, Vymax, fSafe, cota 
    type(NodePtrDT), dimension(:), allocatable :: nodalPoints
    nElem = model%getnElement()
    nNode = model%getnNode()
    nDof = 4
    do i = 1, nNode
       node = model%getNode(i)
       if(node%ptr%dof(1)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-3)
       end if
       if(node%ptr%dof(2)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-2)
       end if
       if(node%ptr%dof(3)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof-1)
       end if
       if(node%ptr%dof(4)%isFixed) then
          nodeID = node%ptr%getID()
          call model%lhs%setDirichlet(nodeID*nDof)
       end if
    end do
    nCond = model%getnCondition()
    do iCond = 1, nCond
       condition = model%getCondition(iCond)
       nNode = condition%getnNode()
       do iNode = 1, nNode
          iNodeID = condition%getNodeID(iNode)
          call model%lhs%setDirichlet(iNodeID*nDof-2)
          call model%lhs%setDirichlet(iNodeID*nDof-1)
       end do
    end do
  end subroutine applyDirichletAndNeumannLHS

  function MA1(v,gamma,E)
    implicit none
    real(rkind), dimension(4,4) :: MA1
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: gamma
    real(rkind)                 :: E
    MA1(1,1) = 0._rkind
    MA1(1,2) = 1._rkind
    MA1(1,3) = 0._rkind
    MA1(1,4) = 0._rkind

    MA1(2,1) = (((gamma-1._rkind)/2._rkind)*(v(1)**2+v(2)**2))-v(1)**2
    MA1(2,2) = (3._rkind-gamma)*v(1)
    MA1(2,3) = -(gamma-1._rkind)*v(2)
    MA1(2,4) = (gamma-1._rkind)

    MA1(3,1) = -v(1)*v(2)
    MA1(3,2) = v(2)
    MA1(3,3) = v(1)
    MA1(3,4) = 0._rkind

    MA1(4,1) = (((gamma-1._rkind)*(v(1)**2+v(2)**2))-(gamma*E))*v(1)
    MA1(4,2) = (gamma*E)-((gamma-1._rkind)*(v(1)**2+v(2)**2)/2._rkind)&
         -((gamma-1._rkind)*v(1)*v(2))
    MA1(4,3) = -(gamma-1._rkind)*v(1)*v(2)
    MA1(4,4) = gamma*v(1)
  end function MA1

  function MA2(v,gamma,E)
    implicit none
    real(rkind), dimension(4,4) :: MA2
    real(rkind), dimension(2)   :: v
    real(rkind)                 :: gamma
    real(rkind)                 :: E
    MA2(1,1) = 0._rkind
    MA2(1,2) = 0._rkind
    MA2(1,3) = 1._rkind
    MA2(1,4) = 0._rkind

    MA2(2,1) = -v(1)*v(2)
    MA2(2,2) = v(2)
    MA2(2,3) = v(1)
    MA2(2,4) = 0._rkind

    MA2(3,1) = (((gamma-1._rkind)/2._rkind)*(v(1)**2+v(2)**2))-v(2)**2
    MA2(3,2) = -(gamma-1._rkind)*v(1)
    MA2(3,3) = (3._rkind-gamma)*v(2)
    MA2(3,4) = (gamma-1._rkind)

    MA2(4,1) = (((gamma-1._rkind)*(v(1)**2+v(2)**2))-(gamma*E))*v(2)
    MA2(4,2) = -(gamma-1._rkind)*v(1)*v(2)
    MA2(4,3) = (gamma*E)-((gamma-1._rkind)*(v(1)**2+v(2)**2)/2._rkind)&
         -((gamma-1._rkind)*v(2)**2)
    MA2(4,4) = gamma*v(2)
  end function MA2

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
    do i = 1, model%getnNode()
       model%dof(i*4-3) = rho
       model%dof(i*4-2) = rho*Vx
       model%dof(i*4-1) = rho*Vy
       model%dof(i*4  ) = rho*(Cv*T+0.5d0*(Vx**2+Vy**2)) 
    end do
  end subroutine initDof
  
end module CFDBuilderAndSolverM
