module DirectBuilderAndSolverM

  use UtilitiesM
  use DebuggerM

  use Element1DPtrM
  use Element2DPtrM

  use SparseKit

  use ThermalModelM

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
    class(ThermalModelDT)          , intent(inout) :: model
    write(*,*) '*** Direct Builder And Solver ***'
    call debugLog('  Assembling stiffness matrix and right hand side vector')
    print'(A)', 'Assembling stiffness matrix and right hand side vector'
    call assembleStiffness(model)
    call model%applySource(model%rhs)
    call model%applyBC2D(model%stiffness, model%rhs)
    call model%applyBC1D(model%stiffness, model%rhs)
    write(*,*) '*** Init Linear Solver ***'
    call solve(model)
  end subroutine buildAndSolve

  subroutine assembleStiffness(model)
    implicit none
    class(ThermalModelDT), intent(inout) :: model
    type(Element1DPtrTYPE)               :: element1D
    type(Element2DPtrTYPE)               :: element2D
    integer(ikind) :: i, j, iElem, nElem1D, nElem2D, nPoint
    real(rkind), dimension(:,:), allocatable :: localStiffness
    nElem1D = model%nLine
    do iElem = 1, nElem1D
       element1D = model%elementList1D%getElement(iElem)
       nPoint = element1D%getnPoint()
       allocate(localStiffness(nPoint,nPoint))
       localStiffness = element1D%getStiffness()
       do i = 1, nPoint
          do j = 1, nPoint
             call model%stiffness%append(val = localStiffness(i,j) &
                  , row = element1D%getPointID(i)                 &
                  , col = element1D%getPointID(j)                 )
          end do
       end do
       deallocate(localStiffness)
    end do
    nElem2D = model%nTriang + model%nQuad
    do iElem = 1, nElem2D
       element2D = model%elementList2D%getElement(iElem)
       nPoint = element2D%getnPoint()
       allocate(localStiffness(nPoint,nPoint))
       localStiffness = element2D%getStiffness()
       do i = 1, nPoint
          do j = 1, nPoint
             call model%stiffness%append(val = localStiffness(i,j) &
                  , row = element2D%getPointID(i)                 &
                  , col = element2D%getPointID(j)                 )
          end do
       end do
       deallocate(localStiffness)
    end do
    call model%stiffness%makeCRS()
  end subroutine assembleStiffness

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

    call UserSolver%solve(model%rhs, model%stiffness &
         , model%dof, data)

    call cpu_time(finish)
    print'(a,e14.7)', 'solver time => ', (finish-start)
    call debuglog('done solving')
  end subroutine solve

end module DirectBuilderAndSolverM
