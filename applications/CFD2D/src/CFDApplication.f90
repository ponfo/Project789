module CFDApplicationM
  
  use UtilitiesM
  use DebuggerM

  use SourceM
  use NodeM
  
  use CFDElementM
  use NormalVelocityM
  use CFDMaterialM
  use CFDModelM

  implicit none

  private
  public :: CFDApplicationDT, cfdApplication

  type :: CFDApplicationDT
     type(NodeDT)          , dimension(:), allocatable :: node
     type(CFDElementDT)    , dimension(:), allocatable :: element
     type(NormalVelocityDT), dimension(:), allocatable :: normalVelocity
     type(SourceDT)        , dimension(:), allocatable :: source
     type(CFDMaterialDT)   , dimension(:), allocatable :: material
     type(CFDModelDT)                                  :: model
   contains
     procedure, public :: init
     procedure, public :: setTransientValues
  end type CFDApplicationDT

  interface cfdApplication
     procedure :: constructor
  end interface cfdApplication

contains

  type(CFDApplicationDT) function constructor&
       (nNode, nElement, nNormalVelocity, nSource, nMaterial, nGauss)
    implicit none
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nNormalVelocity
    integer(ikind), intent(in) :: nSource
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nGauss
    call constructor%init(nNode, nElement, nNormalVelocity, nSource, nMaterial, nGauss)
  end function constructor

  subroutine init(this, nNode, nElement, nNormalVelocity, nSource, nMaterial, nGauss)
    implicit none
    class(CFDApplicationDT), intent(inout) :: this
    integer(ikind)         , intent(in)    :: nNode
    integer(ikind)         , intent(in)    :: nElement
    integer(ikind)         , intent(in)    :: nNormalVelocity
    integer(ikind)         , intent(in)    :: nSource
    integer(ikind)         , intent(in)    :: nMaterial
    integer(ikind)         , intent(in)    :: nGauss
    allocate(this%node(nNode))
    allocate(this%element(nElement))
    allocate(this%normalVelocity(nNormalVelocity))
    allocate(this%source(nSource))
    allocate(this%material(nMaterial))
    call initGeometries(nGauss)
    this%model = cfdModel(                &
           nDof =   4                    &
         , nnz = nElement*256            &
         , id = 1                        &
         , nNode = nNode                 &
         , nElement = nElement           &
         , nCondition = nNormalVelocity  )
  end subroutine init

  subroutine setTransientValues(this, printStep, t0, errorTol, maxIter&
       , fSafe, constant, R, Cv, Vc, gamma, Vx, Vy, T, P, rho, mach)
    implicit none
    class(CFDApplicationDT), intent(inout)  :: this
    integer(ikind)          , intent(in)    :: printStep
    real(rkind)             , intent(in)    :: t0
    real(rkind)             , intent(in)    :: errorTol
    integer(ikind)          , intent(in)    :: maxIter
    real(rkind)             , intent(in)    :: fSafe
    real(rkind)             , intent(in)    :: constant
    real(rkind)             , intent(in)    :: R
    real(rkind)             , intent(in)    :: Cv
    real(rkind)             , intent(in)    :: Vc
    real(rkind)             , intent(in)    :: gamma
    real(rkind)             , intent(in)    :: Vx
    real(rkind)             , intent(in)    :: Vy
    real(rkind)             , intent(in)    :: T
    real(rkind)             , intent(in)    :: P
    real(rkind)             , intent(in)    :: rho
    real(rkind)             , intent(in)    :: mach
    real(rkind), dimension(:), allocatable  :: vector
    allocate(vector(7))
    vector = (/fSafe, constant, R, Cv, Vc, gamma&
         , Vx, Vy, T, P, rho, mach/)
    call this%model%processInfo%setPrintStep(printStep)
    call this%model%processInfo%setT0(t0)
    call this%model%processInfo%setErrorTol(errorTol)
    call this%model%processInfo%setMaxIter(maxIter)
    call this%model%processInfo%setConstants(12, vector)
  end subroutine setTransientValues

end module CFDApplicationM
