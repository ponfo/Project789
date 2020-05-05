module IntegratorPtrM
  use UtilitiesM
  use DebuggerM
  use IntegratorM
  implicit none
  private
  public :: IntegratorPtrDT
  type :: IntegratorPtrDT
     class(IntegratorDT), pointer :: ptr
   contains
     generic           :: associate => associateWithIntegrator  &
                                     , associateWithIntegratorPtr
     procedure, public :: associateWithIntegrator
     procedure, public :: associateWithIntegratorPtr

     procedure, public :: getGaussOrder
     procedure, public :: getIntegTerms
     procedure, public :: getWeight
     procedure, public :: getGPoint
     procedure, public :: getShapeFunc
     procedure, public :: getDShapeFunc
     procedure, public :: getDDShapeFunc
     procedure, public :: getWeightFull
     procedure, public :: getGPointFull
     procedure, public :: getShapeFuncFull
     procedure, public :: getDShapeFuncFull
     procedure, public :: getDDShapeFuncFull
  end type IntegratorPtrDT

contains

  subroutine associateWithIntegrator(this, integrator)
    implicit none
    class(IntegratorPtrDT), intent(inout) :: this
    type(IntegratorDT), target, intent(in) :: integrator
    this%ptr => integrator
  end subroutine associateWithIntegrator

  subroutine associateWithIntegratorPtr(this, integratorPtr)
    implicit none
    class(IntegratorPtrDT), intent(inout) :: this
    type(IntegratorPtrDT), target, intent(in) :: integratorPtr
    this%ptr => integratorPtr%ptr
  end subroutine associateWithIntegratorPtr

  integer(ikind) pure function getGaussOrder(this)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    getGaussOrder = this%ptr%getGaussOrder()
  end function getGaussOrder

  integer(ikind) pure function getIntegTerms(this)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    getIntegTerms = this%ptr%getIntegTerms()
  end function getIntegTerms

  real(rkind) pure function getWeight(this, i)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    integer(ikind)        , intent(in) :: i
    getWeight = this%ptr%getWeight(i)
  end function getWeight

  real(rkind) pure function getGPoint(this, i, j)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    integer(ikind)        , intent(in) :: i
    integer(ikind)        , intent(in) :: j
    getGPoint = this%ptr%getGPoint(i,j)
  end function getGPoint

  real(rkind) pure function getShapeFunc(this, i, j)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    integer(ikind)        , intent(in) :: i
    integer(ikind)        , intent(in) :: j
    getShapeFunc = this%ptr%getShapeFunc(i,j)
  end function getShapeFunc

  real(rkind) pure function getDShapeFunc(this, i, j, k)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    integer(ikind)        , intent(in) :: i
    integer(ikind)        , intent(in) :: j
    integer(ikind)        , intent(in) :: k
    getDShapeFunc = this%ptr%getDShapeFunc(i,j,k)
  end function getDShapeFunc

  real(rkind) pure function getDDShapeFunc(this, i, j, k, l)
    implicit none
    class(IntegratorPtrDT), intent(in) :: this
    integer(ikind)        , intent(in) :: i
    integer(ikind)        , intent(in) :: j
    integer(ikind)        , intent(in) :: k
    integer(ikind)        , intent(in) :: l
    getDDShapeFunc = this%ptr%getDDShapeFunc(i,j,k,l)
  end function getDDShapeFunc

  pure function getWeightFull(this)
    implicit none
    class(IntegratorPtrDT)              , intent(in)  :: this
    real(rkind)           , dimension(:), allocatable :: getWeightFull
    getWeightFull = this%ptr%getWeightFull()
  end function getWeightFull

  pure function getGPointFull(this)
    implicit none
    class(IntegratorPtrDT)                , intent(in)  :: this
    real(rkind)           , dimension(:,:), allocatable :: getGPointFull
    getGPointFull = this%ptr%getGPointFull()
  end function getGPointFull

  pure function getShapeFuncFull(this)
    implicit none
    class(IntegratorPtrDT)                , intent(in)  :: this
    real(rkind)           , dimension(:,:), allocatable :: getShapeFuncFull
    getShapeFuncFull = this%ptr%getShapeFuncFull()
  end function getShapeFuncFull

  pure function getDShapeFuncFull(this)
    implicit none
    class(IntegratorPtrDT)                  , intent(in)  :: this
    real(rkind)           , dimension(:,:,:), allocatable :: getDShapeFuncFull
    getDShapeFuncFull = this%ptr%getDShapeFuncFull()
  end function getDShapeFuncFull

  pure function getDDShapeFuncFull(this)
    implicit none
    class(IntegratorPtrDT)                    , intent(in)  :: this
    real(rkind)           , dimension(:,:,:,:), allocatable :: getDDShapeFuncFull
    getDDShapeFuncFull = this%ptr%getDDShapeFuncFull()
  end function getDDShapeFuncFull
  
end module IntegratorPtrM
