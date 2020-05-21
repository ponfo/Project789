module CFDSchemeM

  use UtilitiesM
  
  use CFDModelM

  use ElementPtrM
  
  use SchemeM

  use ProcessM

  implicit none

  private
  public :: CFDSchemeDT

  type, extends(NewSchemeDT) :: CFDSchemeDT
   contains
     procedure, public :: calculateOutputs
     procedure, nopass :: integrate => integrator
  end type CFDSchemeDT
  
contains

  subroutine calculateOutputs(this, model)
    implicit none
    class(CFDSchemeDT), intent(inout)          :: this
    class(CFDModelDT) , intent(inout)          :: model
    integer(ikind)                             :: iElem, iGauss, nNode
    integer(ikind)                             :: nElem, nGauss, dim
    real(rkind), dimension(:,:,:), allocatable :: localResultMat
    type(ElementPtrDT)                         :: element
    dim = model%getnNode()
    allocate(model%results%velocity(dim,dim))
    allocate(model%results%density(dim))
    allocate(model%results%mach(dim))
    allocate(model%results%pressure(dim))
    allocate(model%results%temperature(dim))
    allocate(model%results%internalEnergy(dim))
    nElem = model%getnElement()
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(model%processInfo, localResultMat)
       model%results%velocity(dim,dim)   = localResultMat()
       model%results%density(dim)        = localResultMat()
       model%results%mach(dim)           = localResultMat()
       model%results%pressure(dim)       = localResultMat()
       model%results%temperature(dim)    = localResultMat()
       model%results%internalEnergy(dim) = localResultMat()
    end do
  end subroutine calculateOutputs
  
  subroutine integrator(this, dt)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
  end subroutine integrator

end module CFDSchemeM



!!$    integer(ikind)                    :: i
!!$    integer(ikind)                    :: n
!!$    n = model%getnNode()
!!$    allocate(model%results%density(n)       , model%results%velocity(n,n) )
!!$    allocate(model%results%internalEnergy(n), model%results%temperature(n))
!!$    allocate(model%results%pressure(n)      , model%results%mach(n)       )
!!$    do i = 1, size(model%dof),4
!!$       model%results%density(i)        = model%dof(i)
!!$       model%results%velocity(i,i)     = (/model%dof(i+1),model%dof(i+2)/)/model%dof(i)
!!$       model%results%internalEnergy(i) = model%dof(i+3)/model%dof(i)
!!$       model%results%temperature(i)    = &
!!$            (model%gamma-1)/R*(model%dof(i+3)-0.5*(model%dof(i+1)**2+model%dof(i+2)**2))
!!$       model%results%pressure(i)       = &
!!$            (model%gamma-1)*model%dof(i)*(model%dof(i+3)-0.5*(model%dof(i+1)**2+model%dof(i+2)**2))
!!$       model%results%mach(i)           = sqrt(model%dof(i+1)**2+model%dof(i+2)**2)         &
!!$            /sqrt(model%gamma*model%R                                                   &
!!$            *(model%gamma-1)/R*(model%dof(i+3)-0.5*(model%dof(i+1)**2+model%dof(i+2)**2)))
!!$    end do
