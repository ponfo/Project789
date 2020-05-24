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
    integer(ikind), dimension(:), allocatable  :: counter
    integer(ikind)                             :: iElem, nNode
    integer(ikind)                             :: nElem, dim, iNode
    integer(ikind)                             :: position
    real(rkind), dimension(:,:,:), allocatable :: localResultMat
    type(ElementPtrDT)                         :: element
    dim = model%getnNode()
    if (allocated(model%results%velocity)) then
       deallocate(model%results%velocity)
       deallocate(model%results%density)
       deallocate(model%results%mach)
       deallocate(model%results%pressure)
       deallocate(model%results%temperature)
       deallocate(model%results%internalEnergy)
       allocate(model%results%velocity(dim,2))
       allocate(model%results%density(dim))
       allocate(model%results%mach(dim))
       allocate(model%results%pressure(dim))
       allocate(model%results%temperature(dim))
       allocate(model%results%internalEnergy(dim))
       allocate(counter(dim))
    else
       allocate(model%results%velocity(dim,2))
       allocate(model%results%density(dim))
       allocate(model%results%mach(dim))
       allocate(model%results%pressure(dim))
       allocate(model%results%temperature(dim))
       allocate(model%results%internalEnergy(dim))
       allocate(counter(dim))
    end if
    counter = 0
    nElem = model%getnElement()
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(model%processInfo, localResultMat)
       do iNode = 1, nNode
          position                              = localResultMat(iNode,1,1)
          model%results%velocity(position,1)     = localResultMat(iNode,2,1)
          model%results%velocity(position,2)     = localResultMat(iNode,3,1)
          model%results%density(position)        = localResultMat(iNode,4,1)
          model%results%mach(position)           = localResultMat(iNode,5,1)
          model%results%pressure(position)       = localResultMat(iNode,6,1)
          model%results%temperature(position)    = localResultMat(iNode,7,1)
          model%results%internalEnergy(position) = localResultMat(iNode,8,1)
          counter(position) = counter(position) + 1
       end do
       deallocate(localResultMat)
    end do
    do iNode = 1, dim
       model%results%velocity(iNode,1)     = model%results%velocity(iNode,1)/counter(iNode)
       model%results%velocity(iNode,2)     = model%results%velocity(iNode,2)/counter(iNode)
       model%results%density(iNode)        = model%results%density(iNode)/counter(iNode) 
       model%results%mach(iNode)           = model%results%mach(iNode)/counter(iNode)
       model%results%pressure(iNode)       = model%results%pressure(iNode)/counter(iNode) 
       model%results%temperature(iNode)    = model%results%temperature(iNode)/counter(iNode)
       model%results%internalEnergy(iNode) = model%results%internalEnergy(iNode)/counter(iNode)
    end do
    deallocate(counter)
  end subroutine calculateOutputs
  
  subroutine integrator(this, dt)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
  end subroutine integrator

end module CFDSchemeM



