module CFDSchemeM

  use UtilitiesM
  
  use ElementPtrM
  use CFDApplicationM
  
  use SchemeM

  use ProcessM

  implicit none

  private
  public :: CFDSchemeDT

  type, extends(NewSchemeDT) :: CFDSchemeDT
   contains
     procedure :: calculateOutputs
  end type CFDSchemeDT
  
contains

  subroutine calculateOutputs(this, app)
    implicit none
    class(CFDSchemeDT)     , intent(inout)     :: this
    class(CFDApplicationDT), intent(inout)     :: app
    integer(ikind)                             :: iNode, iNodeID, dim
    integer(ikind)                             :: iElem, nElem, nNode
    real(rkind), dimension(:,:,:), allocatable :: resultMat
    type(ElementPtrDT)                         :: element
    dim   = app%model%getnNode()
    nElem = app%model%getnElement()
    if (allocated(app%model%results%velocity)) then
       deallocate(app%model%results%velocity)
       deallocate(app%model%results%density)
       deallocate(app%model%results%mach)
       deallocate(app%model%results%pressure)
       deallocate(app%model%results%temperature)
       deallocate(app%model%results%internalEnergy)
       allocate(app%model%results%velocity(dim,2))
       allocate(app%model%results%density(dim))
       allocate(app%model%results%mach(dim))
       allocate(app%model%results%pressure(dim))
       allocate(app%model%results%temperature(dim))
       allocate(app%model%results%internalEnergy(dim))
    else
       allocate(app%model%results%velocity(dim,2))
       allocate(app%model%results%density(dim))
       allocate(app%model%results%mach(dim))
       allocate(app%model%results%pressure(dim))
       allocate(app%model%results%temperature(dim))
       allocate(app%model%results%internalEnergy(dim))
    end if
    app%model%results%density(:)        = 0.d0
    app%model%results%velocity(:,1)     = 0.d0
    app%model%results%velocity(:,2)     = 0.d0
    app%model%results%internalEnergy(:) = 0.d0
    app%model%results%temperature(:)    = 0.d0
    app%model%results%mach(:)           = 0.d0
    app%model%results%pressure(:)       = 0.d0
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(element, nNode, iNode, iNodeID, resultMat) &
    !$OMP SHARED(app, nElem)
    do iElem = 1, nElem
       element = app%model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(app%model%processInfo, resultMat)
       do iNode = 1, nNode
          iNodeID = element%getNodeID(iNode)
          app%model%results%density(iNodeID)        = resultMat(1,iNode,1)
          app%model%results%velocity(iNodeID,1)     = resultMat(2,iNode,1)
          app%model%results%velocity(iNodeID,2)     = resultMat(3,iNode,1)
          app%model%results%internalEnergy(iNodeID) = resultMat(4,iNode,1)
          app%model%results%temperature(iNodeID)    = resultMat(5,iNode,1)
          app%model%results%mach(iNodeID)           = resultMat(6,iNode,1)
          app%model%results%pressure(iNodeID)       = resultMat(7,iNode,1)
       end do
       deallocate(resultMat)
    end do
    !$OMP END PARALLEL DO
  end subroutine calculateOutputs

end module CFDSchemeM



