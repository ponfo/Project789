module CFDSchemeM

  use UtilitiesM
  
  use CFDModelM
  
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
    class(CFDSchemeDT), intent(inout) :: this
    class(CFDModelDT) , intent(inout) :: model
    integer(ikind)                    :: dim, iNode
    real(rkind)                       :: rho, Vx, Vy, T, P, E, M
    real(rkind)                       :: R, Cv, Vc, gamma
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
    else
       allocate(model%results%velocity(dim,2))
       allocate(model%results%density(dim))
       allocate(model%results%mach(dim))
       allocate(model%results%pressure(dim))
       allocate(model%results%temperature(dim))
       allocate(model%results%internalEnergy(dim))
    end if
    model%results%density(:)        = 0._rkind
    model%results%velocity(:,1)     = 0._rkind
    model%results%velocity(:,2)     = 0._rkind
    model%results%internalEnergy(:) = 0._rkind
    model%results%temperature(:)    = 0._rkind
    model%results%mach(:)           = 0._rkind
    model%results%pressure(:)       = 0._rkind
    R     = model%processInfo%getConstants(3)
    Cv    = model%processInfo%getConstants(4)
    Vc    = model%processInfo%getConstants(5)
    gamma = model%processInfo%getConstants(6)
    !$OMP PARALLEL DO PRIVATE(iNode,rho,Vx,Vy,E,T,M,P)
    do iNode = 1, model%getnNode()
       rho = model%dof(1,iNode)
       Vx  = model%dof(2,iNode)/model%dof(1,iNode)
       Vy  = model%dof(3,iNode)/model%dof(1,iNode)
       E   = model%dof(4,iNode)/model%dof(1,iNode)
       T   = (E-0.5d0*(Vx**2+Vy**2))/Cv
       M   = sqrt(Vx**2+Vy**2)/Vc
       P   = rho*R*T
       model%results%density(iNode)        = rho
       model%results%velocity(iNode,1)     = Vx
       model%results%velocity(iNode,2)     = Vy
       model%results%internalEnergy(iNode) = E
       model%results%temperature(iNode)    = T
       model%results%mach(iNode)           = M
       model%results%pressure(iNode)       = P
    end do
    !$OMP END PARALLEL DO
  end subroutine calculateOutputs
  
  subroutine integrator(this, dt)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
  end subroutine integrator

end module CFDSchemeM



