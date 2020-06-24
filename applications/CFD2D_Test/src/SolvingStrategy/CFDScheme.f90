module CFDSchemeM

  use UtilitiesM
  
  use CFDModelM
  use CFDApplicationM
  
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

  subroutine calculateOutputs(this, app)
    implicit none
    class(CFDSchemeDT)     , intent(inout) :: this
    class(CFDApplicationDT), intent(inout) :: app
    integer(ikind)                         :: dim, iNode
    real(rkind)                            :: rho, Vx, Vy, T, P, E, M
    real(rkind)                            :: R, Cv, Vc, gamma
    dim = app%model%getnNode()
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
    app%model%results%density(:)        = 0._rkind
    app%model%results%velocity(:,1)     = 0._rkind
    app%model%results%velocity(:,2)     = 0._rkind
    app%model%results%internalEnergy(:) = 0._rkind
    app%model%results%temperature(:)    = 0._rkind
    app%model%results%mach(:)           = 0._rkind
    app%model%results%pressure(:)       = 0._rkind
    R     = app%model%processInfo%getConstants(3)
    Cv    = app%model%processInfo%getConstants(4)
    Vc    = app%model%processInfo%getConstants(5)
    gamma = app%model%processInfo%getConstants(6)
    !$OMP PARALLEL DO PRIVATE(iNode,rho,Vx,Vy,E,T,M,P)
    do iNode = 1, app%model%getnNode()
       rho = app%model%dof(iNode*4-3)
       Vx  = app%model%dof(iNode*4-2)/rho
       Vy  = app%model%dof(iNode*4-1)/rho
       E   = app%model%dof(iNode*4)/rho
       T   = (E-0.5d0*(Vx**2+Vy**2))/Cv
       M   = sqrt(Vx**2+Vy**2)/Vc
       P   = rho*R*T
       app%model%results%density(iNode)        = rho
       app%model%results%velocity(iNode,1)     = Vx
       app%model%results%velocity(iNode,2)     = Vy
       app%model%results%internalEnergy(iNode) = E
       app%model%results%temperature(iNode)    = T
       app%model%results%mach(iNode)           = M
       app%model%results%pressure(iNode)       = P
    end do
    !$OMP END PARALLEL DO
  end subroutine calculateOutputs
  
  subroutine integrator(this, dt, multi_step)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
    logical            , intent(in) :: multi_step
  end subroutine integrator

end module CFDSchemeM
