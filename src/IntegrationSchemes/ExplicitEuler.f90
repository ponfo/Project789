module ExplicitEulerM

  use UtilitiesM
  
  use ProcessM
 
  use SchemeM
  
  use IntegrandM
  
  implicit none
  
  private
  public :: ExplicitEulerDT
  
  type, extends(NewSchemeDT) :: ExplicitEulerDT
   contains
     procedure, nopass :: integrate
  end type ExplicitEulerDT
  
contains
  
  subroutine integrate(this,dt) 
    class(NewProcessDT), intent(inout) :: this 
    real(rkind), intent(in) :: dt
    select type (this)
    class is (IntegrandDT)
       this = this + this%t(this%state)*dt
       class default
       stop 'integrate: unsupported class.'
    end select
  end subroutine integrate
  
end module ExplicitEulerM
