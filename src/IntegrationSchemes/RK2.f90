module RK2M

  use UtilitiesM

  use ProcessM
  
  use SchemeM
  
  use IntegrandM
  
  implicit none
  
  private
  public :: RK2DT
  
  type, extends(NewSchemeDT) :: RK2DT
   contains
     procedure, nopass :: integrate
  end type RK2DT
  
contains
  
  subroutine integrate(this,dt)
    class(NewProcessDT), intent(inout) :: this
    class(IntegrandDT), allocatable   :: this_half
    real(rkind)       , intent(in)    :: dt
    select type (this)
    class is (IntegrandDT)
       allocate(this_half,source=this)
       this_half = this + this%t(this%state)*(0.5*dt)
       this      = this + this_half*dt
       class default
       stop 'integrate: unsupported class'
    end select
  end subroutine integrate
  
end module RK2M
