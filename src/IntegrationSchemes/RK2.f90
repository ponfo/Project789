module RK2M

  use UtilitiesM

  use ProcessM

  use BaseIntegrandM

  use IntegrandM

  implicit none

  private
  public :: RK2DT

  type, extends(BaseIntegrandDT) :: RK2DT
   contains
     procedure, nopass :: integrate
  end type RK2DT

contains

  subroutine integrate(this,dt,multi_step)
    class(NewProcessDT), intent(inout) :: this
    class(IntegrandDT), allocatable    :: this_half
    real(rkind)       , intent(in)     :: dt
    logical           , intent(in) :: multi_step
    select type (this)
    class is (IntegrandDT)
       allocate(this_half,source=this)
       if (multi_step) then
          if (this%step .eq. 1) then
             allocate(this%previous1, source = this)
             this%previous1 = this
             if (this%step == this%previous1%step) then
                this%previous1 = this
             end if
          else if (this%step .eq. 2) then
             allocate(this%previous2, source = this)
             this%previous2 = this%previous1
             this%previous1 = this
             if (this%step == this%previous1%step) then
                this%previous1 = this
             end if
          else if (this%step .eq. 3) then
             allocate(this%previous3, source = this)
             this%previous3 = this%previous2
             this%previous2 = this%previous1
             this%previous1 = this
             if (this%step == this%previous1%step) then
                this%previous1 = this
             end if
          else
             if (this%step == this%previous1%step) then
                this%previous1 = this
             else
                this%previous3 = this%previous2
                this%previous2 = this%previous1
                this%previous1 = this
             end if
          end if
       end if
       this_half = this + this%t(this%state)*(0.5*dt)
       this      = this + this_half*dt
    class default
       stop 'integrate: unsupported class'
    end select
  end subroutine integrate

end module RK2M
