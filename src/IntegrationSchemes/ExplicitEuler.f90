module ExplicitEulerM

  use UtilitiesM

  use ProcessM

  use BaseIntegrandM

  use IntegrandM

  implicit none

  private
  public :: ExplicitEulerDT

  type, extends(BaseIntegrandDT) :: ExplicitEulerDT
   contains
     procedure, nopass :: integrate
  end type ExplicitEulerDT

contains

  subroutine integrate(this, dt, multi_step) 
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
    logical            , intent(in)    :: multi_step
    select type (this)
    class is (IntegrandDT)
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
       this = this + this%t(this%state)*dt
    class default
       stop 'integrate: unsupported class.'
    end select
  end subroutine integrate

end module ExplicitEulerM
