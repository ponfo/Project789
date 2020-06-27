module AdamsB4M

  use UtilitiesM

  use ProcessM

  use BaseIntegrandM

  use IntegrandM

  implicit none

  private
  public :: AdamsB4DT

  type, extends(BaseIntegrandDT) :: AdamsB4DT
   contains
     procedure, nopass  :: integrate
  end type ADAMSB4DT

contains

  subroutine integrate(this,dt,multi_step)
    class(NewProcessDT), intent(inout) :: this
    class(IntegrandDT) , pointer       :: temp
    real(rkind)        , intent(in)    :: dt
    logical            , intent(in)    :: multi_step
    select type (this)
    class is (IntegrandDT)
       if (this%step .eq. 1) then
          allocate(temp, source = this)
          temp = this
          this = this + (this%t(this%state) * (55.0_rkind/24.0_rkind))*dt
          allocate(this%previous1, source = temp)
          this%previous1 = temp
          deallocate(temp)
       else if (this%step .eq. 2) then
          allocate(temp, source = this)
          temp = this
          this = this + (this%previous1%t(this%previous1%state) * (-59.0_rkind/24.0_rkind)&
               + this%t(this%state)                            * (55.0_rkind/24.0_rkind))*dt
          allocate(this%previous2, source = this)
          this%previous2 = this%previous1
          this%previous1 = temp
          deallocate(temp)
       else if (this%step .eq. 3) then
          allocate(temp, source = this)
          temp = this
          this = this + (this%previous2%t(this%previous2%state) * (37.0_rkind/24.0_rkind) &
               + this%previous1%t(this%previous1%state)         * (-59.0_rkind/24.0_rkind)&
               + this%t(this%state)                            * (55.0_rkind/24.0_rkind))*dt
          allocate(this%previous3, source = this)
          this%previous3 = this%previous2
          this%previous2 = this%previous1
          this%previous1 = temp
          deallocate(temp)
       else
          allocate(temp, source = this)
          temp = this
          this = this + (this%previous3%t(this%previous3%state) * (-9.0_rkind/24.0_rkind) &
               + this%previous2%t(this%previous2%state)         * (37.0_rkind/24.0_rkind) &
               + this%previous1%t(this%previous1%state)         * (-59.0_rkind/24.0_rkind)&
               + this%t(this%state)                            * (55.0_rkind/24.0_rkind))*dt
          this%previous3 = this%previous2
          this%previous2 = this%previous1
          this%previous1 = temp
          deallocate(temp)
       end if
    class default
       stop 'integrate: unsupported class'
    end select
  end subroutine integrate

end module AdamsB4M
