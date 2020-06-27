module RK4M

  use UtilitiesM

  use ProcessM

  use BaseIntegrandM

  use IntegrandM

  implicit none

  private
  public :: RK4DT

  type, extends(BaseIntegrandDT) :: RK4DT
   contains
     procedure, nopass :: integrate
  end type RK4DT

contains

  subroutine integrate(this,dt, multi_step)
    class(NewProcessDT), intent(inout)     :: this
    real(rkind), intent(in)                :: dt
    real(rkind), dimension(:), allocatable :: v1
    real(rkind), dimension(:), allocatable :: v2
    real(rkind), dimension(:), allocatable :: v3
    real(rkind), dimension(:), allocatable :: v4
    class(IntegrandDT), allocatable        :: k1
    class(IntegrandDT), allocatable        :: k2
    class(IntegrandDT), allocatable        :: k3
    class(IntegrandDT), allocatable        :: k4
    logical                   , intent(in) :: multi_step
    select type (this)
    class is (IntegrandDT)
       allocate(k1, source = this)
       allocate(k2, source = this)
       allocate(k3, source = this)
       allocate(k4, source = this)
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
       v1 = this%state
       k1 = this%t(v1)
       v2 = this%state+0.5*dt*k1%state
       k2 = this%t(v2)
       v3 = this%state+0.5*dt*k2%state
       k3 = this%t(v3)
       v4 = this%state+dt*k3%state
       k4 = this%t(v4)
       this=this&
            +(k1+k2*2.d0+k3*2.d0+k4)*(1/6.d0)*dt
       deallocate(k1,k2,k3,k4)
    class default
       stop 'integrate: unsupported class'
    end select
  end subroutine integrate

end module RK4M
