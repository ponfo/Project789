module RK4M

  use UtilitiesM

  use ProcessM
  
  use SchemeM
  
  use IntegrandM
  
  implicit none
  
  private
  public :: RK4DT
  
  type, extends(NewSchemeDT) :: RK4DT
   contains
     procedure, nopass :: integrate
  end type RK4DT
  
contains
  
  subroutine integrate(this,dt)
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
    select type (this)
    class is (IntegrandDT)
       allocate(v1(size(this%state)))
       allocate(v2(size(this%state)))
       allocate(v3(size(this%state)))
       allocate(v4(size(this%state)))
       allocate(k1, source = this)
       allocate(k2, source = this)
       allocate(k3, source = this)
       allocate(k4, source = this)
       v1 = this%state
       k1 = k1%t(v1)
       v2 = this%state+0.5*dt*k1%state
       k2 = this%t(v2)
       v3 = this%state+0.5*dt*k2%state
       k3 = this%t(v3)
       v4 = this%state+dt*k3%state
       k4 = this%t(v4)
       this=this&
            +(k1+k2*2.d0+k3*2.d0+k4)*(1/6.d0)*dt
       deallocate(k1,k2,k3,k4,v1,v2,v3,v4)
    class default
       stop 'integrate: unsupported class'
    end select
  end subroutine integrate
  
end module RK4M
