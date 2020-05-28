module AdamsB4M

  use UtilitiesM

  use ProcessM

  use SchemeM

  use IntegrandM

  implicit none

  private
  public :: AdamsB4DT

  type, extends(NewSchemeDT) :: AdamsB4DT
   contains
     procedure, nopass  :: integrate
  end type ADAMSB4DT

contains

  subroutine integrate(this,dt)
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
    class(IntegrandDT), allocatable    :: k1
    class(IntegrandDT), allocatable    :: k2
    class(IntegrandDT), allocatable    :: k3
    class(IntegrandDT), allocatable    :: k4
    select type (this)
    class is (IntegrandDT)
       allocate(k1, source = this)
       allocate(k2, source = this)
       allocate(k3, source = this)
       allocate(k4, source = this)
       if (this%step .eq. 0) then
          if (allocated(this%values)) deallocate(this%values)
          allocate(this%values(3,size(this%state)))
          this%values = 0.d0
          k1 = this%t(this%state)
          k2 = this%t(this%state) * 0.d0
          k3 = this%t(this%state) * 0.d0
          k4 = this%t(this%state) * 0.d0
          this%values(1,:) = this%state
       else if (this%step == 1) then
          k4 = this%t(this%state) * 0.d0
          k3 = this%t(this%state) * 0.d0
          k2 = this%t(this%values(1,:))
          k1 = this%t(this%state)
          this%values(2,:) = this%state
       else if (this%step == 2) then
          k4 = this%t(this%state) * 0.d0
          k3 = this%t(this%values(1,:))
          k2 = this%t(this%values(2,:))
          k1 = this%t(this%state)
          this%values(3,:) = this%state
       else if (this%step .ge. 3) then
          k4 = this%t(this%values(1,:))
          k3 = this%t(this%values(2,:))
          k2 = this%t(this%values(3,:))
          k1 = this%t(this%state)
          this%values(1,:) = this%values(2,:)
          this%values(2,:) = this%values(3,:)
          this%values(3,:) = this%state
       end if
       this = this + (k4 * (-9.0_rkind/24.0_rkind)&
            + k3 * (37.0_rkind/24.0_rkind)       &
            + k2 * (-59.0_rkind/24.0_rkind)      &
            + k1 * (55.0_rkind/24.0_rkind))*dt 
    class default
       stop 'integrate: unsupported class'
    end select
  end subroutine integrate
  
end module AdamsB4M
