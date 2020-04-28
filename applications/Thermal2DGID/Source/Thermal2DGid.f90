module thermal2DGidM

  use ThermalIODataM

  use ThermalModelM

  use ThermalStrategyM

  use SolvingStrategyM
  
  implicit none

  private
  public :: Thermal2DGidDT

  type Thermal2DGidDT
     type(ThermalIODataDT)   :: ioData
     type(ThermalModelDT)    :: model
     type(ThermalStrategyDT) :: strategy
   contains
     procedure :: initModel
     procedure :: initStategy
     procedure :: solve
     procedure :: writeOutputData
  end type Thermal2DGidDT

contains

  subroutine initModel(this)
    implicit none
    class(Thermal2DGidDT), intent(inout) :: this
    call this%ioData%init(this%model)
  end subroutine initModel

  subroutine initStrategy(this)
    implicit none
    class(Thermal2DGidDT)   , intent(inout) :: this
    class(SolvingStrategyDT), allocatable   :: userStrategy
    userStrategy = SetSolvingStrategy(this%strategy)
    write(*,*) '*** Solving Strategy Allocated ***'
  end subroutine initStrategy

  subroutine solve(this)
    implicit none
    class(Thermal2DGidDT), intent(inout) :: this
    call this%strategy%buildStrategyAndSolve(this%model)
  end subroutine solve
  
  subroutine writeOutputData(this)
    implicit none
    class(Thermal2DGidDT), intent(inout) :: this
    call this%ioData%writeOutputData(this%model)
  end subroutine writeOutputData
  
end module thermal2DGidM
