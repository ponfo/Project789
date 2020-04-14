module UserAssembleStrategyM
  use UtilitiesM
  use DebuggerM
  use AssembleStrategyM

  implicit none

  private
  public :: UserAssembleStrategyDT, userAssembleStrategy

  type :: UserAssembleStrategyDT
     class(AssembleStrategyDT), allocatable :: assembler
   contains
     procedure, public :: init
     procedure, public :: changeStrategy
     procedure, public :: assembleProcedure
  end type UserAssembleStrategyDT

  interface userAssembleStrategy
     procedure :: constructor
  end interface userAssembleStrategy

contains

  type(UserAssembleStrategyDT) function constructor(assembleStrategy)
    implicit none
    class(AssembleStrategyDT), intent(in)    :: assembleStrategy
    call constructor%init(assembleStrategy)
  end function constructor

  subroutine init(this, assembleStrategy)
    implicit none
    class(UserAssembleStrategyDT), intent(inout) :: this
    class(AssembleStrategyDT)    , intent(in)    :: assembleStrategy
    allocate(this%assembler, source = assembleStrategy)
  end subroutine init

end module UserAssembleStrategyM

    

  
