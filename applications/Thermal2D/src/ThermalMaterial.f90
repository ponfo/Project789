module ThermalMaterialM
  
  use UtilitiesM
  use PropertyM
  
  implicit none
  
  private
  public :: ThermalMaterialDT, thermalMaterial
  
  type :: ThermalMaterialDT
     real(rkind), dimension(2) :: conductivity
   contains
     procedure :: init
  end type ThermalMaterialDT
  
  interface thermalMaterial
     procedure :: constructor
  end interface thermalMaterial
  
contains
  
  type(ThermalMaterialDT) function constructor(kx,ky)
    implicit none
    real(rkind), intent(in) :: kx, ky
    call constructor%init(kx, ky)
  end function constructor
  
  subroutine init(this, kx, ky)
    implicit none
    real(rkind), intent(in) :: kx, ky
    class(ThermalMaterialDT), intent(inout) :: this
    this%conductivity = (/kx, ky/)
  end subroutine init
  
end module ThermalMaterialM
