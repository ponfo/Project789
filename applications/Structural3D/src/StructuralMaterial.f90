module StructuralMaterialM
  use UtilitiesM
  use PropertyM
  
  implicit none
  
  private
  public :: StructuralMaterialDT, structuralMaterial
  
  type :: StructuralMaterialDT
     real(rkind)                 :: young
     real(rkind)                 :: poissonCoef
     real(rkind)                 :: thermalCoef
     real(rkind), dimension(6,6) :: d
   contains
     procedure :: init
  end type StructuralMaterialDT
  
  interface structuralMaterial
     procedure :: constructor
  end interface structuralMaterial
  
contains
  
  type(StructuralMaterialDT) function constructor(young, poissonCoef, thermalCoef)
    implicit none
    real(rkind), intent(in) :: young
    real(rkind), intent(in) :: poissonCoef
    real(rkind), intent(in) :: thermalCoef
    call constructor%init(young, poissonCoef, thermalCoef)
  end function constructor

  subroutine init(this, young, poissonCoef, thermalCoef)
    implicit none
    class(StructuralMaterialDT), intent(inout) :: this
    real(rkind)                , intent(in)    :: young
    real(rkind)                , intent(in)    :: poissonCoef
    real(rkind)                , intent(in)    :: thermalCoef
    real(rkind)                                :: factor
    this%young       = young
    this%poissonCoef = poissonCoef
    this%thermalCoef = thermalCoef
    !Deformación plana:
    factor = young/((1+poissonCoef)*(1-2*poissonCoef))
    this%d(1,1) = factor*(1-poissonCoef)
    this%d(2,2) = factor*(1-poissonCoef)
    this%d(3,3) = factor*(1-poissonCoef)
    this%d(1,2) = factor*poissonCoef
    this%d(1,3) = factor*poissonCoef
    this%d(1,2) = factor*poissonCoef
    this%d(2,3) = factor*poissonCoef
    this%d(3,1) = factor*poissonCoef
    this%d(3,2) = factor*poissonCoef
    this%d(4,4) = factor*(1-2*poissonCoef)/2._rkind
    this%d(5,5) = factor*(1-2*poissonCoef)/2._rkind
    this%d(6,6) = factor*(1-2*poissonCoef)/2._rkind
  end subroutine init
  
end module StructuralMaterialM
