module FortranParserPtrM
  use UtilitiesM
  use FortranParser

  implicit none

  private
  public :: FortranParserPtrDT

  type :: FortranParserPtrDT
     class(EquationParser), pointer :: ptr
  end type FortranParserPtrDT

end module FortranParserPtrM
