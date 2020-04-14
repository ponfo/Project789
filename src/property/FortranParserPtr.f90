module FortranParserPtrM
  use UtilitiesM
  use FortranParserM

  implicit none

  private
  public :: FortranParserPtrDT

  type :: FortranParserPtrDT
     class(FortranParser), pointer :: ptr
  end type FortranParserPtrDT

end module FortranParserPtrM
