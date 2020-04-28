module ModelM

  use MeshM

  use SourceM

  use PropertyM

  implicit none

  private
  public :: ModelDT

  type ModelDT
     type(MeshDT) :: mesh
   contains
     
  end type ModelDT

contains
  
end module ModelM
