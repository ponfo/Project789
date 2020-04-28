program main
  use Thermal2DGidM
  implicit none
  type(Thermal2DGidDT) :: thermal2DGid
  call thermal2DGid%initModel()
  call thermal2DGid%initStrategy()
  call thermal2DGid%solve()
  call thermal2DGid%writeOutputData() 
end program main
