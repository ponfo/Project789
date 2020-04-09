program thermal2DGid
  use DataInputM
  use SolverM
  use heatFluxM
  use DataOutputM
  use finishM
  implicit none
  type(IODataTYPE) :: io
  call initThermal(io)
  call solverThermal(io)
  call heatflux(io)
  call dataOutput(io)
  call finishProgram() 
end program thermal2DGid
