program main
  
  use DataInputM
  use CFDApplicationM
  use CFDStrategyM
  use GIDDataOutputM

  implicit none

  type(CFDApplicationDT)   :: application
  type(CFDStrategyDT)      :: cfdStrategy

  call initFEM2D(application)
  call cfdStrategy%buildStrategyAndSolve(application%model)
 
end program main
