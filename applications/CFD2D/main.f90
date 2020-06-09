program main

  use UtilitiesM
  use DebuggerM
  use DataInputM
  use CFDApplicationM
  use CFDStrategyM

  implicit none

  type(CFDApplicationDT)   :: application
  type(CFDStrategyDT)      :: cfdStrategy
  real(rkind)              :: start, finish

  call cpu_time(start)
  call initFEM2D(application)
  call cfdStrategy%buildStrategyAndSolve(application%model)
  call cpu_time(finish)
  call debugLog('    Total time = ', (finish-start))
  print'(A,E14.7)', 'Total time = ', (finish-start)
 
end program main
