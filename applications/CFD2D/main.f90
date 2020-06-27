program main

  use UtilitiesM
  use DebuggerM
  use DataInputM
  use CFDApplicationM
  use CFDStrategyM
  use SolvingStrategyM

  implicit none

  type(CFDApplicationDT)  :: application
  type(CFDStrategyDT)     :: strategy
  type(SolvingStrategyDT) :: solvingStrategy
  real(rkind)             :: start, finish

  call cpu_time(start)
  call initFEM2D(application)
  solvingStrategy = InitSolvingStrategy(strategy, application)
  call solvingStrategy%useStrategy()
  call cpu_time(finish)
  call debugLog('    Total time = ', (finish-start))
  print'(A,E14.7)', 'Total time = ', (finish-start)
 
end program main
