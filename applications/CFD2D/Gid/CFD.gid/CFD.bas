Programa de Navier-Stokes (2D) 
FSAFE      U       V      MACH   TEMP      RHO    PRESS
*GenData(5)    *GenData(6)       0.0     *GenData(7)     *GenData(8)    *GenData(9)     0.0
FMU     FGX     FGY    Q
*GenData(10)   *GenData(11)   *GenData(12) *GenData(13) 
FK      FR        FCv      Gama
*GenData(14)   *GenData(15)   *GenData(16)   *GenData(17)   *GenData(18)
Cte
*GenData(19)
*realformat "%16.6f"
*intformat "%7i"
*Set  Cond  FIX_VELOCITY *nodes 
*set  var NFIXVI(int)=CondNumEntities(int)
*Set  Cond NORMAL_VELOCITY *elems *Canrepeat
*set  var NELNORM(int)=CondNumEntities(int) 
*Set  Cond  FIX_DENSITY *nodes
*set  var NFIXRHO(int)=CondNumEntities(int) 
*Set  Cond FIX_TEMPERATURE *nodes
*set  var NFIXT(int)=CondNumEntities(int)
*npoin   *nelem *NFIXRHO *NFIXVI  *NELNORM  *NFIXT    
NODES
*loop nodes
  *NodesNum *NodesCoord
*end
ELEMENTS
*loop elems
  *elemsnum *elemsConec
*end 
Fix Density 
*Set Cond FIX_DENSITY *nodes
*loop nodes *OnlyInCond
     *nodesNum *Cond(1)
*end nodes
Inflow Velocity
*Set Cond FIX_VELOCITY *nodes
*loop nodes *OnlyInCond
     *nodesNum *Cond(1) *Cond(1)
*end nodes
Normal Velocity
*Set Cond NORMAL_VELOCITY *elems *Canrepeat
*loop elems *OnlyInCond
      *globalnodes
*end elems
Fix Temperature
*Set Cond FIX_TEMPERATURE *nodes
*loop nodes *OnlyInCond
     *nodesNum *Cond(1)
*end nodes

