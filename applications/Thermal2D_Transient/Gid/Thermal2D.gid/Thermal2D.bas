###########################################################
            PROGRAM THERMAL 2D ANALYSIS
###########################################################

##################### PROBLEM DATA ########################

Problem Name: *gendata(Problem_Name)

###################### Mesh Data ##########################
Elements_Number........................: *nelem
Nodes_Number...........................: *npoin
Are_Elements_Quadratic.................: *isQuadratic
*#---------------------------------------------------------
*set var j=0
*set var k=0
*loop elems
*#ElemsTypeName
*if(strcasecmp(ElemsTypeName(),"Triangle")==0)
*set var j=operation(j+1)
*endif
*if(strcasecmp(ElemsTypeName(),"Quadrilateral")==0)
*set var k=operation(k+1)
*endif
*end elems
*#---------------------------------------------------------
Triangular_Elements_Number.............: *j
Rectangular_Elements_Number............: *k
*#---------------------------------------------------------
Materials_Number.......................: *nmats
Gauss_Order............................: *GenData(Gauss_Order)
*#---------------------------------------------------------
*Set Cond Fix_Temperature_On_Lines *nodes
*Set var a = condnumentities
*Set Cond Fix_Temperature_On_Points *nodes
*Set var b = condnumentities
*Set var c = operation(a+b)
*#---------------------------------------------------------
Dirichlet_Conditions_Number............: *c
*#---------------------------------------------------------
*Set Cond Normal_Flux_On_Lines *elems
Normal_Flux_Line_Condition_Elements....: *condnumentities
*#--------------------------------------------------------- REVISAR ESTAS DOS CONDICIONES DE LINEAR, *ELEMS ESTÁ MAL
*Set Cond Convection_On_Lines *elems
Convection_Lines_Condition_Elements....: *condnumentities
*#---------------------------------------------------------
Source_Number_On_Points................: *Gendata(Source_Number_On_Points,int)
*#---------------------------------------------------------
Source_Number_On_Surfaces..............: *Gendata(Source_Number_On_Surfaces,int)
*#---------------------------------------------------------
*Set Cond Source_On_Points *nodes
Points_With_Point_Source...............: *condnumentities
*#---------------------------------------------------------
*Set Cond Source_On_Surfaces *elems
Surfaces_With_Surface_Source...........: *condnumentities
*#---------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coordinates:

  Node   |     X       |       Y        |
----------------------------------------- 
*set elems(all)
*loop nodes
*format "%5i%10.4e%10.4e"
*NodesNum       *NodesCoord(1,real)     *NodesCoord(2,real)
*end nodes

######################## Materials ########################

Materials List:

Material | Thermal conductivity X | Thermal conductivity Y
-----------------------------------------------------------
*loop materials
*format "%5i%10.4e%10.4e"
*matnum  *matprop(Thermal_Conductivity_X)  *matprop(Thermal_Conductivity_Y) *matprop(cp) *matprop(rho)
*end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Element List:

      Element  |      Type      |  Material  |   Nodes  |      Conectivities
-----------------------------------------------------------------------------------
*Set Cond Source_On_Surfaces *elems
*loop elems
*format "%10i%10i%9i%9i%9i"
*elemsnum         *ElemstypeName  *elemsmat  *ElemsNnode  *elemsconec
*end elems

################### Sources On Points ######################

Conditions List:

Source   |  Function
----------------------
*Set var j = 0
*Set cond Source_On_Points *nodes
*for(i=1;i<=Gendata(Source_Number_On_Points,int);i=i+1))
*set var j=operation(j+1)
*loop nodes *OnlyInCond
*if(i==cond(Source_Number_On_Points,real))
*format "%5i%10s"
*j *cond(SourceOP)
*break
*endif
*end
*end for
*Set cond Source_On_Surfaces
*for(i=1;i<=Gendata(Source_Number_On_Surfaces,int);i=i+1))
*set var j=operation(j+1)
*loop elems
*if(i==cond(Source_Number_On_Surfaces,real))
*format "%5i%10s"
*j *cond(SourceOS) 
*break
*endif
*end
*end for

###################### Point Sources ######################

Conditions List:

  Node   |  Source
--------------------------
*Set Cond Source_On_Points *nodes
*loop nodes *OnlyInCond
*format "%5i%5i"
*NodesNum      *cond(Source_Number_On_Points) 
*end

###################### Surface Sources #####################

Conditions List:

  Element   |  Source
--------------------------
*Set Cond Source_On_Surfaces *elems
*loop elems *OnlyInCond
*format "%8i%8i"
*elemsnum  *cond(Source_Number_On_Surfaces)
*end

####################### Temperatures ######################

Conditions List:

  Node    |    Temperature
------------------------------
*Set Cond Fix_Temperature_On_Lines *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Temperature) 
*end
*Set Cond Fix_Temperature_On_Points *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Temperature) 
*end

########################## Normal Fluxs On Lines ##########################

Conditions List:

 Element |       Nodes      |   Normal Flux
--------------------------------------------
*Set Cond Normal_Flux_On_Lines *elems *canrepeat
*loop elems *OnlyInCond
*format "%5i%7i%7i"
*elemsnum  *localnodes  *cond(Flux,real)
*end

########################### Convection On Lines ###########################

Conditions List:

 Element |       Nodes      |   Coeficient  |  Temperature
-------------------------------------------------------------
*Set Cond Convection_On_Lines *elems *canrepeat
*loop elems *OnlyInCond
*elemsnum   *localnodes  *cond(Coeficient)  *cond(Temperature)
*end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
