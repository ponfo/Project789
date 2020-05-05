###########################################################
            PROGRAM STRUCTURAL 2D ANALYSIS
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
*Set Cond Fix_Displacement_X_On_Line *nodes
*Set var a = condnumentities
*Set Cond Fix_Displacement_X_On_Point *nodes
*Set var b = condnumentities
*Set var c = operation(a+b)
*#---------------------------------------------------------
Fix_Displacement_X_Number..............: *c
*#---------------------------------------------------------
*Set Cond Fix_Displacement_Y_On_Line *nodes
*Set var a = condnumentities
*Set Cond Fix_Displacement_Y_On_Point *nodes
*Set var b = condnumentities
*Set var c = operation(a+b)
*#---------------------------------------------------------
Fix_Displacement_Y_Number..............: *c
*#---------------------------------------------------------
Load_Number_On_Points..................: *Gendata(Load_Number_On_Points,int)
*#---------------------------------------------------------
Load_Number_On_Surfaces................: *Gendata(Load_Number_On_Surfaces,int)
*#---------------------------------------------------------
*Set Cond Loads_On_Points *nodes
Points_With_Point_Load.................: *condnumentities
*#---------------------------------------------------------
*Set Cond Loads_On_Surfaces *elems
Surfaces_With_Surface_Load.............: *condnumentities
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

Material | Thermal expansion | Young's modulus | Poisson's ratio |  Area  |  Thickness
------------------------------------------------------------------------------------------
*loop materials
*format "%5i%10.4e%10.4e%10.4e%10.4e%10.4e"
*matnum          *matprop(Thermal_Expansion)         *matprop(Young_Modulus)       *matprop(Poisson's_Ratio)    *matprop(Area)  *matprop(Thickness)
*end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Element List:

      Element  |      Type      |  Material  |   Nodes  |      Conectivities
-----------------------------------------------------------------------------------
*Set Cond Loads_On_Surfaces *elems
*loop elems
*format "%10i%10i%9i%9i%9i"
*elemsnum         *ElemstypeName  *elemsmat  *ElemsNnode  *elemsconec
*end elems

################### Loads On Points ######################

Conditions List:

Load   |  Function
----------------------
*Set var j = 0
*Set cond Loads_On_Points *nodes
*for(i=1;i<=Gendata(Load_Number_On_Points,int);i=i+1))
*set var j=operation(j+1)
*loop nodes *OnlyInCond
*if(i==cond(Load_Number_On_Points,real))
*format "%5i%10s"
*j *cond(LoadXOP) *cond(LoadYOP)
*break
*endif
*end
*end for
*Set cond Loads_On_Surfaces
*for(i=1;i<=Gendata(Load_Number_On_Surfaces,int);i=i+1))
*set var j=operation(j+1)
*loop elems
*if(i==cond(Load_Number_On_Surfaces,real))
*format "%5i%10s"
*j *cond(LoadXOS) *cond(LoadYOS)
*break
*endif
*end
*end for

###################### Point Loads ######################

Conditions List:

  Node   |  Load
--------------------------
*Set Cond Loads_On_Points *nodes
*loop nodes *OnlyInCond
*format "%5i%5i"
*NodesNum      *cond(Load_Number_On_Points) 
*end

###################### Surface Loads #####################

Conditions List:

  Element   |  Load
--------------------------
*Set Cond Loads_On_Surfaces *elems
*loop elems *OnlyInCond
*format "%8i%8i"
*elemsnum  *cond(Load_Number_On_Surfaces)
*end

####################### Displacements X ######################

Conditions List:

  Node    |    Displacement X
--------------------------------
*Set Cond Fix_Displacement_X_On_Line *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Displacement_X) 
*end
*Set Cond Fix_Displacement_X_On_Point *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Displacement_X) 
*end

####################### Displacements Y ######################

Conditions List:

  Node    |    Displacement Y
--------------------------------
*Set Cond Fix_Displacement_Y_On_Line *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Displacement_Y) 
*end
*Set Cond Fix_Displacement_Y_On_Point *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Displacement_Y) 
*end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
