  E  µ   k820309              19.1        Þ5Å^                                                                                                          
       src/SolvingStrategy/Poisson2D.f90 POISSON2DM              POISSON2DDT gen@SETPOISSON2D                                                     
                            @                              
                            @                              
                            @                              
                                                              u #CONSTRUCTOR    &         @   @                                                     #INITIAL_STATE    #STIFFNESS    #RHS 	   #LUMPEDMASSINVERSE 
   #THIS_STRATEGY    #STEP    #POISSON2DDT              
                                                    
              &                                                     
                                                    #SPARSE              
                                 	                   
              &                                                     
                                  
                  #SPARSE              
  @                                                  #NEWSCHEMEDT              
 @                                                       @  @               @             '0                   #NEWPROCESSDT    #QUADRATURE    #STATE    #VALUES    #STEP    #INTEGRATE    #SET_QUADRATURE !   #GET_QUADRATURE %   #T (   #ADD ,   #MULTIPLY 0   #ASSIGN 4                 $                                                          #NEWPROCESSDT                   @  @                              '                      #USEPROCESS    1         À    $                                             #NEWPROCESS_PROCEDURE    #         @     @                                	               #THIS              
                                                    #NEWPROCESSDT                 $                                                         #NEWSCHEMEDT                  @  @                               '                      #INTEGRATE    1         À    $                      $                        #INTEGRATOR_INTERFACE    #         @     @                                	               #THIS    #DT              
                                                    #NEWPROCESSDT              
                                     
                 $                                                          
            &                                                       $                                         È                 
            &                   &                                                         $                                  (            1         À    $                                             #INTEGRATE    #         @     @                                                #MODEL    #DT                                                    0              #INTEGRANDDT              
                                       
      1         À    $                          !                  #SET_QUADRATURE "   #         @     @                           "                    #THIS #   #S $             
                                #     0              #INTEGRANDDT              
                                 $                    #NEWSCHEMEDT    1         À    $                         %                  #GET_QUADRATURE &   &        @    @                         &                            #THIS '   #NEWSCHEMEDT              
                                 '     0             #INTEGRANDDT    1         À    $                          (             	     #TIME_DERIVATIVE )   &        @    @                         )     0                     #THIS *   #DOF +   #INTEGRANDDT              
                                *     0             #INTEGRANDDT              
                                +                   
              &                                           1         À    $                          ,             
     #SYMMETRIC_OPERATOR -   &        @    @                         -     0                     #LHS .   #RHS /   #INTEGRANDDT              
                                .     0             #INTEGRANDDT              
                                /     0             #INTEGRANDDT    1         À    $                          0                  #ASYMMETRIC_OPERATOR 1   &        @    @                         1     0                     #LHS 2   #RHS 3   #INTEGRANDDT              
                                2     0             #INTEGRANDDT              
                                3     
      1         À    $                           4                  #SYMMETRIC_ASSIGNMENT 5   #         @     @                           5     	               #LHS 6   #RHS 7             
                               6     0              #INTEGRANDDT              
                                7     0             #INTEGRANDDT    3                                                           #INTEGRANDDT    #ADD ,   3                                                           #INTEGRANDDT    #MULTIPLY 0   3                                                          |  #INTEGRANDDT    #ASSIGN 4                     @              À                '                   #A 8   #AI 9   #AJ :   #ROWCOUNTER ;   #N <   #NNZ =   #COUNTER >   #TRIPLET ?   #ISCRSDONE D   #INIT E   #UPDATE J   #APPEND O   #MAKECRS U   #APPENDPOSTCRS Y   #CHANGE _   #SETDIRICHLET e   #GET i   #GETNNZ n   #GETN q   #GETA t   #GETAI w   #GETAJ z   #PRINTVALUE }   #PRINTNONZEROS    #PRINTALL    #DELETEROWANDCOL    #FREE    #HANDLEDUPLICATES                D                             8                              
            &                                                      D                             9            H                             &                                                      D                             :                                         &                                                       D                             ;            Ø                             &                                                         D                             <                                D                             =     $                          D                             >     (                          D                              ?     Ø       0             #TRIPLET @                 @  @              D           @     'Ø                    #A A   #ROW B   #COL C                                            A                              
            &                                                                                    B            H                             &                                                                                    C                                         &                                                         D                              D           	      1         À    $                            E             
     #INIT F   #         @     @                            F                    #THIS G   #NNZ H   #ROWS I             
                                G                   #SPARSE              
                                 H                     
                                 I           1         À    $                            J                  #UPDATE K   #         @     @                            K                    #THIS L   #NNZ M   #ROWS N             
                                L                   #SPARSE              
                                 M                     
                                 N           1         À    $                            O                  #APPEND P   #         @     @                            P                    #THIS Q   #VAL R   #ROW S   #COL T             
                                Q                   #SPARSE              
                                 R     
                
                                 S                     
                                 T           1         À    $                            U                  #MAKECRS V   #         @     @                            V                    #THIS W   #SORTROWS X             
                                W                   #SPARSE              
                                 X           1         À    $                            Y                  #APPENDPOSTCRS Z   #         @     @                            Z                    #THIS [   #VAL \   #ROW ]   #COL ^             
                                [                   #SPARSE              
                                 \     
                
                                 ]                     
                                 ^           1         À    $                            _                  #CHANGE `   #         @     @                            `                    #THIS a   #VAL b   #ROW c   #COL d             
                                a                   #SPARSE              
                                 b     
                
                                 c                     
                                 d           1         À    $                            e                  #SETDIRICHLET f   #         @     @                            f                    #THIS g   #ROW h             
                                g                   #SPARSE              
                                 h           1         À    $                           i                  #GET j   %         @   @                           j                    
       #THIS k   #I l   #J m             
                                k                   #SPARSE              
                                 l                     
                                 m           1         À    $                           n              	    #GETNNZ o   %         @   @                           o                           #THIS p             
                                p                   #SPARSE    1         À    $                           q              
    #GETN r   %         @   @                           r                           #THIS s             
                                s                   #SPARSE    1         À    $                           t                  #GETA u   (        `   @                           u                                    
    #THIS v   p          5 8 O#SPARSE     p        U     =       5 8 O#SPARSE     p        U     =                               
                                v                   #SPARSE    1         À    $                           w                  #GETAI x   (        `   @                           x                                        #THIS y   p           5 8 O#SPARSE     p        U     <   n                                         1     5 8 O#SPARSE     p        U     <   n                                          1                               
                                y                   #SPARSE    1         À    $                           z                  #GETAJ {   (        `   @                           {                                        #THIS |   p          5 8 O#SPARSE     p        U     =       5 8 O#SPARSE     p        U     =                               
                                |                   #SPARSE    1         À    $                            }                  #PRINTVALUE ~   #         @     @                            ~                    #THIS    #I    #J    #FILENAME              
                                                   #SPARSE              
                                                      
                                                      
                                                   1 1         À    $                                              #PRINTNONZEROS    #         @     @                                                #THIS    #FILENAME              
                                                   #SPARSE              
                                                   1 1         À    $                                              #PRINTALL    #         @     @                                                #THIS    #FILENAME              
                                                   #SPARSE              
                                                   1 1         À    $                                              #DELETEROWANDCOL    #         @     @                                                #THIS    #ROW    #COL              
                                                   #SPARSE              
                                                      
                                            1         À    $                                              #FREE    #         @     @                                                #THIS              
                                                   #SPARSE    1         À    D                                              #HANDLEDUPLICATES    #         @     @                                                #THIS              
                                                   #SPARSE                     @@               À              '             
      #INTEGRANDDT    #STIFFNESS    #LUMPEDMASSINVERSE    #RHS    #T    #ADD    #MULTIPLY ¢   #ASSIGN ¦   #GETSTATE ª   #USEPROCESS ­                 $                                   0                     #INTEGRANDDT                  $                                         0             #SPARSE                  $                                         @             #SPARSE                $                                         P                
            &                                           1         À    $                                            #DPOISSON2D_DT    &        @    @                                0                     #THIS    #DOF    #INTEGRANDDT              
                                                   #POISSON2DDT              
                                                    
              &                                           1         À    $                                            #ADD_POISSON2D    &        @    @                                0                     #LHS     #RHS ¡   #INTEGRANDDT              
                                                    #POISSON2DDT              
                                 ¡     0             #INTEGRANDDT    1         À    $                          ¢                  #MULTIPLY_POISSON2D £   &        @    @                           £     0                     #LHS ¤   #RHS ¥   #INTEGRANDDT              
                                 ¤                  #POISSON2DDT              
                                 ¥     
      1         À    $                           ¦                  #ASSIGN_POISSON2D §   #         @     @                             §                    #LHS ¨   #RHS ©             
D @                              ¨                   #POISSON2DDT              
                                 ©     0             #INTEGRANDDT    1         À    $                           ª             	 	    #GETSTATE «   (        D    @                            «                                   
    #THIS ¬             &                                                     
                                 ¬                  #POISSON2DDT    1         À    $                           ­             
     #PROCESS ®   #         @     @                             ®                    #THIS ¯             
                                ¯                   #POISSON2DDT           5      fn#fn     Õ   -   b   uapp(POISSON2DM      @   J  UTILITIESM    B  @   J  SPARSEKIT      @   J  SCHEMEM    Â  @   J  INTEGRANDM !     Q       gen@SETPOISSON2D    S  À      CONSTRUCTOR *        a   CONSTRUCTOR%INITIAL_STATE &     T   a   CONSTRUCTOR%STIFFNESS     ó     a   CONSTRUCTOR%RHS .     T   a   CONSTRUCTOR%LUMPEDMASSINVERSE *   Ó  Y   a   CONSTRUCTOR%THIS_STRATEGY !   ,  @   a   CONSTRUCTOR%STEP '   l  ô      INTEGRANDDT+INTEGRANDM 4   `  b   a   INTEGRANDDT%NEWPROCESSDT+INTEGRANDM &   Â  `      NEWPROCESSDT+PROCESSM 1   "  b   a   NEWPROCESSDT%USEPROCESS+PROCESSM .     R      NEWPROCESS_PROCEDURE+PROCESSM 3   Ö  Z   a   NEWPROCESS_PROCEDURE%THIS+PROCESSM 2   0  a   a   INTEGRANDDT%QUADRATURE+INTEGRANDM $     _      NEWSCHEMEDT+SCHEMEM .   ð  b   a   NEWSCHEMEDT%INTEGRATE+SCHEMEM -   R	  Z      INTEGRATOR_INTERFACE+SCHEMEM 2   ¬	  Z   a   INTEGRATOR_INTERFACE%THIS+SCHEMEM 0   
  @   a   INTEGRATOR_INTERFACE%DT+SCHEMEM -   F
     a   INTEGRANDDT%STATE+INTEGRANDM .   Ú
  ¬   a   INTEGRANDDT%VALUES+INTEGRANDM ,     H   a   INTEGRANDDT%STEP+INTEGRANDM 1   Î  W   a   INTEGRANDDT%INTEGRATE+INTEGRANDM %   %  [      INTEGRATE+INTEGRANDM +     Y   a   INTEGRATE%MODEL+INTEGRANDM (   Ù  @   a   INTEGRATE%DT+INTEGRANDM 6     \   a   INTEGRANDDT%SET_QUADRATURE+INTEGRANDM *   u  Y      SET_QUADRATURE+INTEGRANDM /   Î  Y   a   SET_QUADRATURE%THIS+INTEGRANDM ,   '  Y   a   SET_QUADRATURE%S+INTEGRANDM 6     \   a   INTEGRANDDT%GET_QUADRATURE+INTEGRANDM *   Ü  k      GET_QUADRATURE+INTEGRANDM /   G  Y   a   GET_QUADRATURE%THIS+INTEGRANDM )      ]   a   INTEGRANDDT%T+INTEGRANDM +   ý  t      TIME_DERIVATIVE+INTEGRANDM 0   q  Y   a   TIME_DERIVATIVE%THIS+INTEGRANDM /   Ê     a   TIME_DERIVATIVE%DOF+INTEGRANDM +   V  `   a   INTEGRANDDT%ADD+INTEGRANDM .   ¶  s      SYMMETRIC_OPERATOR+INTEGRANDM 2   )  Y   a   SYMMETRIC_OPERATOR%LHS+INTEGRANDM 2     Y   a   SYMMETRIC_OPERATOR%RHS+INTEGRANDM 0   Û  a   a   INTEGRANDDT%MULTIPLY+INTEGRANDM /   <  s      ASYMMETRIC_OPERATOR+INTEGRANDM 3   ¯  Y   a   ASYMMETRIC_OPERATOR%LHS+INTEGRANDM 3     @   a   ASYMMETRIC_OPERATOR%RHS+INTEGRANDM .   H  b   a   INTEGRANDDT%ASSIGN+INTEGRANDM 0   ª  Z      SYMMETRIC_ASSIGNMENT+INTEGRANDM 4     Y   a   SYMMETRIC_ASSIGNMENT%LHS+INTEGRANDM 4   ]  Y   a   SYMMETRIC_ASSIGNMENT%RHS+INTEGRANDM     ¶  Z   p   i@+INTEGRANDM       _   p   i@+INTEGRANDM    o  ]   p   i@| !   Ì  µ      SPARSE+SPARSEKIT %        %   SPARSE%A+SPARSEKIT=A '        %   SPARSE%AI+SPARSEKIT=AI '   ©     %   SPARSE%AJ+SPARSEKIT=AJ 7   =     %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   Ñ  H   %   SPARSE%N+SPARSEKIT=N )     H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   a  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   ©  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "     i      TRIPLET+SPARSEKIT $   o     a   TRIPLET%A+SPARSEKIT &        a   TRIPLET%ROW+SPARSEKIT &        a   TRIPLET%COL+SPARSEKIT 5   +  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   s  R   a   SPARSE%INIT+SPARSEKIT    Å  e      INIT+SPARSEKIT $   *  T   a   INIT%THIS+SPARSEKIT #   ~  @   a   INIT%NNZ+SPARSEKIT $   ¾  @   a   INIT%ROWS+SPARSEKIT (   þ  T   a   SPARSE%UPDATE+SPARSEKIT !   R   e      UPDATE+SPARSEKIT &   ·   T   a   UPDATE%THIS+SPARSEKIT %   !  @   a   UPDATE%NNZ+SPARSEKIT &   K!  @   a   UPDATE%ROWS+SPARSEKIT (   !  T   a   SPARSE%APPEND+SPARSEKIT !   ß!  m      APPEND+SPARSEKIT &   L"  T   a   APPEND%THIS+SPARSEKIT %    "  @   a   APPEND%VAL+SPARSEKIT %   à"  @   a   APPEND%ROW+SPARSEKIT %    #  @   a   APPEND%COL+SPARSEKIT )   `#  U   a   SPARSE%MAKECRS+SPARSEKIT "   µ#  `      MAKECRS+SPARSEKIT '   $  T   a   MAKECRS%THIS+SPARSEKIT +   i$  @   a   MAKECRS%SORTROWS+SPARSEKIT /   ©$  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   %  m      APPENDPOSTCRS+SPARSEKIT -   q%  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   Å%  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   &  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   E&  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   &  T   a   SPARSE%CHANGE+SPARSEKIT !   Ù&  m      CHANGE+SPARSEKIT &   F'  T   a   CHANGE%THIS+SPARSEKIT %   '  @   a   CHANGE%VAL+SPARSEKIT %   Ú'  @   a   CHANGE%ROW+SPARSEKIT %   (  @   a   CHANGE%COL+SPARSEKIT .   Z(  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   ´(  [      SETDIRICHLET+SPARSEKIT ,   )  T   a   SETDIRICHLET%THIS+SPARSEKIT +   c)  @   a   SETDIRICHLET%ROW+SPARSEKIT %   £)  Q   a   SPARSE%GET+SPARSEKIT    ô)  h      GET+SPARSEKIT #   \*  T   a   GET%THIS+SPARSEKIT     °*  @   a   GET%I+SPARSEKIT     ð*  @   a   GET%J+SPARSEKIT (   0+  T   a   SPARSE%GETNNZ+SPARSEKIT !   +  Z      GETNNZ+SPARSEKIT &   Þ+  T   a   GETNNZ%THIS+SPARSEKIT &   2,  R   a   SPARSE%GETN+SPARSEKIT    ,  Z      GETN+SPARSEKIT $   Þ,  T   a   GETN%THIS+SPARSEKIT &   2-  R   a   SPARSE%GETA+SPARSEKIT    -  ö      GETA+SPARSEKIT $   z.  T   a   GETA%THIS+SPARSEKIT '   Î.  S   a   SPARSE%GETAI+SPARSEKIT     !/  h     GETAI+SPARSEKIT %   0  T   a   GETAI%THIS+SPARSEKIT '   Ý0  S   a   SPARSE%GETAJ+SPARSEKIT     01  ö      GETAJ+SPARSEKIT %   &2  T   a   GETAJ%THIS+SPARSEKIT ,   z2  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   Ò2  n      PRINTVALUE+SPARSEKIT *   @3  T   a   PRINTVALUE%THIS+SPARSEKIT '   3  @   a   PRINTVALUE%I+SPARSEKIT '   Ô3  @   a   PRINTVALUE%J+SPARSEKIT .   4  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   `4  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   »4  `      PRINTNONZEROS+SPARSEKIT -   5  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   o5  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   »5  V   a   SPARSE%PRINTALL+SPARSEKIT #   6  `      PRINTALL+SPARSEKIT (   q6  T   a   PRINTALL%THIS+SPARSEKIT ,   Å6  L   a   PRINTALL%FILENAME+SPARSEKIT 1   7  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   n7  d      DELETEROWANDCOL+SPARSEKIT /   Ò7  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   &8  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   f8  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   ¦8  R   a   SPARSE%FREE+SPARSEKIT    ø8  R      FREE+SPARSEKIT $   J9  T   a   FREE%THIS+SPARSEKIT C   9  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   ü9  R      HANDLEDUPLICATES+SPARSEKIT 0   N:  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT    ¢:  Ø       POISSON2DDT (   z;  a   a   POISSON2DDT%INTEGRANDDT &   Û;  \   a   POISSON2DDT%STIFFNESS .   7<  \   a   POISSON2DDT%LUMPEDMASSINVERSE     <     a   POISSON2DDT%RHS    '=  [   a   POISSON2DDT%T    =  t      DPOISSON2D_DT #   ö=  Y   a   DPOISSON2D_DT%THIS "   O>     a   DPOISSON2D_DT%DOF     Û>  [   a   POISSON2DDT%ADD    6?  s      ADD_POISSON2D "   ©?  Y   a   ADD_POISSON2D%LHS "   @  Y   a   ADD_POISSON2D%RHS %   [@  `   a   POISSON2DDT%MULTIPLY #   »@  s      MULTIPLY_POISSON2D '   .A  Y   a   MULTIPLY_POISSON2D%LHS '   A  @   a   MULTIPLY_POISSON2D%RHS #   ÇA  ^   a   POISSON2DDT%ASSIGN !   %B  Z      ASSIGN_POISSON2D %   B  Y   a   ASSIGN_POISSON2D%LHS %   ØB  Y   a   ASSIGN_POISSON2D%RHS %   1C  V   a   POISSON2DDT%GETSTATE    C  ¦      GETSTATE    -D  Y   a   GETSTATE%THIS '   D  U   a   POISSON2DDT%USEPROCESS    ÛD  R      PROCESS    -E  Y   a   PROCESS%THIS 