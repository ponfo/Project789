  ¬1     k820309              19.1        Ub¹^                                                                                                          
       src/solvers/Linear/Preconditioner/UsePreconditioner.f90 USEPRECONDITIONERM              USEPRECONDITIONERDT gen@SETPRECONDITIONER                                                     
                            @                              
                            @                              
                                                              u #CONSTRUCTOR    &         @   @                                                      #PRECONDITIONER    #USEPRECONDITIONERDT              
D @                                                   #PRECONDITIONERDT                  @  @                               '                      #USEPRECONDITIONER    1         À    $                                            #PRECONDITIONER_PROCEDURE 	   #         @     @                          	     	               #THIS 
   #VECTOR    #MATRIX    #SOLUTION    #ARG              
                               
                     #PRECONDITIONERDT              
                                                  
               &                                                     
                                                  #SPARSE              
                                                  
               &                                                     
                                                                 &                                                             @               À                '                   #A    #AI    #AJ    #ROWCOUNTER    #N    #NNZ    #COUNTER    #TRIPLET    #ISCRSDONE    #INIT    #UPDATE "   #APPEND '   #MAKECRS -   #APPENDPOSTCRS 1   #CHANGE 7   #SETDIRICHLET =   #GET A   #GETNNZ F   #GETN I   #GETA L   #GETAI O   #GETAJ R   #PRINTVALUE U   #PRINTNONZEROS [   #PRINTALL _   #DELETEROWANDCOL c   #FREE h   #HANDLEDUPLICATES k               D                                                           
            &                                                      D                                         H                             &                                                      D                                                                      &                                                       D                                         Ø                             &                                                         D                                                             D                                  $                          D                                  (                          D                                   Ø       0             #TRIPLET                  @  @              D                'Ø                    #A    #ROW    #COL                                                                           
            &                                                                                                H                             &                                                                                                                             &                                                         D                                         	      1         À    $                                         
     #INIT    #         @     @                                                #THIS    #NNZ     #ROWS !             
                                                   #SPARSE              
                                                       
                                 !           1         À    $                            "                  #UPDATE #   #         @     @                            #                    #THIS $   #NNZ %   #ROWS &             
                                $                   #SPARSE              
                                 %                     
                                 &           1         À    $                            '                  #APPEND (   #         @     @                            (                    #THIS )   #VAL *   #ROW +   #COL ,             
                                )                   #SPARSE              
                                 *     
                
                                 +                     
                                 ,           1         À    $                            -                  #MAKECRS .   #         @     @                            .                    #THIS /   #SORTROWS 0             
                                /                   #SPARSE              
                                 0           1         À    $                            1                  #APPENDPOSTCRS 2   #         @     @                            2                    #THIS 3   #VAL 4   #ROW 5   #COL 6             
                                3                   #SPARSE              
                                 4     
                
                                 5                     
                                 6           1         À    $                            7                  #CHANGE 8   #         @     @                            8                    #THIS 9   #VAL :   #ROW ;   #COL <             
                                9                   #SPARSE              
                                 :     
                
                                 ;                     
                                 <           1         À    $                            =                  #SETDIRICHLET >   #         @     @                            >                    #THIS ?   #ROW @             
                                ?                   #SPARSE              
                                 @           1         À    $                           A                  #GET B   %         @   @                           B                    
       #THIS C   #I D   #J E             
                                C                   #SPARSE              
                                 D                     
                                 E           1         À    $                           F              	    #GETNNZ G   %         @   @                           G                           #THIS H             
                                H                   #SPARSE    1         À    $                           I              
    #GETN J   %         @   @                           J                           #THIS K             
                                K                   #SPARSE    1         À    $                           L                  #GETA M   (        `   @                           M                                    
    #THIS N   p          5 8 O#SPARSE     p        U            5 8 O#SPARSE     p        U                                    
                                N                   #SPARSE    1         À    $                           O                  #GETAI P   (        `   @                           P                                        #THIS Q   p           5 8 O#SPARSE     p        U        n                                         1     5 8 O#SPARSE     p        U        n                                          1                               
                                Q                   #SPARSE    1         À    $                           R                  #GETAJ S   (        `   @                           S                                        #THIS T   p          5 8 O#SPARSE     p        U            5 8 O#SPARSE     p        U                                    
                                T                   #SPARSE    1         À    $                            U                  #PRINTVALUE V   #         @     @                            V                    #THIS W   #I X   #J Y   #FILENAME Z             
                                W                   #SPARSE              
                                 X                     
                                 Y                     
                               Z                    1 1         À    $                            [                  #PRINTNONZEROS \   #         @     @                            \                    #THIS ]   #FILENAME ^             
                                ]                   #SPARSE              
                               ^                    1 1         À    $                            _                  #PRINTALL `   #         @     @                            `                    #THIS a   #FILENAME b             
                                a                   #SPARSE              
                               b                    1 1         À    $                            c                  #DELETEROWANDCOL d   #         @     @                            d                    #THIS e   #ROW f   #COL g             
                                e                   #SPARSE              
                                 f                     
                                 g           1         À    $                            h                  #FREE i   #         @     @                            i                    #THIS j             
                                j                   #SPARSE    1         À    D                            k                  #HANDLEDUPLICATES l   #         @     @                            l                    #THIS m             
                                m                   #SPARSE                      @               @               '                    #PRECONDITIONERMETHOD n   #INIT o   #CHANGEPRECONDITIONER s   #USE w                $                             n                            #PRECONDITIONERDT    1         À    $                           o                  #INIT p   #         @     @                            p                    #THIS q   #METHOD r             
D                                q                    #USEPRECONDITIONERDT              
                                r                     #PRECONDITIONERDT    1         À    $                            s                  #CHANGEPRECONDITIONER t   #         @     @                             t                    #THIS u   #NEWMETHOD v             
D                                u                    #USEPRECONDITIONERDT              
                                v                     #PRECONDITIONERDT    1         À    $                            w                  #USE x   #         @     @                             x                    #THIS y   #VECTOR z   #MATRIX {   #SOLUTION |   #ARG }             
D @                              y                    #USEPRECONDITIONERDT              
D @                              z                   
               &                                                     
D @                              {                   #SPARSE              
D @                              |                   
               &                                                     
D @                              }                                  &                                                  S      fn#fn (   ó   :   b   uapp(USEPRECONDITIONERM    -  @   J  UTILITIESM    m  @   J  SPARSEKIT     ­  @   J  PRECONDITIONERM &   í  Q       gen@SETPRECONDITIONER    >  }      CONSTRUCTOR +   »  ^   a   CONSTRUCTOR%PRECONDITIONER 1     g      PRECONDITIONERDT+PRECONDITIONERM C     f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   æ        PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   g  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @   Å     a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @   Q  T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   ¥     a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   1     a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM !   ½  µ      SPARSE+SPARSEKIT %   r     %   SPARSE%A+SPARSEKIT=A '   	     %   SPARSE%AI+SPARSEKIT=AI '   	     %   SPARSE%AJ+SPARSEKIT=AJ 7   .
     %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   Â
  H   %   SPARSE%N+SPARSEKIT=N )   
  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   R  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1     ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   ÷  i      TRIPLET+SPARSEKIT $   `     a   TRIPLET%A+SPARSEKIT &   ô     a   TRIPLET%ROW+SPARSEKIT &        a   TRIPLET%COL+SPARSEKIT 5     H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   d  R   a   SPARSE%INIT+SPARSEKIT    ¶  e      INIT+SPARSEKIT $     T   a   INIT%THIS+SPARSEKIT #   o  @   a   INIT%NNZ+SPARSEKIT $   ¯  @   a   INIT%ROWS+SPARSEKIT (   ï  T   a   SPARSE%UPDATE+SPARSEKIT !   C  e      UPDATE+SPARSEKIT &   ¨  T   a   UPDATE%THIS+SPARSEKIT %   ü  @   a   UPDATE%NNZ+SPARSEKIT &   <  @   a   UPDATE%ROWS+SPARSEKIT (   |  T   a   SPARSE%APPEND+SPARSEKIT !   Ð  m      APPEND+SPARSEKIT &   =  T   a   APPEND%THIS+SPARSEKIT %     @   a   APPEND%VAL+SPARSEKIT %   Ñ  @   a   APPEND%ROW+SPARSEKIT %     @   a   APPEND%COL+SPARSEKIT )   Q  U   a   SPARSE%MAKECRS+SPARSEKIT "   ¦  `      MAKECRS+SPARSEKIT '     T   a   MAKECRS%THIS+SPARSEKIT +   Z  @   a   MAKECRS%SORTROWS+SPARSEKIT /     [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   õ  m      APPENDPOSTCRS+SPARSEKIT -   b  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   ¶  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   ö  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   6  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   v  T   a   SPARSE%CHANGE+SPARSEKIT !   Ê  m      CHANGE+SPARSEKIT &   7  T   a   CHANGE%THIS+SPARSEKIT %     @   a   CHANGE%VAL+SPARSEKIT %   Ë  @   a   CHANGE%ROW+SPARSEKIT %     @   a   CHANGE%COL+SPARSEKIT .   K  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   ¥  [      SETDIRICHLET+SPARSEKIT ,      T   a   SETDIRICHLET%THIS+SPARSEKIT +   T  @   a   SETDIRICHLET%ROW+SPARSEKIT %     Q   a   SPARSE%GET+SPARSEKIT    å  h      GET+SPARSEKIT #   M  T   a   GET%THIS+SPARSEKIT     ¡  @   a   GET%I+SPARSEKIT     á  @   a   GET%J+SPARSEKIT (   !  T   a   SPARSE%GETNNZ+SPARSEKIT !   u  Z      GETNNZ+SPARSEKIT &   Ï  T   a   GETNNZ%THIS+SPARSEKIT &   #  R   a   SPARSE%GETN+SPARSEKIT    u  Z      GETN+SPARSEKIT $   Ï  T   a   GETN%THIS+SPARSEKIT &   #  R   a   SPARSE%GETA+SPARSEKIT    u  ö      GETA+SPARSEKIT $   k  T   a   GETA%THIS+SPARSEKIT '   ¿  S   a   SPARSE%GETAI+SPARSEKIT       h     GETAI+SPARSEKIT %   z   T   a   GETAI%THIS+SPARSEKIT '   Î   S   a   SPARSE%GETAJ+SPARSEKIT     !!  ö      GETAJ+SPARSEKIT %   "  T   a   GETAJ%THIS+SPARSEKIT ,   k"  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   Ã"  n      PRINTVALUE+SPARSEKIT *   1#  T   a   PRINTVALUE%THIS+SPARSEKIT '   #  @   a   PRINTVALUE%I+SPARSEKIT '   Å#  @   a   PRINTVALUE%J+SPARSEKIT .   $  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   Q$  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   ¬$  `      PRINTNONZEROS+SPARSEKIT -   %  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   `%  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   ¬%  V   a   SPARSE%PRINTALL+SPARSEKIT #   &  `      PRINTALL+SPARSEKIT (   b&  T   a   PRINTALL%THIS+SPARSEKIT ,   ¶&  L   a   PRINTALL%FILENAME+SPARSEKIT 1   '  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   _'  d      DELETEROWANDCOL+SPARSEKIT /   Ã'  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   (  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   W(  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   (  R   a   SPARSE%FREE+SPARSEKIT    é(  R      FREE+SPARSEKIT $   ;)  T   a   FREE%THIS+SPARSEKIT C   )  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   í)  R      HANDLEDUPLICATES+SPARSEKIT 0   ?*  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT $   *         USEPRECONDITIONERDT 9   *+  f   a   USEPRECONDITIONERDT%PRECONDITIONERMETHOD )   +  R   a   USEPRECONDITIONERDT%INIT    â+  ^      INIT    @,  a   a   INIT%THIS    ¡,  ^   a   INIT%METHOD 9   ÿ,  b   a   USEPRECONDITIONERDT%CHANGEPRECONDITIONER %   a-  a      CHANGEPRECONDITIONER *   Â-  a   a   CHANGEPRECONDITIONER%THIS /   #.  ^   a   CHANGEPRECONDITIONER%NEWMETHOD (   .  Q   a   USEPRECONDITIONERDT%USE    Ò.        USE    S/  a   a   USE%THIS    ´/     a   USE%VECTOR    @0  T   a   USE%MATRIX    0     a   USE%SOLUTION     1     a   USE%ARG 