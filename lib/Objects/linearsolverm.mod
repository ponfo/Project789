  &R  Í   k820309              19.1        ?!­^                                                                                                          
       src/solvers/Linear/LinearSolver.f90 LINEARSOLVERM              LINEARSOLVERDT gen@SETLINEARSOLVER                                                     
                            @                              
                            @                              
                            @                              
                            @                              
                                                              u #ITERATIVECONSTRUCTOR    #DIRECTCONSTRUCTOR 
   &         @   @                                                     #SOLVER    #LINEARSOLVERDT 	             
  @                                                 #ITERATIVELINEARSOLVERDT    &         @   @                           
                          #SOLVER    #LINEARSOLVERDT 	             
  @                                                  #DIRECTLINEARSOLVERDT                  @  @                               '                      #SOLVESYSTEM    1         À    $                                            #DIRECTLINEARSOLVER_PROCEDURE    #         @     @                               	               #THIS    #VECTOR    #MATRIX    #SOLUTION    #ARG              
                                                    #DIRECTLINEARSOLVERDT              
                                                  
               &                                                     
                                                  #SPARSE              
                                                  
               &                                                     
                                                                 &                                                         @  @                             '                    #PRECONDITIONER    #SOLVESYSTEM /                 $                                                         #USEPRECONDITIONERDT                  @  @              D               '                    #PRECONDITIONERMETHOD    #INIT     #CHANGEPRECONDITIONER $   #USE (                $                                                         #PRECONDITIONERDT                  @  @                              '                      #USEPRECONDITIONER    1         À    $                                             #PRECONDITIONER_PROCEDURE    #         @     @                                	               #THIS    #VECTOR    #MATRIX    #SOLUTION    #ARG              
                                                    #PRECONDITIONERDT              
                                                  
               &                                                     
                                                  #SPARSE              
                                                  
               &                                                     
                                                                 &                                           1         À    $                                               #INIT !   #         @     @                            !                    #THIS "   #METHOD #             
                                "                    #USEPRECONDITIONERDT              
                                #                     #PRECONDITIONERDT    1         À    $                            $                  #CHANGEPRECONDITIONER %   #         @     @                            %                    #THIS &   #NEWMETHOD '             
                                &                    #USEPRECONDITIONERDT              
                                '                     #PRECONDITIONERDT    1         À    $                            (                  #USE )   #         @     @                            )                    #THIS *   #VECTOR +   #MATRIX ,   #SOLUTION -   #ARG .             
                                *                    #USEPRECONDITIONERDT              
                                +                   
               &                                                     
                                ,                   #SPARSE              
                                -                   
               &                                                     
                                .                                  &                                           1         À    $                          /                  #ITERATIVELINEARSOLVER_PROCEDURE 0   #         @     @                          0     	               #THIS 1   #VECTOR 2   #MATRIX 3   #SOLUTION 4   #ARG 5             
                               1                    #ITERATIVELINEARSOLVERDT              
                               2                   
               &                                                     
                               3                   #SPARSE              
                               4                   
               &                                                     
                               5                                  &                                                             @               À                '                   #A 6   #AI 7   #AJ 8   #ROWCOUNTER 9   #N :   #NNZ ;   #COUNTER <   #TRIPLET =   #ISCRSDONE B   #INIT C   #UPDATE H   #APPEND M   #MAKECRS S   #APPENDPOSTCRS W   #CHANGE ]   #SETDIRICHLET c   #GET g   #GETNNZ l   #GETN o   #GETA r   #GETAI u   #GETAJ x   #PRINTVALUE {   #PRINTNONZEROS    #PRINTALL    #DELETEROWANDCOL    #FREE    #HANDLEDUPLICATES                D                             6                              
            &                                                      D                             7            H                             &                                                      D                             8                                         &                                                       D                             9            Ø                             &                                                         D                             :                                D                             ;     $                          D                             <     (                          D                              =     Ø       0             #TRIPLET >                 @  @              D           >     'Ø                    #A ?   #ROW @   #COL A                                            ?                              
            &                                                                                    @            H                             &                                                                                    A                                         &                                                         D                              B           	      1         À    $                            C             
     #INIT D   #         @     @                            D                    #THIS E   #NNZ F   #ROWS G             
                                E                   #SPARSE              
                                 F                     
                                 G           1         À    $                            H                  #UPDATE I   #         @     @                            I                    #THIS J   #NNZ K   #ROWS L             
                                J                   #SPARSE              
                                 K                     
                                 L           1         À    $                            M                  #APPEND N   #         @     @                            N                    #THIS O   #VAL P   #ROW Q   #COL R             
                                O                   #SPARSE              
                                 P     
                
                                 Q                     
                                 R           1         À    $                            S                  #MAKECRS T   #         @     @                            T                    #THIS U   #SORTROWS V             
                                U                   #SPARSE              
                                 V           1         À    $                            W                  #APPENDPOSTCRS X   #         @     @                            X                    #THIS Y   #VAL Z   #ROW [   #COL \             
                                Y                   #SPARSE              
                                 Z     
                
                                 [                     
                                 \           1         À    $                            ]                  #CHANGE ^   #         @     @                            ^                    #THIS _   #VAL `   #ROW a   #COL b             
                                _                   #SPARSE              
                                 `     
                
                                 a                     
                                 b           1         À    $                            c                  #SETDIRICHLET d   #         @     @                            d                    #THIS e   #ROW f             
                                e                   #SPARSE              
                                 f           1         À    $                           g                  #GET h   %         @   @                           h                    
       #THIS i   #I j   #J k             
                                i                   #SPARSE              
                                 j                     
                                 k           1         À    $                           l              	    #GETNNZ m   %         @   @                           m                           #THIS n             
                                n                   #SPARSE    1         À    $                           o              
    #GETN p   %         @   @                           p                           #THIS q             
                                q                   #SPARSE    1         À    $                           r                  #GETA s   (        `   @                           s                                    
    #THIS t   p          5 8 O#SPARSE     p        U     ;       5 8 O#SPARSE     p        U     ;                               
                                t                   #SPARSE    1         À    $                           u                  #GETAI v   (        `   @                           v                                        #THIS w   p           5 8 O#SPARSE     p        U     :   n                                         1     5 8 O#SPARSE     p        U     :   n                                          1                               
                                w                   #SPARSE    1         À    $                           x                  #GETAJ y   (        `   @                           y                                        #THIS z   p          5 8 O#SPARSE     p        U     ;       5 8 O#SPARSE     p        U     ;                               
                                z                   #SPARSE    1         À    $                            {                  #PRINTVALUE |   #         @     @                            |                    #THIS }   #I ~   #J    #FILENAME              
                                }                   #SPARSE              
                                 ~                     
                                                      
                                                   1 1         À    $                                              #PRINTNONZEROS    #         @     @                                                #THIS    #FILENAME              
                                                   #SPARSE              
                                                   1 1         À    $                                              #PRINTALL    #         @     @                                                #THIS    #FILENAME              
                                                   #SPARSE              
                                                   1 1         À    $                                              #DELETEROWANDCOL    #         @     @                                                #THIS    #ROW    #COL              
                                                   #SPARSE              
                                                      
                                            1         À    $                                              #FREE    #         @     @                                                #THIS              
                                                   #SPARSE    1         À    D                                              #HANDLEDUPLICATES    #         @     @                                                #THIS              
                                                   #SPARSE                  @  @                              '                      #USEREORDER    1         À    $                                             #REORDERSYSTEM_PROCEDURE    #         @     @                                	               #THIS    #VECTOR    #MATRIX    #SOLUTION    #ARG              
                                                    #REORDERSYSTEMDT              
                                                  
               &                                                     
                                                  #SPARSE              
                                                  
               &                                                     
                                                                 &                                                             @               À         	     '             
      #DIRECTSOLVER    #ITERATIVESOLVER    #REORDER    #INIT °   #ITERATIVEINIT ±   #DIRECTINIT ²   #CHANGESOLVER ¹   #CHANGEITERATIVE º   #CHANGEDIRECT »   #SOLVE Â                $                                                         #DIRECTLINEARSOLVERDT                 $                                                       #ITERATIVELINEARSOLVERDT                  $                                                        #USEREORDERSYSTEMDT                  @  @              @               '                    #REORDERMETHOD     #INIT ¡   #CHANGEMETHOD ¥   #USE ©                $                                                          #REORDERSYSTEMDT    1         À    $                            ¡                  #INIT ¢   #         @     @                            ¢                    #THIS £   #METHOD ¤             
                                £                    #USEREORDERSYSTEMDT              
                                ¤                     #REORDERSYSTEMDT    1         À    $                            ¥                  #CHANGEMETHOD ¦   #         @     @                            ¦                    #THIS §   #NEWMETHOD ¨             
                                §                    #USEREORDERSYSTEMDT              
                                ¨                     #REORDERSYSTEMDT    1         À    $                            ©                  #USE ª   #         @     @                            ª                    #THIS «   #VECTOR ¬   #MATRIX ­   #SOLUTION ®   #ARG ¯             
                                «                    #USEREORDERSYSTEMDT              
                                ¬                   
               &                                                     
                                ­                   #SPARSE              
                                ®                   
               &                                                     
                                ¯                                  &                                           4             $                         @    °                    3             $                         @             u #LINEARSOLVERDT 	   #ITERATIVEINIT ±   #DIRECTINIT ²   1         À    $                           ±                  #ITERATIVEINIT ³   #         @     @                            ³                    #THIS ´   #SOLVER µ             
D                                ´                   #LINEARSOLVERDT 	             
                                 µ                   #ITERATIVELINEARSOLVERDT    1         À    $                           ²                  #DIRECTINIT ¶   #         @     @                            ¶                    #THIS ·   #SOLVER ¸             
D                                ·                   #LINEARSOLVERDT 	             
                                 ¸                    #DIRECTLINEARSOLVERDT    4             $                         @    ¹                    3             $                         @             u #LINEARSOLVERDT 	   #CHANGEITERATIVE º   #CHANGEDIRECT »   1         À    $                            º                  #CHANGEITERATIVE ¼   #         @     @                             ¼                    #THIS ½   #NEWSOLVER ¾             
D                                ½                   #LINEARSOLVERDT 	             
                                ¾                    #ITERATIVELINEARSOLVERDT    1         À    $                            »             	     #CHANGEDIRECT ¿   #         @     @                             ¿                    #THIS À   #NEWSOLVER Á             
D                                À                   #LINEARSOLVERDT 	             
                                Á                     #DIRECTLINEARSOLVERDT    1         À    $                            Â             
     #USESOLVER Ã   #         @     @                             Ã                    #THIS Ä   #VECTOR Å   #MATRIX Æ   #SOLUTION Ç   #ARG È             
D @                              Ä                   #LINEARSOLVERDT 	             
D @                              Å                   
               &                                                     
D @                              Æ                   #SPARSE              
D @                              Ç                   
               &                                                     
D @                              È                                  &                                                  :      fn#fn #   Ú   3   b   uapp(LINEARSOLVERM      @   J  UTILITIESM    M  @   J  SPARSEKIT $     @   J  DIRECTLINEARSOLVERM '   Í  @   J  ITERATIVELINEARSOLVERM "     @   J  USEREORDERSYSTEMM $   M  q       gen@SETLINEARSOLVER %   ¾  p      ITERATIVECONSTRUCTOR ,   .  e   a   ITERATIVECONSTRUCTOR%SOLVER "     p      DIRECTCONSTRUCTOR )     b   a   DIRECTCONSTRUCTOR%SOLVER 9   e  a      DIRECTLINEARSOLVERDT+DIRECTLINEARSOLVERM E   Æ  j   a   DIRECTLINEARSOLVERDT%SOLVESYSTEM+DIRECTLINEARSOLVERM A   0        DIRECTLINEARSOLVER_PROCEDURE+DIRECTLINEARSOLVERM F   ±  b   a   DIRECTLINEARSOLVER_PROCEDURE%THIS+DIRECTLINEARSOLVERM H        a   DIRECTLINEARSOLVER_PROCEDURE%VECTOR+DIRECTLINEARSOLVERM H     T   a   DIRECTLINEARSOLVER_PROCEDURE%MATRIX+DIRECTLINEARSOLVERM J   ó     a   DIRECTLINEARSOLVER_PROCEDURE%SOLUTION+DIRECTLINEARSOLVERM E        a   DIRECTLINEARSOLVER_PROCEDURE%ARG+DIRECTLINEARSOLVERM ?     u      ITERATIVELINEARSOLVERDT+ITERATIVELINEARSOLVERM N     i   a   ITERATIVELINEARSOLVERDT%PRECONDITIONER+ITERATIVELINEARSOLVERM 7   é        USEPRECONDITIONERDT+USEPRECONDITIONERM L   	  f   a   USEPRECONDITIONERDT%PRECONDITIONERMETHOD+USEPRECONDITIONERM 1   æ	  g      PRECONDITIONERDT+PRECONDITIONERM C   M
  f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   ³
        PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   4  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @        a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @     T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   r     a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   þ     a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM <     R   a   USEPRECONDITIONERDT%INIT+USEPRECONDITIONERM (   Ü  ^      INIT+USEPRECONDITIONERM -   :  a   a   INIT%THIS+USEPRECONDITIONERM /     ^   a   INIT%METHOD+USEPRECONDITIONERM L   ù  b   a   USEPRECONDITIONERDT%CHANGEPRECONDITIONER+USEPRECONDITIONERM 8   [  a      CHANGEPRECONDITIONER+USEPRECONDITIONERM =   ¼  a   a   CHANGEPRECONDITIONER%THIS+USEPRECONDITIONERM B     ^   a   CHANGEPRECONDITIONER%NEWMETHOD+USEPRECONDITIONERM ;   {  Q   a   USEPRECONDITIONERDT%USE+USEPRECONDITIONERM '   Ì        USE+USEPRECONDITIONERM ,   M  a   a   USE%THIS+USEPRECONDITIONERM .   ®     a   USE%VECTOR+USEPRECONDITIONERM .   :  T   a   USE%MATRIX+USEPRECONDITIONERM 0        a   USE%SOLUTION+USEPRECONDITIONERM +        a   USE%ARG+USEPRECONDITIONERM K   ¦  m   a   ITERATIVELINEARSOLVERDT%SOLVESYSTEM+ITERATIVELINEARSOLVERM G           ITERATIVELINEARSOLVER_PROCEDURE+ITERATIVELINEARSOLVERM L     e   a   ITERATIVELINEARSOLVER_PROCEDURE%THIS+ITERATIVELINEARSOLVERM N   ù     a   ITERATIVELINEARSOLVER_PROCEDURE%VECTOR+ITERATIVELINEARSOLVERM N     T   a   ITERATIVELINEARSOLVER_PROCEDURE%MATRIX+ITERATIVELINEARSOLVERM P   Ù     a   ITERATIVELINEARSOLVER_PROCEDURE%SOLUTION+ITERATIVELINEARSOLVERM K   e     a   ITERATIVELINEARSOLVER_PROCEDURE%ARG+ITERATIVELINEARSOLVERM !   ñ  µ      SPARSE+SPARSEKIT %   ¦     %   SPARSE%A+SPARSEKIT=A '   :     %   SPARSE%AI+SPARSEKIT=AI '   Î     %   SPARSE%AJ+SPARSEKIT=AJ 7   b     %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   ö  H   %   SPARSE%N+SPARSEKIT=N )   >  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1     H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   Î  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   +  i      TRIPLET+SPARSEKIT $        a   TRIPLET%A+SPARSEKIT &   (     a   TRIPLET%ROW+SPARSEKIT &   ¼     a   TRIPLET%COL+SPARSEKIT 5   P  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &     R   a   SPARSE%INIT+SPARSEKIT    ê  e      INIT+SPARSEKIT $   O  T   a   INIT%THIS+SPARSEKIT #   £  @   a   INIT%NNZ+SPARSEKIT $   ã  @   a   INIT%ROWS+SPARSEKIT (   #   T   a   SPARSE%UPDATE+SPARSEKIT !   w   e      UPDATE+SPARSEKIT &   Ü   T   a   UPDATE%THIS+SPARSEKIT %   0!  @   a   UPDATE%NNZ+SPARSEKIT &   p!  @   a   UPDATE%ROWS+SPARSEKIT (   °!  T   a   SPARSE%APPEND+SPARSEKIT !   "  m      APPEND+SPARSEKIT &   q"  T   a   APPEND%THIS+SPARSEKIT %   Å"  @   a   APPEND%VAL+SPARSEKIT %   #  @   a   APPEND%ROW+SPARSEKIT %   E#  @   a   APPEND%COL+SPARSEKIT )   #  U   a   SPARSE%MAKECRS+SPARSEKIT "   Ú#  `      MAKECRS+SPARSEKIT '   :$  T   a   MAKECRS%THIS+SPARSEKIT +   $  @   a   MAKECRS%SORTROWS+SPARSEKIT /   Î$  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   )%  m      APPENDPOSTCRS+SPARSEKIT -   %  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   ê%  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   *&  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   j&  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   ª&  T   a   SPARSE%CHANGE+SPARSEKIT !   þ&  m      CHANGE+SPARSEKIT &   k'  T   a   CHANGE%THIS+SPARSEKIT %   ¿'  @   a   CHANGE%VAL+SPARSEKIT %   ÿ'  @   a   CHANGE%ROW+SPARSEKIT %   ?(  @   a   CHANGE%COL+SPARSEKIT .   (  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   Ù(  [      SETDIRICHLET+SPARSEKIT ,   4)  T   a   SETDIRICHLET%THIS+SPARSEKIT +   )  @   a   SETDIRICHLET%ROW+SPARSEKIT %   È)  Q   a   SPARSE%GET+SPARSEKIT    *  h      GET+SPARSEKIT #   *  T   a   GET%THIS+SPARSEKIT     Õ*  @   a   GET%I+SPARSEKIT     +  @   a   GET%J+SPARSEKIT (   U+  T   a   SPARSE%GETNNZ+SPARSEKIT !   ©+  Z      GETNNZ+SPARSEKIT &   ,  T   a   GETNNZ%THIS+SPARSEKIT &   W,  R   a   SPARSE%GETN+SPARSEKIT    ©,  Z      GETN+SPARSEKIT $   -  T   a   GETN%THIS+SPARSEKIT &   W-  R   a   SPARSE%GETA+SPARSEKIT    ©-  ö      GETA+SPARSEKIT $   .  T   a   GETA%THIS+SPARSEKIT '   ó.  S   a   SPARSE%GETAI+SPARSEKIT     F/  h     GETAI+SPARSEKIT %   ®0  T   a   GETAI%THIS+SPARSEKIT '   1  S   a   SPARSE%GETAJ+SPARSEKIT     U1  ö      GETAJ+SPARSEKIT %   K2  T   a   GETAJ%THIS+SPARSEKIT ,   2  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   ÷2  n      PRINTVALUE+SPARSEKIT *   e3  T   a   PRINTVALUE%THIS+SPARSEKIT '   ¹3  @   a   PRINTVALUE%I+SPARSEKIT '   ù3  @   a   PRINTVALUE%J+SPARSEKIT .   94  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   4  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   à4  `      PRINTNONZEROS+SPARSEKIT -   @5  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   5  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   à5  V   a   SPARSE%PRINTALL+SPARSEKIT #   66  `      PRINTALL+SPARSEKIT (   6  T   a   PRINTALL%THIS+SPARSEKIT ,   ê6  L   a   PRINTALL%FILENAME+SPARSEKIT 1   67  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   7  d      DELETEROWANDCOL+SPARSEKIT /   ÷7  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   K8  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   8  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   Ë8  R   a   SPARSE%FREE+SPARSEKIT    9  R      FREE+SPARSEKIT $   o9  T   a   FREE%THIS+SPARSEKIT C   Ã9  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   !:  R      HANDLEDUPLICATES+SPARSEKIT 0   s:  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT /   Ç:  `      REORDERSYSTEMDT+REORDERSYSTEMM :   ';  e   a   REORDERSYSTEMDT%USEREORDER+REORDERSYSTEMM 7   ;        REORDERSYSTEM_PROCEDURE+REORDERSYSTEMM <   <  ]   a   REORDERSYSTEM_PROCEDURE%THIS+REORDERSYSTEMM >   j<     a   REORDERSYSTEM_PROCEDURE%VECTOR+REORDERSYSTEMM >   ö<  T   a   REORDERSYSTEM_PROCEDURE%MATRIX+REORDERSYSTEMM @   J=     a   REORDERSYSTEM_PROCEDURE%SOLUTION+REORDERSYSTEMM ;   Ö=     a   REORDERSYSTEM_PROCEDURE%ARG+REORDERSYSTEMM    b>  õ       LINEARSOLVERDT ,   W?  j   a   LINEARSOLVERDT%DIRECTSOLVER /   Á?  m   a   LINEARSOLVERDT%ITERATIVESOLVER '   .@  h   a   LINEARSOLVERDT%REORDER 5   @        USEREORDERSYSTEMDT+USEREORDERSYSTEMM C   A  e   a   USEREORDERSYSTEMDT%REORDERMETHOD+USEREORDERSYSTEMM :   A  R   a   USEREORDERSYSTEMDT%INIT+USEREORDERSYSTEMM '   ÕA  ^      INIT+USEREORDERSYSTEMM ,   3B  `   a   INIT%THIS+USEREORDERSYSTEMM .   B  ]   a   INIT%METHOD+USEREORDERSYSTEMM B   ðB  Z   a   USEREORDERSYSTEMDT%CHANGEMETHOD+USEREORDERSYSTEMM /   JC  a      CHANGEMETHOD+USEREORDERSYSTEMM 4   «C  `   a   CHANGEMETHOD%THIS+USEREORDERSYSTEMM 9   D  ]   a   CHANGEMETHOD%NEWMETHOD+USEREORDERSYSTEMM 9   hD  Q   a   USEREORDERSYSTEMDT%USE+USEREORDERSYSTEMM &   ¹D        USE+USEREORDERSYSTEMM +   :E  `   a   USE%THIS+USEREORDERSYSTEMM -   E     a   USE%VECTOR+USEREORDERSYSTEMM -   &F  T   a   USE%MATRIX+USEREORDERSYSTEMM /   zF     a   USE%SOLUTION+USEREORDERSYSTEMM *   G     a   USE%ARG+USEREORDERSYSTEMM $   G  H   a   LINEARSOLVERDT%INIT    ÚG  w   `   gen@INIT -   QH  [   a   LINEARSOLVERDT%ITERATIVEINIT    ¬H  ^      ITERATIVEINIT #   
I  \   a   ITERATIVEINIT%THIS %   fI  e   a   ITERATIVEINIT%SOLVER *   ËI  X   a   LINEARSOLVERDT%DIRECTINIT    #J  ^      DIRECTINIT     J  \   a   DIRECTINIT%THIS "   ÝJ  b   a   DIRECTINIT%SOLVER ,   ?K  H   a   LINEARSOLVERDT%CHANGESOLVER !   K  {   `   gen@CHANGESOLVER /   L  ]   a   LINEARSOLVERDT%CHANGEITERATIVE     _L  a      CHANGEITERATIVE %   ÀL  \   a   CHANGEITERATIVE%THIS *   M  e   a   CHANGEITERATIVE%NEWSOLVER ,   M  Z   a   LINEARSOLVERDT%CHANGEDIRECT    ÛM  a      CHANGEDIRECT "   <N  \   a   CHANGEDIRECT%THIS '   N  b   a   CHANGEDIRECT%NEWSOLVER %   úN  W   a   LINEARSOLVERDT%SOLVE    QO        USESOLVER    ÒO  \   a   USESOLVER%THIS !   .P     a   USESOLVER%VECTOR !   ºP  T   a   USESOLVER%MATRIX #   Q     a   USESOLVER%SOLUTION    Q     a   USESOLVER%ARG 