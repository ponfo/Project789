  �8  �   k820309              19.1        �ٳ^                                                                                                          
       src/solvers/Linear/Iterative/Solvers/IterativeLinearSolverMethod.f90 ITERATIVELINEARSOLVERMETHODM              ITERATIVELINEARSOLVERMETHODDT                                                     
                            @                              
                            @                              
                      �  @               �                '                   #A    #AI    #AJ    #ROWCOUNTER    #N 	   #NNZ 
   #COUNTER    #TRIPLET    #ISCRSDONE    #INIT    #UPDATE    #APPEND    #MAKECRS "   #APPENDPOSTCRS &   #CHANGE ,   #SETDIRICHLET 2   #GET 6   #GETNNZ ;   #GETN >   #GETA A   #GETAI D   #GETAJ G   #PRINTVALUE J   #PRINTNONZEROS P   #PRINTALL T   #DELETEROWANDCOL X   #FREE ]   #HANDLEDUPLICATES `              � D                                                           
            &                                                     � D                                         H                             &                                                     � D                                         �                             &                                                      � D                                         �                             &                                                        � D                             	                               � D                             
     $                         � D                                  (                         � D                                   �       0             #TRIPLET                  @  @              D                '�                    #A    #ROW    #COL               �                                                            
            &                                                      �                                          H                             &                                                      �                                          �                             &                                                        � D                                         	      1         �   � $                      �                   
     #INIT    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                        #UPDATE    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                        #APPEND    #         @     @                                                #THIS    #VAL    #ROW     #COL !             
                                                   #SPARSE              
                                      
                
                                                       
                                 !           1         �   � $                      �      "                  #MAKECRS #   #         @     @                            #                    #THIS $   #SORTROWS %             
                                $                   #SPARSE              
                                 %           1         �   � $                      �      &                  #APPENDPOSTCRS '   #         @     @                            '                    #THIS (   #VAL )   #ROW *   #COL +             
                                (                   #SPARSE              
                                 )     
                
                                 *                     
                                 +           1         �   � $                      �      ,                  #CHANGE -   #         @     @                            -                    #THIS .   #VAL /   #ROW 0   #COL 1             
                                .                   #SPARSE              
                                 /     
                
                                 0                     
                                 1           1         �   � $                      �      2                  #SETDIRICHLET 3   #         @     @                            3                    #THIS 4   #ROW 5             
                                4                   #SPARSE              
                                 5           1         �   � $                     �      6                  #GET 7   %         @   @                           7                    
       #THIS 8   #I 9   #J :             
                                8                   #SPARSE              
                                 9                     
                                 :           1         �   � $                     �      ;              	    #GETNNZ <   %         @   @                           <                           #THIS =             
                                =                   #SPARSE    1         �   � $                     �      >              
    #GETN ?   %         @   @                           ?                           #THIS @             
                                @                   #SPARSE    1         �   � $                     �      A                  #GETA B   (        `   @                           B                                    
    #THIS C   p          5 8 O#SPARSE     p        U     
       5 8 O#SPARSE     p        U     
                               
                                C                   #SPARSE    1         �   � $                     �      D                  #GETAI E   (        `   @                           E                                        #THIS F   p           5 8 O#SPARSE     p        U     	   n                                         1     5 8 O#SPARSE     p        U     	   n                                          1                               
                                F                   #SPARSE    1         �   � $                     �      G                  #GETAJ H   (        `   @                           H                                        #THIS I   p          5 8 O#SPARSE     p        U     
       5 8 O#SPARSE     p        U     
                               
                                I                   #SPARSE    1         �   � $                      �      J                  #PRINTVALUE K   #         @     @                            K                    #THIS L   #I M   #J N   #FILENAME O             
                                L                   #SPARSE              
                                 M                     
                                 N                     
                               O                    1 1         �   � $                      �      P                  #PRINTNONZEROS Q   #         @     @                            Q                    #THIS R   #FILENAME S             
                                R                   #SPARSE              
                               S                    1 1         �   � $                      �      T                  #PRINTALL U   #         @     @                            U                    #THIS V   #FILENAME W             
                                V                   #SPARSE              
                               W                    1 1         �   � $                      �      X                  #DELETEROWANDCOL Y   #         @     @                            Y                    #THIS Z   #ROW [   #COL \             
                                Z                   #SPARSE              
                                 [                     
                                 \           1         �   � $                      �      ]                  #FREE ^   #         @     @                            ^                    #THIS _             
                                _                   #SPARSE    1         �   � D                      �      `                  #HANDLEDUPLICATES a   #         @     @                            a                    #THIS b             
                                b                   #SPARSE                  @  @              D          c     '�                    #PRECONDITIONERMETHOD d   #INIT m   #CHANGEPRECONDITIONER q   #USE u               � $                             d                            #PRECONDITIONERDT e                 @  @                         e     '                      #USEPRECONDITIONER f   1         �   � $                      �     f                  #PRECONDITIONER_PROCEDURE g   #         @     @                           g     	               #THIS h   #VECTOR i   #MATRIX j   #SOLUTION k   #ARG l             
                               h                     #PRECONDITIONERDT e             
                               i                   
               &                                                     
                               j                   #SPARSE              
                               k                   
               &                                                     
                               l                                  &                                           1         �   � $                      �      m                  #INIT n   #         @     @                            n                    #THIS o   #METHOD p             
                                o     �               #USEPRECONDITIONERDT c             
                                p                     #PRECONDITIONERDT e   1         �   � $                      �      q                  #CHANGEPRECONDITIONER r   #         @     @                            r                    #THIS s   #NEWMETHOD t             
                                s     �               #USEPRECONDITIONERDT c             
                                t                     #PRECONDITIONERDT e   1         �   � $                      �      u                  #USE v   #         @     @                            v                    #THIS w   #VECTOR x   #MATRIX y   #SOLUTION z   #ARG {             
                                w     �               #USEPRECONDITIONERDT c             
                                x                   
               &                                                     
                                y                   #SPARSE              
                                z                   
               &                                                     
                                {                                  &                                                             @               �         |     '�                    #ITERATIVELINEARSOLVERDT }   #SOLVESYSTEM �                � $                              }     �                      #ITERATIVELINEARSOLVERDT ~                 @  @               �         ~     '�                    #PRECONDITIONER    #SOLVESYSTEM �                � $                                   �                      #USEPRECONDITIONERDT c   1         �   � $                      �     �                  #ITERATIVELINEARSOLVER_PROCEDURE �   #         @     @                           �     	               #THIS �   #VECTOR �   #MATRIX �   #SOLUTION �   #ARG �             
                               �     �               #ITERATIVELINEARSOLVERDT ~             
                               �                   
               &                                                     
                               �                   #SPARSE              
                               �                   
               &                                                     
                               �                                  &                                           1         �   � $                      �     �                  #METHOD �   #         @     @                             �                    #THIS �   #VECTOR �   #MATRIX �   #SOLUTION �   #ARG �             
                                �     �               #ITERATIVELINEARSOLVERMETHODDT |             
                                �                   
               &                                                     
                                �                   #SPARSE              
                                �                   
               &                                                     
                                �                                  &                                              �   j      fn#fn 2   
  .   b   uapp(ITERATIVELINEARSOLVERMETHODM    8  @   J  UTILITIESM    x  @   J  SPARSEKIT '   �  @   J  ITERATIVELINEARSOLVERM !   �  �      SPARSE+SPARSEKIT %   �  �   %   SPARSE%A+SPARSEKIT=A '   A  �   %   SPARSE%AI+SPARSEKIT=AI '   �  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   i  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �  H   %   SPARSE%N+SPARSEKIT=N )   E  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   �  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   �  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   2  i      TRIPLET+SPARSEKIT $   �  �   a   TRIPLET%A+SPARSEKIT &   /  �   a   TRIPLET%ROW+SPARSEKIT &   �  �   a   TRIPLET%COL+SPARSEKIT 5   W	  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   �	  R   a   SPARSE%INIT+SPARSEKIT    �	  e      INIT+SPARSEKIT $   V
  T   a   INIT%THIS+SPARSEKIT #   �
  @   a   INIT%NNZ+SPARSEKIT $   �
  @   a   INIT%ROWS+SPARSEKIT (   *  T   a   SPARSE%UPDATE+SPARSEKIT !   ~  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %   7  @   a   UPDATE%NNZ+SPARSEKIT &   w  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !     m      APPEND+SPARSEKIT &   x  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %     @   a   APPEND%ROW+SPARSEKIT %   L  @   a   APPEND%COL+SPARSEKIT )   �  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   A  T   a   MAKECRS%THIS+SPARSEKIT +   �  @   a   MAKECRS%SORTROWS+SPARSEKIT /   �  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   0  m      APPENDPOSTCRS+SPARSEKIT -   �  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   1  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   q  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   �  T   a   SPARSE%CHANGE+SPARSEKIT !     m      CHANGE+SPARSEKIT &   r  T   a   CHANGE%THIS+SPARSEKIT %   �  @   a   CHANGE%VAL+SPARSEKIT %     @   a   CHANGE%ROW+SPARSEKIT %   F  @   a   CHANGE%COL+SPARSEKIT .   �  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   �  [      SETDIRICHLET+SPARSEKIT ,   ;  T   a   SETDIRICHLET%THIS+SPARSEKIT +   �  @   a   SETDIRICHLET%ROW+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT       h      GET+SPARSEKIT #   �  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT       @   a   GET%J+SPARSEKIT (   \  T   a   SPARSE%GETNNZ+SPARSEKIT !   �  Z      GETNNZ+SPARSEKIT &   
  T   a   GETNNZ%THIS+SPARSEKIT &   ^  R   a   SPARSE%GETN+SPARSEKIT    �  Z      GETN+SPARSEKIT $   
  T   a   GETN%THIS+SPARSEKIT &   ^  R   a   SPARSE%GETA+SPARSEKIT    �  �      GETA+SPARSEKIT $   �  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     M  h     GETAI+SPARSEKIT %   �  T   a   GETAI%THIS+SPARSEKIT '   	  S   a   SPARSE%GETAJ+SPARSEKIT     \  �      GETAJ+SPARSEKIT %   R  T   a   GETAJ%THIS+SPARSEKIT ,   �  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *   l  T   a   PRINTVALUE%THIS+SPARSEKIT '   �  @   a   PRINTVALUE%I+SPARSEKIT '      @   a   PRINTVALUE%J+SPARSEKIT .   @  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   �  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �  `      PRINTNONZEROS+SPARSEKIT -   G   T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   �   L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �   V   a   SPARSE%PRINTALL+SPARSEKIT #   =!  `      PRINTALL+SPARSEKIT (   �!  T   a   PRINTALL%THIS+SPARSEKIT ,   �!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   ="  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   �"  d      DELETEROWANDCOL+SPARSEKIT /   �"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   R#  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   �#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   �#  R   a   SPARSE%FREE+SPARSEKIT    $$  R      FREE+SPARSEKIT $   v$  T   a   FREE%THIS+SPARSEKIT C   �$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   (%  R      HANDLEDUPLICATES+SPARSEKIT 0   z%  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT 7   �%  �      USEPRECONDITIONERDT+USEPRECONDITIONERM L   e&  f   a   USEPRECONDITIONERDT%PRECONDITIONERMETHOD+USEPRECONDITIONERM 1   �&  g      PRECONDITIONERDT+PRECONDITIONERM C   2'  f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   �'  �      PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   (  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @   w(  �   a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @   )  T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   W)  �   a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   �)  �   a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM <   o*  R   a   USEPRECONDITIONERDT%INIT+USEPRECONDITIONERM (   �*  ^      INIT+USEPRECONDITIONERM -   +  a   a   INIT%THIS+USEPRECONDITIONERM /   �+  ^   a   INIT%METHOD+USEPRECONDITIONERM L   �+  b   a   USEPRECONDITIONERDT%CHANGEPRECONDITIONER+USEPRECONDITIONERM 8   @,  a      CHANGEPRECONDITIONER+USEPRECONDITIONERM =   �,  a   a   CHANGEPRECONDITIONER%THIS+USEPRECONDITIONERM B   -  ^   a   CHANGEPRECONDITIONER%NEWMETHOD+USEPRECONDITIONERM ;   `-  Q   a   USEPRECONDITIONERDT%USE+USEPRECONDITIONERM '   �-  �      USE+USEPRECONDITIONERM ,   2.  a   a   USE%THIS+USEPRECONDITIONERM .   �.  �   a   USE%VECTOR+USEPRECONDITIONERM .   /  T   a   USE%MATRIX+USEPRECONDITIONERM 0   s/  �   a   USE%SOLUTION+USEPRECONDITIONERM +   �/  �   a   USE%ARG+USEPRECONDITIONERM .   �0  ~       ITERATIVELINEARSOLVERMETHODDT F   	1  m   a   ITERATIVELINEARSOLVERMETHODDT%ITERATIVELINEARSOLVERDT ?   v1  u      ITERATIVELINEARSOLVERDT+ITERATIVELINEARSOLVERM N   �1  i   a   ITERATIVELINEARSOLVERDT%PRECONDITIONER+ITERATIVELINEARSOLVERM K   T2  m   a   ITERATIVELINEARSOLVERDT%SOLVESYSTEM+ITERATIVELINEARSOLVERM G   �2  �      ITERATIVELINEARSOLVER_PROCEDURE+ITERATIVELINEARSOLVERM L   B3  e   a   ITERATIVELINEARSOLVER_PROCEDURE%THIS+ITERATIVELINEARSOLVERM N   �3  �   a   ITERATIVELINEARSOLVER_PROCEDURE%VECTOR+ITERATIVELINEARSOLVERM N   34  T   a   ITERATIVELINEARSOLVER_PROCEDURE%MATRIX+ITERATIVELINEARSOLVERM P   �4  �   a   ITERATIVELINEARSOLVER_PROCEDURE%SOLUTION+ITERATIVELINEARSOLVERM K   5  �   a   ITERATIVELINEARSOLVER_PROCEDURE%ARG+ITERATIVELINEARSOLVERM :   �5  T   a   ITERATIVELINEARSOLVERMETHODDT%SOLVESYSTEM    �5  �      METHOD    t6  k   a   METHOD%THIS    �6  �   a   METHOD%VECTOR    k7  T   a   METHOD%MATRIX     �7  �   a   METHOD%SOLUTION    K8  �   a   METHOD%ARG 