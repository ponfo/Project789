  �4  �   k820309              19.1        ?��^                                                                                                          
       src/solvers/Linear/Iterative/IterativeLinearSolver.f90 ITERATIVELINEARSOLVERM              ITERATIVELINEARSOLVERDT                                                     
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
                                b                   #SPARSE                  @  @                         c     '                      #USEPRECONDITIONER d   1         �   � $                      �     d                  #PRECONDITIONER_PROCEDURE e   #         @     @                           e     	               #THIS f   #VECTOR g   #MATRIX h   #SOLUTION i   #ARG j             
                               f                     #PRECONDITIONERDT c             
                               g                   
               &                                                     
                               h                   #SPARSE              
                               i                   
               &                                                     
                               j                                  &                                                             @               �         k     '�                    #PRECONDITIONER l   #SOLVESYSTEM ~                � $                              l     �                      #USEPRECONDITIONERDT m                 @  @              @          m     '�                    #PRECONDITIONERMETHOD n   #INIT o   #CHANGEPRECONDITIONER s   #USE w               � $                             n                            #PRECONDITIONERDT c   1         �   � $                      �      o                  #INIT p   #         @     @                            p                    #THIS q   #METHOD r             
                                q     �               #USEPRECONDITIONERDT m             
                                r                     #PRECONDITIONERDT c   1         �   � $                      �      s                  #CHANGEPRECONDITIONER t   #         @     @                            t                    #THIS u   #NEWMETHOD v             
                                u     �               #USEPRECONDITIONERDT m             
                                v                     #PRECONDITIONERDT c   1         �   � $                      �      w                  #USE x   #         @     @                            x                    #THIS y   #VECTOR z   #MATRIX {   #SOLUTION |   #ARG }             
                                y     �               #USEPRECONDITIONERDT m             
                                z                   
               &                                                     
                                {                   #SPARSE              
                                |                   
               &                                                     
                                }                                  &                                           1         �   � $                      �     ~                  #ITERATIVELINEARSOLVER_PROCEDURE    #         @     @                                 	               #THIS �   #VECTOR �   #MATRIX �   #SOLUTION �   #ARG �             
                               �     �               #ITERATIVELINEARSOLVERDT k             
                               �                   
               &                                                     
                               �                   #SPARSE              
                               �                   
               &                                                     
                               �                                  &                                              �   V      fn#fn ,   �   (   b   uapp(ITERATIVELINEARSOLVERM      @   J  UTILITIESM    ^  @   J  SPARSEKIT #   �  @   J  USEPRECONDITIONERM !   �  �      SPARSE+SPARSEKIT %   �  �   %   SPARSE%A+SPARSEKIT=A '   '  �   %   SPARSE%AI+SPARSEKIT=AI '   �  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   O  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �  H   %   SPARSE%N+SPARSEKIT=N )   +  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   s  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   �  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "     i      TRIPLET+SPARSEKIT $   �  �   a   TRIPLET%A+SPARSEKIT &     �   a   TRIPLET%ROW+SPARSEKIT &   �  �   a   TRIPLET%COL+SPARSEKIT 5   =	  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   �	  R   a   SPARSE%INIT+SPARSEKIT    �	  e      INIT+SPARSEKIT $   <
  T   a   INIT%THIS+SPARSEKIT #   �
  @   a   INIT%NNZ+SPARSEKIT $   �
  @   a   INIT%ROWS+SPARSEKIT (     T   a   SPARSE%UPDATE+SPARSEKIT !   d  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %     @   a   UPDATE%NNZ+SPARSEKIT &   ]  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &   ^  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   2  @   a   APPEND%COL+SPARSEKIT )   r  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   '  T   a   MAKECRS%THIS+SPARSEKIT +   {  @   a   MAKECRS%SORTROWS+SPARSEKIT /   �  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (     m      APPENDPOSTCRS+SPARSEKIT -   �  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,     @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   W  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   �  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &   X  T   a   CHANGE%THIS+SPARSEKIT %   �  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %   ,  @   a   CHANGE%COL+SPARSEKIT .   l  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   �  [      SETDIRICHLET+SPARSEKIT ,   !  T   a   SETDIRICHLET%THIS+SPARSEKIT +   u  @   a   SETDIRICHLET%ROW+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT      h      GET+SPARSEKIT #   n  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT       @   a   GET%J+SPARSEKIT (   B  T   a   SPARSE%GETNNZ+SPARSEKIT !   �  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   D  R   a   SPARSE%GETN+SPARSEKIT    �  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &   D  R   a   SPARSE%GETA+SPARSEKIT    �  �      GETA+SPARSEKIT $   �  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     3  h     GETAI+SPARSEKIT %   �  T   a   GETAI%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAJ+SPARSEKIT     B  �      GETAJ+SPARSEKIT %   8  T   a   GETAJ%THIS+SPARSEKIT ,   �  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *   R  T   a   PRINTVALUE%THIS+SPARSEKIT '   �  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   &  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   r  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �  `      PRINTNONZEROS+SPARSEKIT -   -   T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   �   L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �   V   a   SPARSE%PRINTALL+SPARSEKIT #   #!  `      PRINTALL+SPARSEKIT (   �!  T   a   PRINTALL%THIS+SPARSEKIT ,   �!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   #"  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   �"  d      DELETEROWANDCOL+SPARSEKIT /   �"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   8#  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   x#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   �#  R   a   SPARSE%FREE+SPARSEKIT    
$  R      FREE+SPARSEKIT $   \$  T   a   FREE%THIS+SPARSEKIT C   �$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   %  R      HANDLEDUPLICATES+SPARSEKIT 0   `%  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT 1   �%  g      PRECONDITIONERDT+PRECONDITIONERM C   &  f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   �&  �      PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   '  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @   `'  �   a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @   �'  T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   @(  �   a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   �(  �   a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM (   X)  u       ITERATIVELINEARSOLVERDT 7   �)  i   a   ITERATIVELINEARSOLVERDT%PRECONDITIONER 7   6*  �      USEPRECONDITIONERDT+USEPRECONDITIONERM L   �*  f   a   USEPRECONDITIONERDT%PRECONDITIONERMETHOD+USEPRECONDITIONERM <   3+  R   a   USEPRECONDITIONERDT%INIT+USEPRECONDITIONERM (   �+  ^      INIT+USEPRECONDITIONERM -   �+  a   a   INIT%THIS+USEPRECONDITIONERM /   D,  ^   a   INIT%METHOD+USEPRECONDITIONERM L   �,  b   a   USEPRECONDITIONERDT%CHANGEPRECONDITIONER+USEPRECONDITIONERM 8   -  a      CHANGEPRECONDITIONER+USEPRECONDITIONERM =   e-  a   a   CHANGEPRECONDITIONER%THIS+USEPRECONDITIONERM B   �-  ^   a   CHANGEPRECONDITIONER%NEWMETHOD+USEPRECONDITIONERM ;   $.  Q   a   USEPRECONDITIONERDT%USE+USEPRECONDITIONERM '   u.  �      USE+USEPRECONDITIONERM ,   �.  a   a   USE%THIS+USEPRECONDITIONERM .   W/  �   a   USE%VECTOR+USEPRECONDITIONERM .   �/  T   a   USE%MATRIX+USEPRECONDITIONERM 0   70  �   a   USE%SOLUTION+USEPRECONDITIONERM +   �0  �   a   USE%ARG+USEPRECONDITIONERM 4   O1  m   a   ITERATIVELINEARSOLVERDT%SOLVESYSTEM 0   �1  �      ITERATIVELINEARSOLVER_PROCEDURE 5   =2  e   a   ITERATIVELINEARSOLVER_PROCEDURE%THIS 7   �2  �   a   ITERATIVELINEARSOLVER_PROCEDURE%VECTOR 7   .3  T   a   ITERATIVELINEARSOLVER_PROCEDURE%MATRIX 9   �3  �   a   ITERATIVELINEARSOLVER_PROCEDURE%SOLUTION 4   4  �   a   ITERATIVELINEARSOLVER_PROCEDURE%ARG 