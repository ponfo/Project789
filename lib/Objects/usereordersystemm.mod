  q1  �   k820309              19.1        ��^                                                                                                          
       src/solvers/Linear/Reorder/UseReorderSystem.f90 USEREORDERSYSTEMM              USEREORDERSYSTEMDT gen@SETREORDERSYSTEM                                                     
                            @                              
                            @                              
                                                              u #CONSTRUCTOR    &         @   @                                �                      #METHOD    #USEREORDERSYSTEMDT              
D @                                                   #REORDERSYSTEMDT                  @  @                               '                      #USEREORDER    1         �   � $                     �                       #REORDERSYSTEM_PROCEDURE 	   #         @     @                          	     	               #THIS 
   #VECTOR    #MATRIX    #SOLUTION    #ARG              
                               
                     #REORDERSYSTEMDT              
                                                  
               &                                                     
                                                  #SPARSE              
                                                  
               &                                                     
                                                                 &                                                          �  @               �                '                   #A    #AI    #AJ    #ROWCOUNTER    #N    #NNZ    #COUNTER    #TRIPLET    #ISCRSDONE    #INIT    #UPDATE "   #APPEND '   #MAKECRS -   #APPENDPOSTCRS 1   #CHANGE 7   #SETDIRICHLET =   #GET A   #GETNNZ F   #GETN I   #GETA L   #GETAI O   #GETAJ R   #PRINTVALUE U   #PRINTNONZEROS [   #PRINTALL _   #DELETEROWANDCOL c   #FREE h   #HANDLEDUPLICATES k              � D                                                           
            &                                                     � D                                         H                             &                                                     � D                                         �                             &                                                      � D                                         �                             &                                                        � D                                                            � D                                  $                         � D                                  (                         � D                                   �       0             #TRIPLET                  @  @              D                '�                    #A    #ROW    #COL               �                                                            
            &                                                      �                                          H                             &                                                      �                                          �                             &                                                        � D                                         	      1         �   � $                      �                   
     #INIT    #         @     @                                                #THIS    #NNZ     #ROWS !             
                                                   #SPARSE              
                                                       
                                 !           1         �   � $                      �      "                  #UPDATE #   #         @     @                            #                    #THIS $   #NNZ %   #ROWS &             
                                $                   #SPARSE              
                                 %                     
                                 &           1         �   � $                      �      '                  #APPEND (   #         @     @                            (                    #THIS )   #VAL *   #ROW +   #COL ,             
                                )                   #SPARSE              
                                 *     
                
                                 +                     
                                 ,           1         �   � $                      �      -                  #MAKECRS .   #         @     @                            .                    #THIS /   #SORTROWS 0             
                                /                   #SPARSE              
                                 0           1         �   � $                      �      1                  #APPENDPOSTCRS 2   #         @     @                            2                    #THIS 3   #VAL 4   #ROW 5   #COL 6             
                                3                   #SPARSE              
                                 4     
                
                                 5                     
                                 6           1         �   � $                      �      7                  #CHANGE 8   #         @     @                            8                    #THIS 9   #VAL :   #ROW ;   #COL <             
                                9                   #SPARSE              
                                 :     
                
                                 ;                     
                                 <           1         �   � $                      �      =                  #SETDIRICHLET >   #         @     @                            >                    #THIS ?   #ROW @             
                                ?                   #SPARSE              
                                 @           1         �   � $                     �      A                  #GET B   %         @   @                           B                    
       #THIS C   #I D   #J E             
                                C                   #SPARSE              
                                 D                     
                                 E           1         �   � $                     �      F              	    #GETNNZ G   %         @   @                           G                           #THIS H             
                                H                   #SPARSE    1         �   � $                     �      I              
    #GETN J   %         @   @                           J                           #THIS K             
                                K                   #SPARSE    1         �   � $                     �      L                  #GETA M   (        `   @                           M                                    
    #THIS N   p          5 8 O#SPARSE     p        U            5 8 O#SPARSE     p        U                                    
                                N                   #SPARSE    1         �   � $                     �      O                  #GETAI P   (        `   @                           P                                        #THIS Q   p           5 8 O#SPARSE     p        U        n                                         1     5 8 O#SPARSE     p        U        n                                          1                               
                                Q                   #SPARSE    1         �   � $                     �      R                  #GETAJ S   (        `   @                           S                                        #THIS T   p          5 8 O#SPARSE     p        U            5 8 O#SPARSE     p        U                                    
                                T                   #SPARSE    1         �   � $                      �      U                  #PRINTVALUE V   #         @     @                            V                    #THIS W   #I X   #J Y   #FILENAME Z             
                                W                   #SPARSE              
                                 X                     
                                 Y                     
                               Z                    1 1         �   � $                      �      [                  #PRINTNONZEROS \   #         @     @                            \                    #THIS ]   #FILENAME ^             
                                ]                   #SPARSE              
                               ^                    1 1         �   � $                      �      _                  #PRINTALL `   #         @     @                            `                    #THIS a   #FILENAME b             
                                a                   #SPARSE              
                               b                    1 1         �   � $                      �      c                  #DELETEROWANDCOL d   #         @     @                            d                    #THIS e   #ROW f   #COL g             
                                e                   #SPARSE              
                                 f                     
                                 g           1         �   � $                      �      h                  #FREE i   #         @     @                            i                    #THIS j             
                                j                   #SPARSE    1         �   � D                      �      k                  #HANDLEDUPLICATES l   #         @     @                            l                    #THIS m             
                                m                   #SPARSE                      @               @               '�                    #REORDERMETHOD n   #INIT o   #CHANGEMETHOD s   #USE w               � $                             n                            #REORDERSYSTEMDT    1         �   � $                     �      o                  #INIT p   #         @     @                            p                    #THIS q   #METHOD r             
D                                q     �               #USEREORDERSYSTEMDT              
                                r                     #REORDERSYSTEMDT    1         �   � $                      �      s                  #CHANGEMETHOD t   #         @     @                             t                    #THIS u   #NEWMETHOD v             
D                                u     �               #USEREORDERSYSTEMDT              
                                v                     #REORDERSYSTEMDT    1         �   � $                      �      w                  #USE x   #         @     @                             x                    #THIS y   #VECTOR z   #MATRIX {   #SOLUTION |   #ARG }             
D @                              y     �               #USEREORDERSYSTEMDT              
D @                              z                   
               &                                                     
D @                              {                   #SPARSE              
D @                              |                   
               &                                                     
D @                              }                                  &                                              �   J      fn#fn '   �   8   b   uapp(USEREORDERSYSTEMM    "  @   J  UTILITIESM    b  @   J  SPARSEKIT    �  @   J  REORDERSYSTEMM %   �  Q       gen@SETREORDERSYSTEM    3  t      CONSTRUCTOR #   �  ]   a   CONSTRUCTOR%METHOD /     `      REORDERSYSTEMDT+REORDERSYSTEMM :   d  e   a   REORDERSYSTEMDT%USEREORDER+REORDERSYSTEMM 7   �  �      REORDERSYSTEM_PROCEDURE+REORDERSYSTEMM <   J  ]   a   REORDERSYSTEM_PROCEDURE%THIS+REORDERSYSTEMM >   �  �   a   REORDERSYSTEM_PROCEDURE%VECTOR+REORDERSYSTEMM >   3  T   a   REORDERSYSTEM_PROCEDURE%MATRIX+REORDERSYSTEMM @   �  �   a   REORDERSYSTEM_PROCEDURE%SOLUTION+REORDERSYSTEMM ;     �   a   REORDERSYSTEM_PROCEDURE%ARG+REORDERSYSTEMM !   �  �      SPARSE+SPARSEKIT %   T  �   %   SPARSE%A+SPARSEKIT=A '   �  �   %   SPARSE%AI+SPARSEKIT=AI '   |	  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   
  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �
  H   %   SPARSE%N+SPARSEKIT=N )   �
  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   4  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   |  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   �  i      TRIPLET+SPARSEKIT $   B  �   a   TRIPLET%A+SPARSEKIT &   �  �   a   TRIPLET%ROW+SPARSEKIT &   j  �   a   TRIPLET%COL+SPARSEKIT 5   �  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   F  R   a   SPARSE%INIT+SPARSEKIT    �  e      INIT+SPARSEKIT $   �  T   a   INIT%THIS+SPARSEKIT #   Q  @   a   INIT%NNZ+SPARSEKIT $   �  @   a   INIT%ROWS+SPARSEKIT (   �  T   a   SPARSE%UPDATE+SPARSEKIT !   %  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %   �  @   a   UPDATE%NNZ+SPARSEKIT &     @   a   UPDATE%ROWS+SPARSEKIT (   ^  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &     T   a   APPEND%THIS+SPARSEKIT %   s  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   �  @   a   APPEND%COL+SPARSEKIT )   3  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   �  T   a   MAKECRS%THIS+SPARSEKIT +   <  @   a   MAKECRS%SORTROWS+SPARSEKIT /   |  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   �  m      APPENDPOSTCRS+SPARSEKIT -   D  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,     @   a   APPENDPOSTCRS%COL+SPARSEKIT (   X  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &     T   a   CHANGE%THIS+SPARSEKIT %   m  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %   �  @   a   CHANGE%COL+SPARSEKIT .   -  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   �  [      SETDIRICHLET+SPARSEKIT ,   �  T   a   SETDIRICHLET%THIS+SPARSEKIT +   6  @   a   SETDIRICHLET%ROW+SPARSEKIT %   v  Q   a   SPARSE%GET+SPARSEKIT    �  h      GET+SPARSEKIT #   /  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT     �  @   a   GET%J+SPARSEKIT (     T   a   SPARSE%GETNNZ+SPARSEKIT !   W  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &     R   a   SPARSE%GETN+SPARSEKIT    W  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &     R   a   SPARSE%GETA+SPARSEKIT    W  �      GETA+SPARSEKIT $   M  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     �  h     GETAI+SPARSEKIT %   \   T   a   GETAI%THIS+SPARSEKIT '   �   S   a   SPARSE%GETAJ+SPARSEKIT     !  �      GETAJ+SPARSEKIT %   �!  T   a   GETAJ%THIS+SPARSEKIT ,   M"  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �"  n      PRINTVALUE+SPARSEKIT *   #  T   a   PRINTVALUE%THIS+SPARSEKIT '   g#  @   a   PRINTVALUE%I+SPARSEKIT '   �#  @   a   PRINTVALUE%J+SPARSEKIT .   �#  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   3$  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �$  `      PRINTNONZEROS+SPARSEKIT -   �$  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   B%  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �%  V   a   SPARSE%PRINTALL+SPARSEKIT #   �%  `      PRINTALL+SPARSEKIT (   D&  T   a   PRINTALL%THIS+SPARSEKIT ,   �&  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �&  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   A'  d      DELETEROWANDCOL+SPARSEKIT /   �'  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   �'  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   9(  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   y(  R   a   SPARSE%FREE+SPARSEKIT    �(  R      FREE+SPARSEKIT $   )  T   a   FREE%THIS+SPARSEKIT C   q)  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �)  R      HANDLEDUPLICATES+SPARSEKIT 0   !*  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT #   u*  �       USEREORDERSYSTEMDT 1   �*  e   a   USEREORDERSYSTEMDT%REORDERMETHOD (   b+  R   a   USEREORDERSYSTEMDT%INIT    �+  ^      INIT    ,  `   a   INIT%THIS    r,  ]   a   INIT%METHOD 0   �,  Z   a   USEREORDERSYSTEMDT%CHANGEMETHOD    )-  a      CHANGEMETHOD "   �-  `   a   CHANGEMETHOD%THIS '   �-  ]   a   CHANGEMETHOD%NEWMETHOD '   G.  Q   a   USEREORDERSYSTEMDT%USE    �.  �      USE    /  `   a   USE%THIS    y/  �   a   USE%VECTOR    0  T   a   USE%MATRIX    Y0  �   a   USE%SOLUTION    �0  �   a   USE%ARG 