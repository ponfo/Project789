  �1  �   k820309              19.1        >��^                                                                                                          
       src/solvers/Linear/Iterative/Preconditioner/UsePreconditioner.f90 USEPRECONDITIONERM              USEPRECONDITIONERDT gen@SETPRECONDITIONER                                                     
                            @                              
                            @                              
                                                              u #CONSTRUCTOR    &         @   @                                �                      #PRECONDITIONER    #USEPRECONDITIONERDT              
D @                                                   #PRECONDITIONERDT                  @  @                               '                      #USEPRECONDITIONER    1         �   � $                     �                       #PRECONDITIONER_PROCEDURE 	   #         @     @                          	     	               #THIS 
   #VECTOR    #MATRIX    #SOLUTION    #ARG              
                               
                     #PRECONDITIONERDT              
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
                                m                   #SPARSE                      @               @               '�                    #PRECONDITIONERMETHOD n   #INIT o   #CHANGEPRECONDITIONER s   #USE w               � $                             n                            #PRECONDITIONERDT    1         �   � $                     �      o                  #INIT p   #         @     @                            p                    #THIS q   #METHOD r             
D                                q     �               #USEPRECONDITIONERDT              
                                r                     #PRECONDITIONERDT    1         �   � $                      �      s                  #CHANGEPRECONDITIONER t   #         @     @                             t                    #THIS u   #NEWMETHOD v             
D                                u     �               #USEPRECONDITIONERDT              
                                v                     #PRECONDITIONERDT    1         �   � $                      �      w                  #USE x   #         @     @                             x                    #THIS y   #VECTOR z   #MATRIX {   #SOLUTION |   #ARG }             
D @                              y     �               #USEPRECONDITIONERDT              
D @                              z                   
               &                                                     
D @                              {                   #SPARSE              
D @                              |                   
               &                                                     
D @                              }                                  &                                              �   ]      fn#fn (   �   :   b   uapp(USEPRECONDITIONERM    7  @   J  UTILITIESM    w  @   J  SPARSEKIT     �  @   J  PRECONDITIONERM &   �  Q       gen@SETPRECONDITIONER    H  }      CONSTRUCTOR +   �  ^   a   CONSTRUCTOR%PRECONDITIONER 1   #  g      PRECONDITIONERDT+PRECONDITIONERM C   �  f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   �  �      PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   q  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @   �  �   a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @   [  T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   �  �   a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   ;  �   a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM !   �  �      SPARSE+SPARSEKIT %   |  �   %   SPARSE%A+SPARSEKIT=A '   	  �   %   SPARSE%AI+SPARSEKIT=AI '   �	  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   8
  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �
  H   %   SPARSE%N+SPARSEKIT=N )     H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   \  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   �  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "     i      TRIPLET+SPARSEKIT $   j  �   a   TRIPLET%A+SPARSEKIT &   �  �   a   TRIPLET%ROW+SPARSEKIT &   �  �   a   TRIPLET%COL+SPARSEKIT 5   &  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   n  R   a   SPARSE%INIT+SPARSEKIT    �  e      INIT+SPARSEKIT $   %  T   a   INIT%THIS+SPARSEKIT #   y  @   a   INIT%NNZ+SPARSEKIT $   �  @   a   INIT%ROWS+SPARSEKIT (   �  T   a   SPARSE%UPDATE+SPARSEKIT !   M  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %     @   a   UPDATE%NNZ+SPARSEKIT &   F  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &   G  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %     @   a   APPEND%COL+SPARSEKIT )   [  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '     T   a   MAKECRS%THIS+SPARSEKIT +   d  @   a   MAKECRS%SORTROWS+SPARSEKIT /   �  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   �  m      APPENDPOSTCRS+SPARSEKIT -   l  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,      @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   @  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   �  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &   A  T   a   CHANGE%THIS+SPARSEKIT %   �  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %     @   a   CHANGE%COL+SPARSEKIT .   U  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   �  [      SETDIRICHLET+SPARSEKIT ,   
  T   a   SETDIRICHLET%THIS+SPARSEKIT +   ^  @   a   SETDIRICHLET%ROW+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT    �  h      GET+SPARSEKIT #   W  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT     �  @   a   GET%J+SPARSEKIT (   +  T   a   SPARSE%GETNNZ+SPARSEKIT !     Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   -  R   a   SPARSE%GETN+SPARSEKIT      Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &   -  R   a   SPARSE%GETA+SPARSEKIT      �      GETA+SPARSEKIT $   u  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT       h     GETAI+SPARSEKIT %   �   T   a   GETAI%THIS+SPARSEKIT '   �   S   a   SPARSE%GETAJ+SPARSEKIT     +!  �      GETAJ+SPARSEKIT %   !"  T   a   GETAJ%THIS+SPARSEKIT ,   u"  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �"  n      PRINTVALUE+SPARSEKIT *   ;#  T   a   PRINTVALUE%THIS+SPARSEKIT '   �#  @   a   PRINTVALUE%I+SPARSEKIT '   �#  @   a   PRINTVALUE%J+SPARSEKIT .   $  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   [$  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �$  `      PRINTNONZEROS+SPARSEKIT -   %  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   j%  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �%  V   a   SPARSE%PRINTALL+SPARSEKIT #   &  `      PRINTALL+SPARSEKIT (   l&  T   a   PRINTALL%THIS+SPARSEKIT ,   �&  L   a   PRINTALL%FILENAME+SPARSEKIT 1   '  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   i'  d      DELETEROWANDCOL+SPARSEKIT /   �'  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   !(  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   a(  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   �(  R   a   SPARSE%FREE+SPARSEKIT    �(  R      FREE+SPARSEKIT $   E)  T   a   FREE%THIS+SPARSEKIT C   �)  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �)  R      HANDLEDUPLICATES+SPARSEKIT 0   I*  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT $   �*  �       USEPRECONDITIONERDT 9   4+  f   a   USEPRECONDITIONERDT%PRECONDITIONERMETHOD )   �+  R   a   USEPRECONDITIONERDT%INIT    �+  ^      INIT    J,  a   a   INIT%THIS    �,  ^   a   INIT%METHOD 9   	-  b   a   USEPRECONDITIONERDT%CHANGEPRECONDITIONER %   k-  a      CHANGEPRECONDITIONER *   �-  a   a   CHANGEPRECONDITIONER%THIS /   -.  ^   a   CHANGEPRECONDITIONER%NEWMETHOD (   �.  Q   a   USEPRECONDITIONERDT%USE    �.  �      USE    ]/  a   a   USE%THIS    �/  �   a   USE%VECTOR    J0  T   a   USE%MATRIX    �0  �   a   USE%SOLUTION    *1  �   a   USE%ARG 