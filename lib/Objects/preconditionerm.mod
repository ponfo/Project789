  )  k   k820309              19.1        V!�^                                                                                                          
       src/solvers/Linear/Preconditioner/Preconditioner.f90 PRECONDITIONERM              PRECONDITIONERDT                                                     
                            @                              
                      �  @               �                '                   #A    #AI    #AJ    #ROWCOUNTER    #N    #NNZ 	   #COUNTER 
   #TRIPLET    #ISCRSDONE    #INIT    #UPDATE    #APPEND    #MAKECRS !   #APPENDPOSTCRS %   #CHANGE +   #SETDIRICHLET 1   #GET 5   #GETNNZ :   #GETN =   #GETA @   #GETAI C   #GETAJ F   #PRINTVALUE I   #PRINTNONZEROS O   #PRINTALL S   #DELETEROWANDCOL W   #FREE \   #HANDLEDUPLICATES _              � D                                                           
            &                                                     � D                                         H                             &                                                     � D                                         �                             &                                                      � D                                         �                             &                                                        � D                                                            � D                             	     $                         � D                             
     (                         � D                                   �       0             #TRIPLET                  @  @              D                '�                    #A    #ROW    #COL               �                                                            
            &                                                      �                                          H                             &                                                      �                                          �                             &                                                        � D                                         	      1         �   � $                      �                   
     #INIT    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                        #UPDATE    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                        #APPEND    #         @     @                                                #THIS    #VAL    #ROW    #COL               
                                                   #SPARSE              
                                      
                
                                                      
                                             1         �   � $                      �      !                  #MAKECRS "   #         @     @                            "                    #THIS #   #SORTROWS $             
                                #                   #SPARSE              
                                 $           1         �   � $                      �      %                  #APPENDPOSTCRS &   #         @     @                            &                    #THIS '   #VAL (   #ROW )   #COL *             
                                '                   #SPARSE              
                                 (     
                
                                 )                     
                                 *           1         �   � $                      �      +                  #CHANGE ,   #         @     @                            ,                    #THIS -   #VAL .   #ROW /   #COL 0             
                                -                   #SPARSE              
                                 .     
                
                                 /                     
                                 0           1         �   � $                      �      1                  #SETDIRICHLET 2   #         @     @                            2                    #THIS 3   #ROW 4             
                                3                   #SPARSE              
                                 4           1         �   � $                     �      5                  #GET 6   %         @   @                           6                    
       #THIS 7   #I 8   #J 9             
                                7                   #SPARSE              
                                 8                     
                                 9           1         �   � $                     �      :              	    #GETNNZ ;   %         @   @                           ;                           #THIS <             
                                <                   #SPARSE    1         �   � $                     �      =              
    #GETN >   %         @   @                           >                           #THIS ?             
                                ?                   #SPARSE    1         �   � $                     �      @                  #GETA A   (        `   @                           A                                    
    #THIS B   p          5 8 O#SPARSE     p        U     	       5 8 O#SPARSE     p        U     	                               
                                B                   #SPARSE    1         �   � $                     �      C                  #GETAI D   (        `   @                           D                                        #THIS E   p           5 8 O#SPARSE     p        U        n                                         1     5 8 O#SPARSE     p        U        n                                          1                               
                                E                   #SPARSE    1         �   � $                     �      F                  #GETAJ G   (        `   @                           G                                        #THIS H   p          5 8 O#SPARSE     p        U     	       5 8 O#SPARSE     p        U     	                               
                                H                   #SPARSE    1         �   � $                      �      I                  #PRINTVALUE J   #         @     @                            J                    #THIS K   #I L   #J M   #FILENAME N             
                                K                   #SPARSE              
                                 L                     
                                 M                     
                               N                    1 1         �   � $                      �      O                  #PRINTNONZEROS P   #         @     @                            P                    #THIS Q   #FILENAME R             
                                Q                   #SPARSE              
                               R                    1 1         �   � $                      �      S                  #PRINTALL T   #         @     @                            T                    #THIS U   #FILENAME V             
                                U                   #SPARSE              
                               V                    1 1         �   � $                      �      W                  #DELETEROWANDCOL X   #         @     @                            X                    #THIS Y   #ROW Z   #COL [             
                                Y                   #SPARSE              
                                 Z                     
                                 [           1         �   � $                      �      \                  #FREE ]   #         @     @                            ]                    #THIS ^             
                                ^                   #SPARSE    1         �   � D                      �      _                  #HANDLEDUPLICATES `   #         @     @                            `                    #THIS a             
                                a                   #SPARSE                      @                          b     '                      #USEPRECONDITIONER c   1         �   � $                      �     c                  #PRECONDITIONER_PROCEDURE d   #         @     @                            d     	               #THIS e   #VECTOR f   #MATRIX g   #SOLUTION h   #ARG i             
                               e                     #PRECONDITIONERDT b             
                               f                   
               &                                                     
                               g                   #SPARSE              
                               h                   
               &                                                     
                               i                                  &                                              �   M      fn#fn %   �   !   b   uapp(PRECONDITIONERM      @   J  UTILITIESM    N  @   J  SPARSEKIT !   �  �      SPARSE+SPARSEKIT %   C  �   %   SPARSE%A+SPARSEKIT=A '   �  �   %   SPARSE%AI+SPARSEKIT=AI '   k  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   �  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �  H   %   SPARSE%N+SPARSEKIT=N )   �  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   #  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   k  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   �  i      TRIPLET+SPARSEKIT $   1  �   a   TRIPLET%A+SPARSEKIT &   �  �   a   TRIPLET%ROW+SPARSEKIT &   Y  �   a   TRIPLET%COL+SPARSEKIT 5   �  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   5	  R   a   SPARSE%INIT+SPARSEKIT    �	  e      INIT+SPARSEKIT $   �	  T   a   INIT%THIS+SPARSEKIT #   @
  @   a   INIT%NNZ+SPARSEKIT $   �
  @   a   INIT%ROWS+SPARSEKIT (   �
  T   a   SPARSE%UPDATE+SPARSEKIT !     e      UPDATE+SPARSEKIT &   y  T   a   UPDATE%THIS+SPARSEKIT %   �  @   a   UPDATE%NNZ+SPARSEKIT &     @   a   UPDATE%ROWS+SPARSEKIT (   M  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &     T   a   APPEND%THIS+SPARSEKIT %   b  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   �  @   a   APPEND%COL+SPARSEKIT )   "  U   a   SPARSE%MAKECRS+SPARSEKIT "   w  `      MAKECRS+SPARSEKIT '   �  T   a   MAKECRS%THIS+SPARSEKIT +   +  @   a   MAKECRS%SORTROWS+SPARSEKIT /   k  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   �  m      APPENDPOSTCRS+SPARSEKIT -   3  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,     @   a   APPENDPOSTCRS%COL+SPARSEKIT (   G  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &     T   a   CHANGE%THIS+SPARSEKIT %   \  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %   �  @   a   CHANGE%COL+SPARSEKIT .     Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   v  [      SETDIRICHLET+SPARSEKIT ,   �  T   a   SETDIRICHLET%THIS+SPARSEKIT +   %  @   a   SETDIRICHLET%ROW+SPARSEKIT %   e  Q   a   SPARSE%GET+SPARSEKIT    �  h      GET+SPARSEKIT #     T   a   GET%THIS+SPARSEKIT     r  @   a   GET%I+SPARSEKIT     �  @   a   GET%J+SPARSEKIT (   �  T   a   SPARSE%GETNNZ+SPARSEKIT !   F  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   �  R   a   SPARSE%GETN+SPARSEKIT    F  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &   �  R   a   SPARSE%GETA+SPARSEKIT    F  �      GETA+SPARSEKIT $   <  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     �  h     GETAI+SPARSEKIT %   K  T   a   GETAI%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAJ+SPARSEKIT     �  �      GETAJ+SPARSEKIT %   �  T   a   GETAJ%THIS+SPARSEKIT ,   <  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *     T   a   PRINTVALUE%THIS+SPARSEKIT '   V  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   �  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   "  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   }  `      PRINTNONZEROS+SPARSEKIT -   �  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   1   L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   }   V   a   SPARSE%PRINTALL+SPARSEKIT #   �   `      PRINTALL+SPARSEKIT (   3!  T   a   PRINTALL%THIS+SPARSEKIT ,   �!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �!  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   0"  d      DELETEROWANDCOL+SPARSEKIT /   �"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   �"  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   (#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   h#  R   a   SPARSE%FREE+SPARSEKIT    �#  R      FREE+SPARSEKIT $   $  T   a   FREE%THIS+SPARSEKIT C   `$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �$  R      HANDLEDUPLICATES+SPARSEKIT 0   %  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT !   d%  g       PRECONDITIONERDT 3   �%  f   a   PRECONDITIONERDT%USEPRECONDITIONER )   1&  �      PRECONDITIONER_PROCEDURE .   �&  ^   a   PRECONDITIONER_PROCEDURE%THIS 0   '  �   a   PRECONDITIONER_PROCEDURE%VECTOR 0   �'  T   a   PRECONDITIONER_PROCEDURE%MATRIX 2   �'  �   a   PRECONDITIONER_PROCEDURE%SOLUTION -   |(  �   a   PRECONDITIONER_PROCEDURE%ARG 