  �(  k   k820309              19.1        �3�^                                                                                                          
       src/solvers/NonLinear/NonLinearSolvers.f90 NONLINEARSOLVERSM              NONLINEARSOLVERSDT                                                     
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
                                a                   #SPARSE                      @                          b     '                      #USESOLVER c   1         �   � $                      �     c                  #NONLINEARSOLVERS_PROCEDURE d   #         @     @                            d     	               #THIS e   #MATRIX f   #VECTOR g   #SOLUTION h   #ARG i             
                               e                     #NONLINEARSOLVERSDT b             
                               f                   #SPARSE              
                               g                   
               &                                                     
                               h                   
               &                                                     
                               i                                  &                                              �   E      fn#fn '   �   #   b   uapp(NONLINEARSOLVERSM      @   J  UTILITIESM    H  @   J  SPARSEKIT !   �  �      SPARSE+SPARSEKIT %   =  �   %   SPARSE%A+SPARSEKIT=A '   �  �   %   SPARSE%AI+SPARSEKIT=AI '   e  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   �  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �  H   %   SPARSE%N+SPARSEKIT=N )   �  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1     H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   e  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   �  i      TRIPLET+SPARSEKIT $   +  �   a   TRIPLET%A+SPARSEKIT &   �  �   a   TRIPLET%ROW+SPARSEKIT &   S  �   a   TRIPLET%COL+SPARSEKIT 5   �  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   /	  R   a   SPARSE%INIT+SPARSEKIT    �	  e      INIT+SPARSEKIT $   �	  T   a   INIT%THIS+SPARSEKIT #   :
  @   a   INIT%NNZ+SPARSEKIT $   z
  @   a   INIT%ROWS+SPARSEKIT (   �
  T   a   SPARSE%UPDATE+SPARSEKIT !     e      UPDATE+SPARSEKIT &   s  T   a   UPDATE%THIS+SPARSEKIT %   �  @   a   UPDATE%NNZ+SPARSEKIT &     @   a   UPDATE%ROWS+SPARSEKIT (   G  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &     T   a   APPEND%THIS+SPARSEKIT %   \  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   �  @   a   APPEND%COL+SPARSEKIT )     U   a   SPARSE%MAKECRS+SPARSEKIT "   q  `      MAKECRS+SPARSEKIT '   �  T   a   MAKECRS%THIS+SPARSEKIT +   %  @   a   MAKECRS%SORTROWS+SPARSEKIT /   e  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   �  m      APPENDPOSTCRS+SPARSEKIT -   -  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,     @   a   APPENDPOSTCRS%COL+SPARSEKIT (   A  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &     T   a   CHANGE%THIS+SPARSEKIT %   V  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %   �  @   a   CHANGE%COL+SPARSEKIT .     Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   p  [      SETDIRICHLET+SPARSEKIT ,   �  T   a   SETDIRICHLET%THIS+SPARSEKIT +     @   a   SETDIRICHLET%ROW+SPARSEKIT %   _  Q   a   SPARSE%GET+SPARSEKIT    �  h      GET+SPARSEKIT #     T   a   GET%THIS+SPARSEKIT     l  @   a   GET%I+SPARSEKIT     �  @   a   GET%J+SPARSEKIT (   �  T   a   SPARSE%GETNNZ+SPARSEKIT !   @  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   �  R   a   SPARSE%GETN+SPARSEKIT    @  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &   �  R   a   SPARSE%GETA+SPARSEKIT    @  �      GETA+SPARSEKIT $   6  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     �  h     GETAI+SPARSEKIT %   E  T   a   GETAI%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAJ+SPARSEKIT     �  �      GETAJ+SPARSEKIT %   �  T   a   GETAJ%THIS+SPARSEKIT ,   6  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *   �  T   a   PRINTVALUE%THIS+SPARSEKIT '   P  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   �  L   a   PRINTVALUE%FILENAME+SPARSEKIT /     [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   w  `      PRINTNONZEROS+SPARSEKIT -   �  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   +   L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   w   V   a   SPARSE%PRINTALL+SPARSEKIT #   �   `      PRINTALL+SPARSEKIT (   -!  T   a   PRINTALL%THIS+SPARSEKIT ,   �!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �!  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   *"  d      DELETEROWANDCOL+SPARSEKIT /   �"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   �"  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   "#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   b#  R   a   SPARSE%FREE+SPARSEKIT    �#  R      FREE+SPARSEKIT $   $  T   a   FREE%THIS+SPARSEKIT C   Z$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �$  R      HANDLEDUPLICATES+SPARSEKIT 0   
%  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT #   ^%  _       NONLINEARSOLVERSDT -   �%  h   a   NONLINEARSOLVERSDT%USESOLVER +   %&  �      NONLINEARSOLVERS_PROCEDURE 0   �&  `   a   NONLINEARSOLVERS_PROCEDURE%THIS 2   '  T   a   NONLINEARSOLVERS_PROCEDURE%MATRIX 2   Z'  �   a   NONLINEARSOLVERS_PROCEDURE%VECTOR 4   �'  �   a   NONLINEARSOLVERS_PROCEDURE%SOLUTION /   r(  �   a   NONLINEARSOLVERS_PROCEDURE%ARG 