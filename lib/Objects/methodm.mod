  "-  u   k820309              19.1        �1�^                                                                                                          
       src/solvers/NonLinear/Method.f90 METHODM              METHODDT                                                     
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
                                b                   #SPARSE                      @                          c     '                      #NONLINEARSOLVERSDT d   #USESOLVER m                � $                              d                            #NONLINEARSOLVERSDT e                 @  @                          e     '                      #USESOLVER f   1         �   � $                      �     f                  #NONLINEARSOLVERS_PROCEDURE g   #         @     @                           g     	               #THIS h   #MATRIX i   #VECTOR j   #SOLUTION k   #ARG l             
                               h                     #NONLINEARSOLVERSDT e             
                               i                   #SPARSE              
                               j                   
               &                                                     
                               k                   
               &                                                     
                               l                                  &                                           1         �   � $                      �     m                  #METHOD n   #         @     @                             n                    #THIS o   #MATRIX p   #VECTOR q   #SOLUTION r   #ARG s             
                                o                     #METHODDT c             
                                p                   #SPARSE              
                                q                   
               &                                                     
                                r                   
               &                                                     
                                s                                  &                                              �   1      fn#fn    �      b   uapp(METHODM    �   @   J  UTILITIESM    *  @   J  SPARSEKIT "   j  @   J  NONLINEARSOLVERSM !   �  �      SPARSE+SPARSEKIT %   _  �   %   SPARSE%A+SPARSEKIT=A '   �  �   %   SPARSE%AI+SPARSEKIT=AI '   �  �   %   SPARSE%AJ+SPARSEKIT=AJ 7     �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �  H   %   SPARSE%N+SPARSEKIT=N )   �  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   ?  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   �  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   �  i      TRIPLET+SPARSEKIT $   M  �   a   TRIPLET%A+SPARSEKIT &   �  �   a   TRIPLET%ROW+SPARSEKIT &   u  �   a   TRIPLET%COL+SPARSEKIT 5   		  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   Q	  R   a   SPARSE%INIT+SPARSEKIT    �	  e      INIT+SPARSEKIT $   
  T   a   INIT%THIS+SPARSEKIT #   \
  @   a   INIT%NNZ+SPARSEKIT $   �
  @   a   INIT%ROWS+SPARSEKIT (   �
  T   a   SPARSE%UPDATE+SPARSEKIT !   0  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %   �  @   a   UPDATE%NNZ+SPARSEKIT &   )  @   a   UPDATE%ROWS+SPARSEKIT (   i  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &   *  T   a   APPEND%THIS+SPARSEKIT %   ~  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   �  @   a   APPEND%COL+SPARSEKIT )   >  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   �  T   a   MAKECRS%THIS+SPARSEKIT +   G  @   a   MAKECRS%SORTROWS+SPARSEKIT /   �  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   �  m      APPENDPOSTCRS+SPARSEKIT -   O  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   #  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   c  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &   $  T   a   CHANGE%THIS+SPARSEKIT %   x  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %   �  @   a   CHANGE%COL+SPARSEKIT .   8  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   �  [      SETDIRICHLET+SPARSEKIT ,   �  T   a   SETDIRICHLET%THIS+SPARSEKIT +   A  @   a   SETDIRICHLET%ROW+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT    �  h      GET+SPARSEKIT #   :  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT     �  @   a   GET%J+SPARSEKIT (     T   a   SPARSE%GETNNZ+SPARSEKIT !   b  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &     R   a   SPARSE%GETN+SPARSEKIT    b  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &     R   a   SPARSE%GETA+SPARSEKIT    b  �      GETA+SPARSEKIT $   X  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     �  h     GETAI+SPARSEKIT %   g  T   a   GETAI%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAJ+SPARSEKIT       �      GETAJ+SPARSEKIT %     T   a   GETAJ%THIS+SPARSEKIT ,   X  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *     T   a   PRINTVALUE%THIS+SPARSEKIT '   r  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   �  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   >  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �  `      PRINTNONZEROS+SPARSEKIT -   �  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   M   L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �   V   a   SPARSE%PRINTALL+SPARSEKIT #   �   `      PRINTALL+SPARSEKIT (   O!  T   a   PRINTALL%THIS+SPARSEKIT ,   �!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �!  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   L"  d      DELETEROWANDCOL+SPARSEKIT /   �"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   #  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   D#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   �#  R   a   SPARSE%FREE+SPARSEKIT    �#  R      FREE+SPARSEKIT $   ($  T   a   FREE%THIS+SPARSEKIT C   |$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �$  R      HANDLEDUPLICATES+SPARSEKIT 0   ,%  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT    �%  w       METHODDT ,   �%  h   a   METHODDT%NONLINEARSOLVERSDT 5   _&  _      NONLINEARSOLVERSDT+NONLINEARSOLVERSM ?   �&  h   a   NONLINEARSOLVERSDT%USESOLVER+NONLINEARSOLVERSM =   &'  �      NONLINEARSOLVERS_PROCEDURE+NONLINEARSOLVERSM B   �'  `   a   NONLINEARSOLVERS_PROCEDURE%THIS+NONLINEARSOLVERSM D   (  T   a   NONLINEARSOLVERS_PROCEDURE%MATRIX+NONLINEARSOLVERSM D   [(  �   a   NONLINEARSOLVERS_PROCEDURE%VECTOR+NONLINEARSOLVERSM F   �(  �   a   NONLINEARSOLVERS_PROCEDURE%SOLUTION+NONLINEARSOLVERSM A   s)  �   a   NONLINEARSOLVERS_PROCEDURE%ARG+NONLINEARSOLVERSM #   �)  T   a   METHODDT%USESOLVER    S*  �      METHOD    �*  V   a   METHOD%THIS    *+  T   a   METHOD%MATRIX    ~+  �   a   METHOD%VECTOR     
,  �   a   METHOD%SOLUTION    �,  �   a   METHOD%ARG 