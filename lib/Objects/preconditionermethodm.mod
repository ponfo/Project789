  x-  u   k820309              19.1        ��^                                                                                                          
       src/solvers/Linear/Iterative/Preconditioner/PreconditionerMethod.f90 PRECONDITIONERMETHODM              PRECONDITIONERMETHODDT                                                     
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
                                b                   #SPARSE                      @                          c     '                      #PRECONDITIONERDT d   #USEPRECONDITIONER m                � $                              d                            #PRECONDITIONERDT e                 @  @                          e     '                      #USEPRECONDITIONER f   1         �   � $                      �     f                  #PRECONDITIONER_PROCEDURE g   #         @     @                           g     	               #THIS h   #VECTOR i   #MATRIX j   #SOLUTION k   #ARG l             
                               h                     #PRECONDITIONERDT e             
                               i                   
               &                                                     
                               j                   #SPARSE              
                               k                   
               &                                                     
                               l                                  &                                           1         �   � $                      �     m                  #METHOD n   #         @     @                             n                    #THIS o   #VECTOR p   #MATRIX q   #SOLUTION r   #ARG s             
                                o                     #PRECONDITIONERMETHODDT c             
                                p                   
               &                                                     
                                q                   #SPARSE              
                                r                   
               &                                                     
                                s                                  &                                              �   c      fn#fn +     '   b   uapp(PRECONDITIONERMETHODM    *  @   J  UTILITIESM    j  @   J  SPARSEKIT     �  @   J  PRECONDITIONERM !   �  �      SPARSE+SPARSEKIT %   �  �   %   SPARSE%A+SPARSEKIT=A '   3  �   %   SPARSE%AI+SPARSEKIT=AI '   �  �   %   SPARSE%AJ+SPARSEKIT=AJ 7   [  �   %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   �  H   %   SPARSE%N+SPARSEKIT=N )   7  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1     H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   �  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "   $  i      TRIPLET+SPARSEKIT $   �  �   a   TRIPLET%A+SPARSEKIT &   !  �   a   TRIPLET%ROW+SPARSEKIT &   �  �   a   TRIPLET%COL+SPARSEKIT 5   I	  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   �	  R   a   SPARSE%INIT+SPARSEKIT    �	  e      INIT+SPARSEKIT $   H
  T   a   INIT%THIS+SPARSEKIT #   �
  @   a   INIT%NNZ+SPARSEKIT $   �
  @   a   INIT%ROWS+SPARSEKIT (     T   a   SPARSE%UPDATE+SPARSEKIT !   p  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %   )  @   a   UPDATE%NNZ+SPARSEKIT &   i  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &   j  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   >  @   a   APPEND%COL+SPARSEKIT )   ~  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   3  T   a   MAKECRS%THIS+SPARSEKIT +   �  @   a   MAKECRS%SORTROWS+SPARSEKIT /   �  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (   "  m      APPENDPOSTCRS+SPARSEKIT -   �  T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   �  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,   #  @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   c  @   a   APPENDPOSTCRS%COL+SPARSEKIT (   �  T   a   SPARSE%CHANGE+SPARSEKIT !   �  m      CHANGE+SPARSEKIT &   d  T   a   CHANGE%THIS+SPARSEKIT %   �  @   a   CHANGE%VAL+SPARSEKIT %   �  @   a   CHANGE%ROW+SPARSEKIT %   8  @   a   CHANGE%COL+SPARSEKIT .   x  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   �  [      SETDIRICHLET+SPARSEKIT ,   -  T   a   SETDIRICHLET%THIS+SPARSEKIT +   �  @   a   SETDIRICHLET%ROW+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT      h      GET+SPARSEKIT #   z  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT       @   a   GET%J+SPARSEKIT (   N  T   a   SPARSE%GETNNZ+SPARSEKIT !   �  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   P  R   a   SPARSE%GETN+SPARSEKIT    �  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT &   P  R   a   SPARSE%GETA+SPARSEKIT    �  �      GETA+SPARSEKIT $   �  T   a   GETA%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAI+SPARSEKIT     ?  h     GETAI+SPARSEKIT %   �  T   a   GETAI%THIS+SPARSEKIT '   �  S   a   SPARSE%GETAJ+SPARSEKIT     N  �      GETAJ+SPARSEKIT %   D  T   a   GETAJ%THIS+SPARSEKIT ,   �  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *   ^  T   a   PRINTVALUE%THIS+SPARSEKIT '   �  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   2  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   ~  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �  `      PRINTNONZEROS+SPARSEKIT -   9   T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   �   L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �   V   a   SPARSE%PRINTALL+SPARSEKIT #   /!  `      PRINTALL+SPARSEKIT (   �!  T   a   PRINTALL%THIS+SPARSEKIT ,   �!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   /"  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   �"  d      DELETEROWANDCOL+SPARSEKIT /   �"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   D#  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   �#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   �#  R   a   SPARSE%FREE+SPARSEKIT    $  R      FREE+SPARSEKIT $   h$  T   a   FREE%THIS+SPARSEKIT C   �$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   %  R      HANDLEDUPLICATES+SPARSEKIT 0   l%  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT '   �%  }       PRECONDITIONERMETHODDT 8   =&  f   a   PRECONDITIONERMETHODDT%PRECONDITIONERDT 1   �&  g      PRECONDITIONERDT+PRECONDITIONERM C   
'  f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   p'  �      PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   �'  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @   O(  �   a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @   �(  T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   /)  �   a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   �)  �   a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM 9   G*  T   a   PRECONDITIONERMETHODDT%USEPRECONDITIONER    �*  �      METHOD    +  d   a   METHOD%THIS    �+  �   a   METHOD%VECTOR    ,  T   a   METHOD%MATRIX     `,  �   a   METHOD%SOLUTION    �,  �   a   METHOD%ARG 