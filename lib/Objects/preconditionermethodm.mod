  n-  u   k820309              19.1        ÚIĆ^                                                                                                          
       src/solvers/Linear/Preconditioner/PreconditionerMethod.f90 PRECONDITIONERMETHODM              PRECONDITIONERMETHODDT                                                     
                            @                              
                            @                              
                         @               À                '                   #A    #AI    #AJ    #ROWCOUNTER    #N 	   #NNZ 
   #COUNTER    #TRIPLET    #ISCRSDONE    #INIT    #UPDATE    #APPEND    #MAKECRS "   #APPENDPOSTCRS &   #CHANGE ,   #SETDIRICHLET 2   #GET 6   #GETNNZ ;   #GETN >   #GETA A   #GETAI D   #GETAJ G   #PRINTVALUE J   #PRINTNONZEROS P   #PRINTALL T   #DELETEROWANDCOL X   #FREE ]   #HANDLEDUPLICATES `               D                                                           
            &                                                      D                                         H                             &                                                      D                                                                      &                                                       D                                         Ű                             &                                                         D                             	                                D                             
     $                          D                                  (                          D                                   Ű       0             #TRIPLET                  @  @              D                'Ű                    #A    #ROW    #COL                                                                           
            &                                                                                                H                             &                                                                                                                             &                                                         D                                         	      1         À    $                                         
     #INIT    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         À    $                                              #UPDATE    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         À    $                                              #APPEND    #         @     @                                                #THIS    #VAL    #ROW     #COL !             
                                                   #SPARSE              
                                      
                
                                                       
                                 !           1         À    $                            "                  #MAKECRS #   #         @     @                            #                    #THIS $   #SORTROWS %             
                                $                   #SPARSE              
                                 %           1         À    $                            &                  #APPENDPOSTCRS '   #         @     @                            '                    #THIS (   #VAL )   #ROW *   #COL +             
                                (                   #SPARSE              
                                 )     
                
                                 *                     
                                 +           1         À    $                            ,                  #CHANGE -   #         @     @                            -                    #THIS .   #VAL /   #ROW 0   #COL 1             
                                .                   #SPARSE              
                                 /     
                
                                 0                     
                                 1           1         À    $                            2                  #SETDIRICHLET 3   #         @     @                            3                    #THIS 4   #ROW 5             
                                4                   #SPARSE              
                                 5           1         À    $                           6                  #GET 7   %         @   @                           7                    
       #THIS 8   #I 9   #J :             
                                8                   #SPARSE              
                                 9                     
                                 :           1         À    $                           ;              	    #GETNNZ <   %         @   @                           <                           #THIS =             
                                =                   #SPARSE    1         À    $                           >              
    #GETN ?   %         @   @                           ?                           #THIS @             
                                @                   #SPARSE    1         À    $                           A                  #GETA B   (        `   @                           B                                    
    #THIS C   p          5 8 O#SPARSE     p        U     
       5 8 O#SPARSE     p        U     
                               
                                C                   #SPARSE    1         À    $                           D                  #GETAI E   (        `   @                           E                                        #THIS F   p           5 8 O#SPARSE     p        U     	   n                                         1     5 8 O#SPARSE     p        U     	   n                                          1                               
                                F                   #SPARSE    1         À    $                           G                  #GETAJ H   (        `   @                           H                                        #THIS I   p          5 8 O#SPARSE     p        U     
       5 8 O#SPARSE     p        U     
                               
                                I                   #SPARSE    1         À    $                            J                  #PRINTVALUE K   #         @     @                            K                    #THIS L   #I M   #J N   #FILENAME O             
                                L                   #SPARSE              
                                 M                     
                                 N                     
                               O                    1 1         À    $                            P                  #PRINTNONZEROS Q   #         @     @                            Q                    #THIS R   #FILENAME S             
                                R                   #SPARSE              
                               S                    1 1         À    $                            T                  #PRINTALL U   #         @     @                            U                    #THIS V   #FILENAME W             
                                V                   #SPARSE              
                               W                    1 1         À    $                            X                  #DELETEROWANDCOL Y   #         @     @                            Y                    #THIS Z   #ROW [   #COL \             
                                Z                   #SPARSE              
                                 [                     
                                 \           1         À    $                            ]                  #FREE ^   #         @     @                            ^                    #THIS _             
                                _                   #SPARSE    1         À    D                            `                  #HANDLEDUPLICATES a   #         @     @                            a                    #THIS b             
                                b                   #SPARSE                      @                          c     '                      #PRECONDITIONERDT d   #USEPRECONDITIONER m                 $                              d                            #PRECONDITIONERDT e                 @  @                          e     '                      #USEPRECONDITIONER f   1         À    $                           f                  #PRECONDITIONER_PROCEDURE g   #         @     @                           g     	               #THIS h   #VECTOR i   #MATRIX j   #SOLUTION k   #ARG l             
                               h                     #PRECONDITIONERDT e             
                               i                   
               &                                                     
                               j                   #SPARSE              
                               k                   
               &                                                     
                               l                                  &                                           1         À    $                           m                  #METHOD n   #         @     @                             n                    #THIS o   #VECTOR p   #MATRIX q   #SOLUTION r   #ARG s             
                                o                     #PRECONDITIONERMETHODDT c             
                                p                   
               &                                                     
                                q                   #SPARSE              
                                r                   
               &                                                     
                                s                                  &                                                  Y      fn#fn +   ù   '   b   uapp(PRECONDITIONERMETHODM       @   J  UTILITIESM    `  @   J  SPARSEKIT        @   J  PRECONDITIONERM !   à  ”      SPARSE+SPARSEKIT %        %   SPARSE%A+SPARSEKIT=A '   )     %   SPARSE%AI+SPARSEKIT=AI '   œ     %   SPARSE%AJ+SPARSEKIT=AJ 7   Q     %   SPARSE%ROWCOUNTER+SPARSEKIT=ROWCOUNTER %   ć  H   %   SPARSE%N+SPARSEKIT=N )   -  H   %   SPARSE%NNZ+SPARSEKIT=NNZ 1   u  H   %   SPARSE%COUNTER+SPARSEKIT=COUNTER 1   œ  ]   %   SPARSE%TRIPLET+SPARSEKIT=TRIPLET "     i      TRIPLET+SPARSEKIT $        a   TRIPLET%A+SPARSEKIT &        a   TRIPLET%ROW+SPARSEKIT &   «     a   TRIPLET%COL+SPARSEKIT 5   ?	  H   %   SPARSE%ISCRSDONE+SPARSEKIT=ISCRSDONE &   	  R   a   SPARSE%INIT+SPARSEKIT    Ù	  e      INIT+SPARSEKIT $   >
  T   a   INIT%THIS+SPARSEKIT #   
  @   a   INIT%NNZ+SPARSEKIT $   Ò
  @   a   INIT%ROWS+SPARSEKIT (     T   a   SPARSE%UPDATE+SPARSEKIT !   f  e      UPDATE+SPARSEKIT &   Ë  T   a   UPDATE%THIS+SPARSEKIT %     @   a   UPDATE%NNZ+SPARSEKIT &   _  @   a   UPDATE%ROWS+SPARSEKIT (     T   a   SPARSE%APPEND+SPARSEKIT !   ó  m      APPEND+SPARSEKIT &   `  T   a   APPEND%THIS+SPARSEKIT %   Ž  @   a   APPEND%VAL+SPARSEKIT %   ô  @   a   APPEND%ROW+SPARSEKIT %   4  @   a   APPEND%COL+SPARSEKIT )   t  U   a   SPARSE%MAKECRS+SPARSEKIT "   É  `      MAKECRS+SPARSEKIT '   )  T   a   MAKECRS%THIS+SPARSEKIT +   }  @   a   MAKECRS%SORTROWS+SPARSEKIT /   œ  [   a   SPARSE%APPENDPOSTCRS+SPARSEKIT (     m      APPENDPOSTCRS+SPARSEKIT -     T   a   APPENDPOSTCRS%THIS+SPARSEKIT ,   Ù  @   a   APPENDPOSTCRS%VAL+SPARSEKIT ,     @   a   APPENDPOSTCRS%ROW+SPARSEKIT ,   Y  @   a   APPENDPOSTCRS%COL+SPARSEKIT (     T   a   SPARSE%CHANGE+SPARSEKIT !   í  m      CHANGE+SPARSEKIT &   Z  T   a   CHANGE%THIS+SPARSEKIT %   ź  @   a   CHANGE%VAL+SPARSEKIT %   î  @   a   CHANGE%ROW+SPARSEKIT %   .  @   a   CHANGE%COL+SPARSEKIT .   n  Z   a   SPARSE%SETDIRICHLET+SPARSEKIT '   È  [      SETDIRICHLET+SPARSEKIT ,   #  T   a   SETDIRICHLET%THIS+SPARSEKIT +   w  @   a   SETDIRICHLET%ROW+SPARSEKIT %   ·  Q   a   SPARSE%GET+SPARSEKIT      h      GET+SPARSEKIT #   p  T   a   GET%THIS+SPARSEKIT     Ä  @   a   GET%I+SPARSEKIT       @   a   GET%J+SPARSEKIT (   D  T   a   SPARSE%GETNNZ+SPARSEKIT !     Z      GETNNZ+SPARSEKIT &   ò  T   a   GETNNZ%THIS+SPARSEKIT &   F  R   a   SPARSE%GETN+SPARSEKIT      Z      GETN+SPARSEKIT $   ò  T   a   GETN%THIS+SPARSEKIT &   F  R   a   SPARSE%GETA+SPARSEKIT      ö      GETA+SPARSEKIT $     T   a   GETA%THIS+SPARSEKIT '   â  S   a   SPARSE%GETAI+SPARSEKIT     5  h     GETAI+SPARSEKIT %     T   a   GETAI%THIS+SPARSEKIT '   ń  S   a   SPARSE%GETAJ+SPARSEKIT     D  ö      GETAJ+SPARSEKIT %   :  T   a   GETAJ%THIS+SPARSEKIT ,     X   a   SPARSE%PRINTVALUE+SPARSEKIT %   æ  n      PRINTVALUE+SPARSEKIT *   T  T   a   PRINTVALUE%THIS+SPARSEKIT '   š  @   a   PRINTVALUE%I+SPARSEKIT '   è  @   a   PRINTVALUE%J+SPARSEKIT .   (  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   t  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   Ï  `      PRINTNONZEROS+SPARSEKIT -   /   T   a   PRINTNONZEROS%THIS+SPARSEKIT 1      L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   Ï   V   a   SPARSE%PRINTALL+SPARSEKIT #   %!  `      PRINTALL+SPARSEKIT (   !  T   a   PRINTALL%THIS+SPARSEKIT ,   Ù!  L   a   PRINTALL%FILENAME+SPARSEKIT 1   %"  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   "  d      DELETEROWANDCOL+SPARSEKIT /   æ"  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   :#  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   z#  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   ș#  R   a   SPARSE%FREE+SPARSEKIT    $  R      FREE+SPARSEKIT $   ^$  T   a   FREE%THIS+SPARSEKIT C   Č$  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   %  R      HANDLEDUPLICATES+SPARSEKIT 0   b%  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT '   ¶%  }       PRECONDITIONERMETHODDT 8   3&  f   a   PRECONDITIONERMETHODDT%PRECONDITIONERDT 1   &  g      PRECONDITIONERDT+PRECONDITIONERM C    '  f   a   PRECONDITIONERDT%USEPRECONDITIONER+PRECONDITIONERM 9   f'        PRECONDITIONER_PROCEDURE+PRECONDITIONERM >   ç'  ^   a   PRECONDITIONER_PROCEDURE%THIS+PRECONDITIONERM @   E(     a   PRECONDITIONER_PROCEDURE%VECTOR+PRECONDITIONERM @   Ń(  T   a   PRECONDITIONER_PROCEDURE%MATRIX+PRECONDITIONERM B   %)     a   PRECONDITIONER_PROCEDURE%SOLUTION+PRECONDITIONERM =   ±)     a   PRECONDITIONER_PROCEDURE%ARG+PRECONDITIONERM 9   =*  T   a   PRECONDITIONERMETHODDT%USEPRECONDITIONER    *        METHOD    +  d   a   METHOD%THIS    v+     a   METHOD%VECTOR    ,  T   a   METHOD%MATRIX     V,     a   METHOD%SOLUTION    â,     a   METHOD%ARG 