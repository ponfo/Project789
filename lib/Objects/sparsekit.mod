  �2  �   k820309              19.1        ���^                                                                                                          
       src/lib/SparseKit.f90 SPARSEKIT       	       SPARSE i@ i@ i@ gen@SPARSE gen@TRANSPOSE gen@NORM gen@TRACE gen@ID                                                     
                                                           
                                                                  #SPARSE_SPARSE_PROD    #SPARSE_VECT_PROD    #COEF_SPARSE_PROD    &         @    @X                                                      #A    #B    #SPARSE              
  @                                                #SPARSE              
  @                                                #SPARSE    (        `   @X                                                                
    #MAT    #VECT 	   p          H r 
     7
S
O
 p        j            j                                      H r 
     7
S
O
 p        j            j                                                              
                                                   #SPARSE           0  
 @                              	                   
              &                                           &         @    @X                                                      #COEF    #MAT    #SPARSE              
                                      
                
                                                    #SPARSE                                                               #SPARSE_SPARSE_ADD    &         @    @X                                                      #A    #B    #SPARSE              
  @                                                #SPARSE              
                                                   #SPARSE                                                               #SPARSE_SPARSE_SUB    &         @    @X                                                      #A    #B    #SPARSE              
  @                                                #SPARSE              
                                                   #SPARSE                                                           u #CONSTRUCTOR    &         @   @                                                    #NNZ    #ROWS    #SPARSE              
  @                                                   
  @                                                  @    X                                    u #TRANSPOSE             @   X                                     u #NORM             @   X                                     u #TRACE             @    X                                     u #ID                   �  @               �                '                   #A    #AI    #AJ    #ROWCOUNTER    #N    #NNZ     #COUNTER !   #TRIPLET "   #ISCRSDONE '   #INIT (   #UPDATE -   #APPEND 2   #MAKECRS 8   #APPENDPOSTCRS <   #CHANGE B   #SETDIRICHLET H   #GET L   #GETNNZ Q   #GETN T   #GETA W   #GETAI Z   #GETAJ ]   #PRINTVALUE `   #PRINTNONZEROS f   #PRINTALL j   #DELETEROWANDCOL n   #FREE s   #HANDLEDUPLICATES v              � D                                                           
            &                                                     � D                                         H                             &                                                     � D                                         �                             &                                                      � D                                         �                             &                                                       � D                                                           � D                                   $                         � D                             !     (                         � D                              "     �       0             #TRIPLET #                 @  @              @           #     '�                    #A $   #ROW %   #COL &              �                              $                              
            &                                                      �                              %            H                             &                                                      �                              &            �                             &                                                        � D                              '           	      1         �   � $                     �      (             
     #INIT )   #         @     @                            )                    #THIS *   #NNZ +   #ROWS ,             
D                                *                   #SPARSE              
                                 +                     
                                 ,           1         �   � $                      �      -                  #UPDATE .   #         @     @                             .                    #THIS /   #NNZ 0   #ROWS 1             
D                                /                   #SPARSE              
                                 0                     
                                 1           1         �   � $                     �      2                  #APPEND 3   #         @     @                            3                    #THIS 4   #VAL 5   #ROW 6   #COL 7             
D @                              4                   #SPARSE              
  @                              5     
                
  @                              6                     
  @                              7           1         �   � $                     �      8                  #MAKECRS 9   #         @     @                            9                    #THIS :   #SORTROWS ;             
D @                              :                   #SPARSE              
 @                               ;           1         �   � $                     �      <                  #APPENDPOSTCRS =   #         @     @                            =                    #THIS >   #VAL ?   #ROW @   #COL A             
D                                >                   #SPARSE              
                                 ?     
                
                                 @                     
                                 A           1         �   � $                      �      B                  #CHANGE C   #         @     @                             C                    #THIS D   #VAL E   #ROW F   #COL G             
D                                D                   #SPARSE              
                                 E     
                
                                 F                     
                                 G           1         �   � $                      �      H                  #SETDIRICHLET I   #         @     @                             I                    #THIS J   #ROW K             
D                                J                   #SPARSE              
                                 K           1         �   � $                    �      L                  #GET M   %         @   @                           M                    
       #THIS N   #I O   #J P             
                                N                   #SPARSE              
                                 O                     
                                 P           1         �   � $                     �      Q              	    #GETNNZ R   %         @   @                            R                           #THIS S             
                                S                   #SPARSE    1         �   � $                     �      T              
    #GETN U   %         @   @                            U                           #THIS V             
                                V                   #SPARSE    1         �   � $                     �      W                  #GETA X   (        `   @                            X                                    
    #THIS Y   p          5 8 O#SPARSE     p        U             5 8 O#SPARSE     p        U                                    
                                Y                   #SPARSE    1         �   � $                     �      Z                  #GETAI [   (        `   @                            [                                        #THIS \   p           5 8 O#SPARSE     p        U        n                                         1     5 8 O#SPARSE     p        U        n                                          1                              
                                \                   #SPARSE    1         �   � $                     �      ]                  #GETAJ ^   (        `   @                            ^                                        #THIS _   p          5 8 O#SPARSE     p        U             5 8 O#SPARSE     p        U                                    
                                _                   #SPARSE    1         �   � $                      �      `                  #PRINTVALUE a   #         @     @                             a                    #THIS b   #I c   #J d   #FILENAME e             
D @                              b                   #SPARSE              
  @                              c                     
  @                              d                     
 @                             e                    1 1         �   � $                      �      f                  #PRINTNONZEROS g   #         @     @                             g                    #THIS h   #FILENAME i             
                                h                   #SPARSE              
 @                             i                    1 1         �   � $                      �      j                  #PRINTALL k   #         @     @                             k                    #THIS l   #FILENAME m             
D @                              l                   #SPARSE              
 @                             m                    1 1         �   � $                      �      n                  #DELETEROWANDCOL o   #         @     @                             o                    #THIS p   #ROW q   #COL r             
D                                p                   #SPARSE              
                                 q                     
                                 r           1         �   � $                     �      s                  #FREE t   #         @     @                            t                    #THIS u             
D                                u                   #SPARSE    1         �   � D                     �      v                  #HANDLEDUPLICATES w   #         @     @                            w                    #THIS x             
D                                x                   #SPARSE    &         @    X                                                     #A y   #SPARSE             
  @                              y                  #SPARSE    %         @   X                                                
       #A z             
                                 z                  #SPARSE    %         @   X                                                
       #A {             
                                 {                  #SPARSE    &         @    X                                                      #N |   #SPARSE              
  @                              |                         @                           
     SIZE    �   (      fn#fn    �   Y   b   uapp(SPARSEKIT    !  @   J  UTILITIESM    a  @   J  QUICKSORTM    �  �      i@ #   %  j      SPARSE_SPARSE_PROD %   �  T   a   SPARSE_SPARSE_PROD%A %   �  T   a   SPARSE_SPARSE_PROD%B !   7  w     SPARSE_VECT_PROD %   �  T   a   SPARSE_VECT_PROD%MAT &     �   a   SPARSE_VECT_PROD%VECT !   �  o      COEF_SPARSE_PROD &   �  @   a   COEF_SPARSE_PROD%COEF %   =  T   a   COEF_SPARSE_PROD%MAT    �  W      i@ "   �  j      SPARSE_SPARSE_ADD $   R  T   a   SPARSE_SPARSE_ADD%A $   �  T   a   SPARSE_SPARSE_ADD%B    �  W      i@ "   Q  j      SPARSE_SPARSE_SUB $   �  T   a   SPARSE_SPARSE_SUB%A $   	  T   a   SPARSE_SPARSE_SUB%B    c	  Q       gen@SPARSE    �	  o      CONSTRUCTOR     #
  @   a   CONSTRUCTOR%NNZ !   c
  @   a   CONSTRUCTOR%ROWS    �
  O       gen@TRANSPOSE    �
  J       gen@NORM    <  K       gen@TRACE    �  H       gen@ID    �  �      SPARSE    �  �   !   SPARSE%A      �   !   SPARSE%AI    �  �   !   SPARSE%AJ "   @  �   !   SPARSE%ROWCOUNTER    �  H   !   SPARSE%N      H   !   SPARSE%NNZ    d  H   !   SPARSE%COUNTER    �  ]   !   SPARSE%TRIPLET    	  i       TRIPLET    r  �   a   TRIPLET%A      �   a   TRIPLET%ROW    �  �   a   TRIPLET%COL !   .  H   !   SPARSE%ISCRSDONE    v  R   a   SPARSE%INIT    �  e      INIT    -  T   a   INIT%THIS    �  @   a   INIT%NNZ    �  @   a   INIT%ROWS      T   a   SPARSE%UPDATE    U  e      UPDATE    �  T   a   UPDATE%THIS      @   a   UPDATE%NNZ    N  @   a   UPDATE%ROWS    �  T   a   SPARSE%APPEND    �  m      APPEND    O  T   a   APPEND%THIS    �  @   a   APPEND%VAL    �  @   a   APPEND%ROW    #  @   a   APPEND%COL    c  U   a   SPARSE%MAKECRS    �  `      MAKECRS      T   a   MAKECRS%THIS !   l  @   a   MAKECRS%SORTROWS %   �  [   a   SPARSE%APPENDPOSTCRS      m      APPENDPOSTCRS #   t  T   a   APPENDPOSTCRS%THIS "   �  @   a   APPENDPOSTCRS%VAL "     @   a   APPENDPOSTCRS%ROW "   H  @   a   APPENDPOSTCRS%COL    �  T   a   SPARSE%CHANGE    �  m      CHANGE    I  T   a   CHANGE%THIS    �  @   a   CHANGE%VAL    �  @   a   CHANGE%ROW      @   a   CHANGE%COL $   ]  Z   a   SPARSE%SETDIRICHLET    �  [      SETDIRICHLET "     T   a   SETDIRICHLET%THIS !   f  @   a   SETDIRICHLET%ROW    �  Q   a   SPARSE%GET    �  h      GET    _  T   a   GET%THIS    �  @   a   GET%I    �  @   a   GET%J    3   T   a   SPARSE%GETNNZ    �   Z      GETNNZ    �   T   a   GETNNZ%THIS    5!  R   a   SPARSE%GETN    �!  Z      GETN    �!  T   a   GETN%THIS    5"  R   a   SPARSE%GETA    �"  �      GETA    }#  T   a   GETA%THIS    �#  S   a   SPARSE%GETAI    $$  h     GETAI    �%  T   a   GETAI%THIS    �%  S   a   SPARSE%GETAJ    3&  �      GETAJ    )'  T   a   GETAJ%THIS "   }'  X   a   SPARSE%PRINTVALUE    �'  n      PRINTVALUE     C(  T   a   PRINTVALUE%THIS    �(  @   a   PRINTVALUE%I    �(  @   a   PRINTVALUE%J $   )  L   a   PRINTVALUE%FILENAME %   c)  [   a   SPARSE%PRINTNONZEROS    �)  `      PRINTNONZEROS #   *  T   a   PRINTNONZEROS%THIS '   r*  L   a   PRINTNONZEROS%FILENAME     �*  V   a   SPARSE%PRINTALL    +  `      PRINTALL    t+  T   a   PRINTALL%THIS "   �+  L   a   PRINTALL%FILENAME '   ,  ]   a   SPARSE%DELETEROWANDCOL     q,  d      DELETEROWANDCOL %   �,  T   a   DELETEROWANDCOL%THIS $   )-  @   a   DELETEROWANDCOL%ROW $   i-  @   a   DELETEROWANDCOL%COL    �-  R   a   SPARSE%FREE    �-  R      FREE    M.  T   a   FREE%THIS (   �.  ^   !   SPARSE%HANDLEDUPLICATES !   �.  R      HANDLEDUPLICATES &   Q/  T   a   HANDLEDUPLICATES%THIS    �/  c       TRANSPOSE    0  T   a   TRANSPOSE%A    \0  W       NORM    �0  T   a   NORM%A    1  W       TRACE    ^1  T   a   TRACE%A    �1  c       ID    2  @   a   ID%N &   U2  =      SPARSE_VECT_PROD%SIZE 