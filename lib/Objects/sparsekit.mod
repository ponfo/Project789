  Ķ5     k820309              19.1        nę^                                                                                                          
       src/lib/SparseKit.f90 SPARSEKIT              SPARSE i@ i@ i@ gen@SPARSE gen@TRANSPOSE gen@NORM gen@TRACE gen@ID gen@INVERSE gen@INVERSELUMPED                                                     
                                                           
                                                           
                                                                  #SPARSE_SPARSE_PROD    #SPARSE_VECT_PROD    #COEF_SPARSE_PROD    &         @    @X                                                      #A    #B    #SPARSE              
  @                                                #SPARSE              
  @                                                #SPARSE    (        `   @X                                                                
    #MAT 	   #VECT 
   p          H r      7
S
O
 p        j            j                                      H r      7
S
O
 p        j            j                                                              
                                 	                  #SPARSE           0  
 @                              
                   
              &                                           &         @    @X                                                      #COEF    #MAT    #SPARSE              
                                      
                
                                                    #SPARSE                                                               #SPARSE_SPARSE_ADD    &         @    @X                                                      #A    #B    #SPARSE              
  @                                                #SPARSE              
                                                   #SPARSE                                                               #SPARSE_SPARSE_SUB    &         @    @X                                                      #A    #B    #SPARSE              
  @                                                #SPARSE              
                                                   #SPARSE                                                           u #CONSTRUCTOR    &         @   @                                                    #NNZ    #ROWS    #SPARSE              
  @                                                   
  @                                                  @    X                                    u #TRANSPOSE             @   X                                     u #NORM             @   X                                     u #TRACE             @    X                                     u #ID             @   X                                     u #INVERSE             @   X                                     u #INVERSELUMPED                   @  @                                '                    #DUMMY                                                                                    @               Ā                '                   #A     #AI !   #AJ "   #ROWCOUNTER #   #N $   #NNZ %   #COUNTER &   #TRIPLET '   #ISCRSDONE ,   #INIT -   #UPDATE 2   #APPEND 7   #MAKECRS =   #APPENDPOSTCRS A   #CHANGE G   #SETDIRICHLET M   #GET Q   #GETNNZ V   #GETN Y   #GETA \   #GETAI _   #GETAJ b   #PRINTVALUE e   #PRINTNONZEROS k   #PRINTALL o   #DELETEROWANDCOL s   #FREE x   #HANDLEDUPLICATES {               D                                                            
            &                                                      D                             !            H                             &                                                      D                             "                                         &                                                       D                             #            Ø                             &                                                        D                             $                               D                             %     $                          D                             &     (                          D                              '     Ø       0             #TRIPLET (                 @  @              @           (     'Ø                    #A )   #ROW *   #COL +                                            )                              
            &                                                                                    *            H                             &                                                                                    +                                         &                                                         D                              ,           	      1         Ā    $                           -             
     #INIT .   #         @     @                            .                    #THIS /   #NNZ 0   #ROWS 1             
D                                /                   #SPARSE              
                                 0                     
                                 1           1         Ā    $                            2                  #UPDATE 3   #         @     @                             3                    #THIS 4   #NNZ 5   #ROWS 6             
D                                4                   #SPARSE              
                                 5                     
                                 6           1         Ā    $                           7                  #APPEND 8   #         @     @                            8                    #THIS 9   #VAL :   #ROW ;   #COL <             
D @                              9                   #SPARSE              
  @                              :     
                
  @                              ;                     
  @                              <           1         Ā    $                           =                  #MAKECRS >   #         @     @                            >                    #THIS ?   #SORTROWS @             
D @                              ?                   #SPARSE              
 @                               @           1         Ā    $                           A                  #APPENDPOSTCRS B   #         @     @                            B                    #THIS C   #VAL D   #ROW E   #COL F             
D                                C                   #SPARSE              
                                 D     
                
                                 E                     
                                 F           1         Ā    $                            G                  #CHANGE H   #         @     @                             H                    #THIS I   #VAL J   #ROW K   #COL L             
D                                I                   #SPARSE              
                                 J     
                
                                 K                     
                                 L           1         Ā    $                            M                  #SETDIRICHLET N   #         @     @                             N                    #THIS O   #ROW P             
D                                O                   #SPARSE              
                                 P           1         Ā    $                          Q                  #GET R   %         @   @                          R                    
       #THIS S   #I T   #J U             
                                 S                  #SPARSE              
                                 T                     
                                 U           1         Ā    $                           V              	    #GETNNZ W   %         @   @                           W                           #THIS X             
                                 X                  #SPARSE    1         Ā    $                          Y              
    #GETN Z   %         @   @                          Z                           #THIS [             
                                 [                  #SPARSE    1         Ā    $                          \                  #GETA ]   (        `   @                          ]                                    
    #THIS ^   p          5 8 O#SPARSE     p        U     %       5 8 O#SPARSE     p        U     %                              
                                 ^                  #SPARSE    1         Ā    $                          _                  #GETAI `   (        `   @                          `                                        #THIS a   p           5 8 O#SPARSE     p        U     $   n                                         1     5 8 O#SPARSE     p        U     $   n                                          1                              
                                 a                  #SPARSE    1         Ā    $                          b                  #GETAJ c   (        `   @                          c                                        #THIS d   p          5 8 O#SPARSE     p        U     %       5 8 O#SPARSE     p        U     %                              
                                 d                  #SPARSE    1         Ā    $                            e                  #PRINTVALUE f   #         @     @                             f                    #THIS g   #I h   #J i   #FILENAME j             
                                g                   #SPARSE              
  @                              h                     
  @                              i                     
 @                             j                    1 1         Ā    $                            k                  #PRINTNONZEROS l   #         @     @                             l                    #THIS m   #FILENAME n             
                                m                   #SPARSE              
 @                             n                    1 1         Ā    $                            o                  #PRINTALL p   #         @     @                             p                    #THIS q   #FILENAME r             
                                q                   #SPARSE              
 @                             r                    1 1         Ā    $                            s                  #DELETEROWANDCOL t   #         @     @                             t                    #THIS u   #ROW v   #COL w             
D                                u                   #SPARSE              
                                 v                     
                                 w           1         Ā    $                           x                  #FREE y   #         @     @                            y                    #THIS z             
D                                z                   #SPARSE    1         Ā    D                           {                  #HANDLEDUPLICATES |   #         @     @                            |                    #THIS }             
D @                              }                   #SPARSE    &         @    X                                                     #A ~   #SPARSE             
  @                              ~                  #SPARSE    %         @   X                                                
       #A              
                                                   #SPARSE    %         @   X                                                
       #A              
                                                   #SPARSE    &         @    X                                                      #N    #SPARSE              
  @                                         &         @   X                                                      #A    #SPARSE              
                                                   #SPARSE    &         @   X                                                      #MATRIX    #SPARSE              
                                                   #SPARSE                  @                                SIZE        (      fn#fn    Č   w   b   uapp(SPARSEKIT    ?  @   J  UTILITIESM      @   J  QUICKSORTM    ŋ  @   J  MKL_PARDISO    ĸ        i@ #     j      SPARSE_SPARSE_PROD %   í  T   a   SPARSE_SPARSE_PROD%A %   A  T   a   SPARSE_SPARSE_PROD%B !     w     SPARSE_VECT_PROD %     T   a   SPARSE_VECT_PROD%MAT &   `     a   SPARSE_VECT_PROD%VECT !   ė  o      COEF_SPARSE_PROD &   [  @   a   COEF_SPARSE_PROD%COEF %     T   a   COEF_SPARSE_PROD%MAT    ï  W      i@ "   F  j      SPARSE_SPARSE_ADD $   °  T   a   SPARSE_SPARSE_ADD%A $     T   a   SPARSE_SPARSE_ADD%B    X  W      i@ "   Ŋ  j      SPARSE_SPARSE_SUB $   	  T   a   SPARSE_SPARSE_SUB%A $   m	  T   a   SPARSE_SPARSE_SUB%B    Á	  Q       gen@SPARSE    
  o      CONSTRUCTOR     
  @   a   CONSTRUCTOR%NNZ !   Á
  @   a   CONSTRUCTOR%ROWS      O       gen@TRANSPOSE    P  J       gen@NORM      K       gen@TRACE    å  H       gen@ID    -  M       gen@INVERSE "   z  S       gen@INVERSELUMPED 7   Í  [      MKL_PARDISO_HANDLE+MKL_PARDISO_PRIVATE =   (  H   a   MKL_PARDISO_HANDLE%DUMMY+MKL_PARDISO_PRIVATE    p  ĩ      SPARSE    %     !   SPARSE%A    đ     !   SPARSE%AI    M     !   SPARSE%AJ "   á     !   SPARSE%ROWCOUNTER    u  H   !   SPARSE%N    ―  H   !   SPARSE%NNZ      H   !   SPARSE%COUNTER    M  ]   !   SPARSE%TRIPLET    Š  i       TRIPLET         a   TRIPLET%A    §     a   TRIPLET%ROW    ;     a   TRIPLET%COL !   Ï  H   !   SPARSE%ISCRSDONE      R   a   SPARSE%INIT    i  e      INIT    Î  T   a   INIT%THIS    "  @   a   INIT%NNZ    b  @   a   INIT%ROWS    Ē  T   a   SPARSE%UPDATE    ö  e      UPDATE    [  T   a   UPDATE%THIS    Ŋ  @   a   UPDATE%NNZ    ï  @   a   UPDATE%ROWS    /  T   a   SPARSE%APPEND      m      APPEND    ð  T   a   APPEND%THIS    D  @   a   APPEND%VAL      @   a   APPEND%ROW    Ä  @   a   APPEND%COL      U   a   SPARSE%MAKECRS    Y  `      MAKECRS    đ  T   a   MAKECRS%THIS !     @   a   MAKECRS%SORTROWS %   M  [   a   SPARSE%APPENDPOSTCRS    Ļ  m      APPENDPOSTCRS #     T   a   APPENDPOSTCRS%THIS "   i  @   a   APPENDPOSTCRS%VAL "   Đ  @   a   APPENDPOSTCRS%ROW "   é  @   a   APPENDPOSTCRS%COL    )  T   a   SPARSE%CHANGE    }  m      CHANGE    ę  T   a   CHANGE%THIS    >  @   a   CHANGE%VAL    ~  @   a   CHANGE%ROW    ū  @   a   CHANGE%COL $   þ  Z   a   SPARSE%SETDIRICHLET    X  [      SETDIRICHLET "   ģ  T   a   SETDIRICHLET%THIS !      @   a   SETDIRICHLET%ROW    G   Q   a   SPARSE%GET       h      GET     !  T   a   GET%THIS    T!  @   a   GET%I    !  @   a   GET%J    Ô!  T   a   SPARSE%GETNNZ    ("  Z      GETNNZ    "  T   a   GETNNZ%THIS    Ö"  R   a   SPARSE%GETN    (#  Z      GETN    #  T   a   GETN%THIS    Ö#  R   a   SPARSE%GETA    ($  ö      GETA    %  T   a   GETA%THIS    r%  S   a   SPARSE%GETAI    Å%  h     GETAI    -'  T   a   GETAI%THIS    '  S   a   SPARSE%GETAJ    Ô'  ö      GETAJ    Ę(  T   a   GETAJ%THIS "   )  X   a   SPARSE%PRINTVALUE    v)  n      PRINTVALUE     ä)  T   a   PRINTVALUE%THIS    8*  @   a   PRINTVALUE%I    x*  @   a   PRINTVALUE%J $   ļ*  L   a   PRINTVALUE%FILENAME %   +  [   a   SPARSE%PRINTNONZEROS    _+  `      PRINTNONZEROS #   ŋ+  T   a   PRINTNONZEROS%THIS '   ,  L   a   PRINTNONZEROS%FILENAME     _,  V   a   SPARSE%PRINTALL    ĩ,  `      PRINTALL    -  T   a   PRINTALL%THIS "   i-  L   a   PRINTALL%FILENAME '   ĩ-  ]   a   SPARSE%DELETEROWANDCOL     .  d      DELETEROWANDCOL %   v.  T   a   DELETEROWANDCOL%THIS $   Ę.  @   a   DELETEROWANDCOL%ROW $   
/  @   a   DELETEROWANDCOL%COL    J/  R   a   SPARSE%FREE    /  R      FREE    î/  T   a   FREE%THIS (   B0  ^   !   SPARSE%HANDLEDUPLICATES !    0  R      HANDLEDUPLICATES &   ō0  T   a   HANDLEDUPLICATES%THIS    F1  c       TRANSPOSE    Đ1  T   a   TRANSPOSE%A    ý1  W       NORM    T2  T   a   NORM%A    Ļ2  W       TRACE    ĸ2  T   a   TRACE%A    S3  c       ID    ķ3  @   a   ID%N    ö3  c       INVERSE    Y4  T   a   INVERSE%A    ­4  h       INVERSELUMPED %   5  T   a   INVERSELUMPED%MATRIX &   i5  =      SPARSE_VECT_PROD%SIZE 