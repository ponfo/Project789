  dZ  Č   k820309              19.1        5!­^                                                                                                          
       src/lib/Quicksort.f90 QUICKSORTM                                                     
                                                             u #SWAPINT    #SWAPREAL                                                          u #QUICKSORT1VECTI    #QUICKSORT1VECTR    #QUICKSORT2VECTII    #QUICKSORT2VECTIR    #QUICKSORT2VECTRR    #QUICKSORT2VECTRI 	   #QUICKSORT3VECTIII 
   #QUICKSORT3VECTIIR    #QUICKSORT3VECTIRI    #QUICKSORT3VECTIRR    #QUICKSORT3VECTRRR    #QUICKSORT3VECTRRI    #QUICKSORT3VECTRIR    #QUICKSORT3VECTRII    #QUICKSORT4VECTIIII    #QUICKSORT4VECTIIIR    #QUICKSORT4VECTIIRI    #QUICKSORT4VECTIIRR    #QUICKSORT4VECTIRII    #QUICKSORT4VECTIRIR    #QUICKSORT4VECTIRRI    #QUICKSORT4VECTIRRR    #QUICKSORT4VECTRIII    #QUICKSORT4VECTRIIR    #QUICKSORT4VECTRIRI    #QUICKSORT4VECTRIRR    #QUICKSORT4VECTRRII    #QUICKSORT4VECTRRIR    #QUICKSORT4VECTRRRI     #QUICKSORT4VECTRRRR !                                                "                                                       8                                             #                                                       4#         @      X                                                #A $   #B %             
D                                $                      
D                                %            #         @      X                                                #A &   #B '             
D                                &     
                 
D                                '     
       #         @     X                                                #VECT1 (   #LOW )   #HIGH *             
D @                              (                                  &                                                     
  @                              )                     
  @                              *           #         @     X                                                #VECT1 +   #LOW ,   #HIGH -             
D @                              +                   
               &                                                     
  @                              ,                     
  @                              -           #         @     X                                                #VECT1 .   #VECT2 /   #LOW 0   #HIGH 1             
D @                              .                                  &                                                     
D @                              /                                  &                                                     
  @                              0                     
  @                              1           #         @     X                                                #VECT1 2   #VECT2 3   #LOW 4   #HIGH 5             
D @                              2                                  &                                                     
D @                              3                   
               &                                                     
  @                              4                     
  @                              5           #         @     X                                                #VECT1 6   #VECT2 7   #LOW 8   #HIGH 9             
D @                              6                   
 	              &                                                     
D @                              7                   
 
              &                                                     
  @                              8                     
  @                              9           #         @     X                            	                    #VECT1 :   #VECT2 ;   #LOW <   #HIGH =             
D @                              :                   
               &                                                     
D @                              ;                                  &                                                     
  @                              <                     
  @                              =           #         @     X                            
                    #VECT1 >   #VECT2 ?   #VECT3 @   #LOW A   #HIGH B             
D @                              >                                  &                                                     
D @                              ?                                  &                                                     
D @                              @                                  &                                                     
  @                              A                     
  @                              B           #         @     X                                                #VECT1 C   #VECT2 D   #VECT3 E   #LOW F   #HIGH G             
D @                              C                                  &                                                     
D @                              D                                  &                                                     
D @                              E                   
               &                                                     
  @                              F                     
  @                              G           #         @     X                                                #VECT1 H   #VECT2 I   #VECT3 J   #LOW K   #HIGH L             
D @                              H                                  &                                                     
D @                              I                   
               &                                                     
D @                              J                                  &                                                     
  @                              K                     
  @                              L           #         @     X                                                #VECT1 M   #VECT2 N   #VECT3 O   #LOW P   #HIGH Q             
D @                              M                                  &                                                     
D @                              N                   
               &                                                     
D @                              O                   
               &                                                     
  @                              P                     
  @                              Q           #         @     X                                                #VECT1 R   #VECT2 S   #VECT3 T   #LOW U   #HIGH V             
D @                              R                   
               &                                                     
D @                              S                   
               &                                                     
D @                              T                   
               &                                                     
  @                              U                     
  @                              V           #         @     X                                                #VECT1 W   #VECT2 X   #VECT3 Y   #LOW Z   #HIGH [             
D @                              W                   
               &                                                     
D @                              X                   
               &                                                     
D @                              Y                                  &                                                     
  @                              Z                     
  @                              [           #         @     X                                                #VECT1 \   #VECT2 ]   #VECT3 ^   #LOW _   #HIGH `             
D @                              \                   
               &                                                     
D @                              ]                                  &                                                     
D @                              ^                   
               &                                                     
  @                              _                     
  @                              `           #         @     X                                                #VECT1 a   #VECT2 b   #VECT3 c   #LOW d   #HIGH e             
D @                              a                   
 "              &                                                     
D @                              b                                   &                                                     
D @                              c                    !              &                                                     
  @                              d                     
  @                              e           #         @     X                                                #VECT1 f   #VECT2 g   #VECT3 h   #VECT4 i   #LOW j   #HIGH k             
D @                              f                    #              &                                                     
D @                              g                    $              &                                                     
D @                              h                    %              &                                                     
D @                              i                    &              &                                                     
  @                              j                     
  @                              k           #         @     X                                                #VECT1 l   #VECT2 m   #VECT3 n   #VECT4 o   #LOW p   #HIGH q             
D @                              l                    '              &                                                     
D @                              m                    (              &                                                     
D @                              n                    )              &                                                     
D @                              o                   
 *              &                                                     
  @                              p                     
  @                              q           #         @     X                                                #VECT1 r   #VECT2 s   #VECT3 t   #VECT4 u   #LOW v   #HIGH w             
D @                              r                    +              &                                                     
D @                              s                    ,              &                                                     
D @                              t                   
 .              &                                                     
D @                              u                    -              &                                                     
  @                              v                     
  @                              w           #         @     X                                                #VECT1 x   #VECT2 y   #VECT3 z   #VECT4 {   #LOW |   #HIGH }             
D @                              x                    /              &                                                     
D @                              y                    0              &                                                     
D @                              z                   
 1              &                                                     
D @                              {                   
 2              &                                                     
  @                              |                     
  @                              }           #         @     X                                                #VECT1 ~   #VECT2    #VECT3    #VECT4    #LOW    #HIGH              
D @                              ~                    ?              &                                                     
D @                                                 
 B              &                                                     
D @                                                  @              &                                                     
D @                                                  A              &                                                     
  @                                                   
  @                                         #         @     X                                                #VECT1    #VECT2    #VECT3    #VECT4    #LOW    #HIGH              
D @                                                  ;              &                                                     
D @                                                 
 =              &                                                     
D @                                                  <              &                                                     
D @                                                 
 >              &                                                     
  @                                                   
  @                                         #         @     X                                                #VECT1    #VECT2    #VECT3    #VECT4    #LOW    #HIGH              
D @                                                  7              &                                                     
D @                                                 
 9              &                                                     
D @                                                 
 :              &                                                     
D @                                                  8              &                                                     
  @                                                   
  @                                         #         @     X                                                #VECT1    #VECT2    #VECT3    #VECT4    #LOW    #HIGH              
D @                                                  3              &                                                     
D @                                                 
 4              &                                                     
D @                                                 
 5              &                                                     
D @                                                 
 6              &                                                     
  @                                                   
  @                                         #         @     X                                                #VECT1    #VECT2    #VECT3    #VECT4    #LOW    #HIGH              
D @                                                 
 b              &                                                     
D @                                                  _              &                                                     
D @                                                  `              &                                                     
D @                                                  a              &                                                     
  @                                                   
  @                                         #         @     X                                                #VECT1    #VECT2    #VECT3    #VECT4    #LOW     #HIGH Ą             
D @                                                 
 ]              &                                                     
D @                                                  [              &                                                     
D @                                                  \              &                                                     
D @                                                 
 ^              &                                                     
  @                                                    
  @                              Ą           #         @     X                                                #VECT1 ˘   #VECT2 Ł   #VECT3 ¤   #VECT4 Ľ   #LOW Ś   #HIGH §             
D @                              ˘                   
 Y              &                                                     
D @                              Ł                    W              &                                                     
D @                              ¤                   
 Z              &                                                     
D @                              Ľ                    X              &                                                     
  @                              Ś                     
  @                              §           #         @     X                                                #VECT1 ¨   #VECT2 Š   #VECT3 Ş   #VECT4 Ť   #LOW Ź   #HIGH ­             
D @                              ¨                   
 T              &                                                     
D @                              Š                    S              &                                                     
D @                              Ş                   
 U              &                                                     
D @                              Ť                   
 V              &                                                     
  @                              Ź                     
  @                              ­           #         @     X                                                #VECT1 Ž   #VECT2 Ż   #VECT3 °   #VECT4 ą   #LOW ˛   #HIGH ł             
D @                              Ž                   
 Q              &                                                     
D @                              Ż                   
 R              &                                                     
D @                              °                    O              &                                                     
D @                              ą                    P              &                                                     
  @                              ˛                     
  @                              ł           #         @     X                                                #VECT1 ´   #VECT2 ľ   #VECT3 ś   #VECT4 ˇ   #LOW ¸   #HIGH š             
D @                              ´                   
 L              &                                                     
D @                              ľ                   
 M              &                                                     
D @                              ś                    K              &                                                     
D @                              ˇ                   
 N              &                                                     
  @                              ¸                     
  @                              š           #         @     X                                                 #VECT1 ş   #VECT2 ť   #VECT3 ź   #VECT4 ˝   #LOW ž   #HIGH ż             
D @                              ş                   
 H              &                                                     
D @                              ť                   
 I              &                                                     
D @                              ź                   
 J              &                                                     
D @                              ˝                    G              &                                                     
  @                              ž                     
  @                              ż           #         @     X                            !                    #VECT1 Ŕ   #VECT2 Á   #VECT3 Â   #VECT4 Ă   #LOW Ä   #HIGH Ĺ             
D @                              Ŕ                   
 C              &                                                     
D @                              Á                   
 D              &                                                     
D @                              Â                   
 E              &                                                     
D @                              Ă                   
 F              &                                                     
  @                              Ä                     
  @                              Ĺ                  )      fn#fn    É   @   J   UTILITIESM    	  [       gen@SWAP    d  ú      gen@QUICKSORT !   ^  q       RKIND+UTILITIESM !   Ď  q       IKIND+UTILITIESM    @  V       SWAPINT      @   a   SWAPINT%A    Ö  @   a   SWAPINT%B      V       SWAPREAL    l  @   a   SWAPREAL%A    Ź  @   a   SWAPREAL%B     ě  f       QUICKSORT1VECTI &   R     a   QUICKSORT1VECTI%VECT1 $   Ţ  @   a   QUICKSORT1VECTI%LOW %     @   a   QUICKSORT1VECTI%HIGH     ^  f       QUICKSORT1VECTR &   Ä     a   QUICKSORT1VECTR%VECT1 $   P	  @   a   QUICKSORT1VECTR%LOW %   	  @   a   QUICKSORT1VECTR%HIGH !   Đ	  q       QUICKSORT2VECTII '   A
     a   QUICKSORT2VECTII%VECT1 '   Í
     a   QUICKSORT2VECTII%VECT2 %   Y  @   a   QUICKSORT2VECTII%LOW &     @   a   QUICKSORT2VECTII%HIGH !   Ů  q       QUICKSORT2VECTIR '   J     a   QUICKSORT2VECTIR%VECT1 '   Ö     a   QUICKSORT2VECTIR%VECT2 %   b  @   a   QUICKSORT2VECTIR%LOW &   ˘  @   a   QUICKSORT2VECTIR%HIGH !   â  q       QUICKSORT2VECTRR '   S     a   QUICKSORT2VECTRR%VECT1 '   ß     a   QUICKSORT2VECTRR%VECT2 %   k  @   a   QUICKSORT2VECTRR%LOW &   Ť  @   a   QUICKSORT2VECTRR%HIGH !   ë  q       QUICKSORT2VECTRI '   \     a   QUICKSORT2VECTRI%VECT1 '   č     a   QUICKSORT2VECTRI%VECT2 %   t  @   a   QUICKSORT2VECTRI%LOW &   ´  @   a   QUICKSORT2VECTRI%HIGH "   ô  |       QUICKSORT3VECTIII (   p     a   QUICKSORT3VECTIII%VECT1 (   ü     a   QUICKSORT3VECTIII%VECT2 (        a   QUICKSORT3VECTIII%VECT3 &     @   a   QUICKSORT3VECTIII%LOW '   T  @   a   QUICKSORT3VECTIII%HIGH "     |       QUICKSORT3VECTIIR (        a   QUICKSORT3VECTIIR%VECT1 (        a   QUICKSORT3VECTIIR%VECT2 (   (     a   QUICKSORT3VECTIIR%VECT3 &   ´  @   a   QUICKSORT3VECTIIR%LOW '   ô  @   a   QUICKSORT3VECTIIR%HIGH "   4  |       QUICKSORT3VECTIRI (   °     a   QUICKSORT3VECTIRI%VECT1 (   <     a   QUICKSORT3VECTIRI%VECT2 (   Č     a   QUICKSORT3VECTIRI%VECT3 &   T  @   a   QUICKSORT3VECTIRI%LOW '     @   a   QUICKSORT3VECTIRI%HIGH "   Ô  |       QUICKSORT3VECTIRR (   P     a   QUICKSORT3VECTIRR%VECT1 (   Ü     a   QUICKSORT3VECTIRR%VECT2 (   h     a   QUICKSORT3VECTIRR%VECT3 &   ô  @   a   QUICKSORT3VECTIRR%LOW '   4  @   a   QUICKSORT3VECTIRR%HIGH "   t  |       QUICKSORT3VECTRRR (   đ     a   QUICKSORT3VECTRRR%VECT1 (   |     a   QUICKSORT3VECTRRR%VECT2 (        a   QUICKSORT3VECTRRR%VECT3 &     @   a   QUICKSORT3VECTRRR%LOW '   Ô  @   a   QUICKSORT3VECTRRR%HIGH "     |       QUICKSORT3VECTRRI (        a   QUICKSORT3VECTRRI%VECT1 (         a   QUICKSORT3VECTRRI%VECT2 (   ¨      a   QUICKSORT3VECTRRI%VECT3 &   4!  @   a   QUICKSORT3VECTRRI%LOW '   t!  @   a   QUICKSORT3VECTRRI%HIGH "   ´!  |       QUICKSORT3VECTRIR (   0"     a   QUICKSORT3VECTRIR%VECT1 (   ź"     a   QUICKSORT3VECTRIR%VECT2 (   H#     a   QUICKSORT3VECTRIR%VECT3 &   Ô#  @   a   QUICKSORT3VECTRIR%LOW '   $  @   a   QUICKSORT3VECTRIR%HIGH "   T$  |       QUICKSORT3VECTRII (   Đ$     a   QUICKSORT3VECTRII%VECT1 (   \%     a   QUICKSORT3VECTRII%VECT2 (   č%     a   QUICKSORT3VECTRII%VECT3 &   t&  @   a   QUICKSORT3VECTRII%LOW '   ´&  @   a   QUICKSORT3VECTRII%HIGH #   ô&         QUICKSORT4VECTIIII )   {'     a   QUICKSORT4VECTIIII%VECT1 )   (     a   QUICKSORT4VECTIIII%VECT2 )   (     a   QUICKSORT4VECTIIII%VECT3 )   )     a   QUICKSORT4VECTIIII%VECT4 '   Ť)  @   a   QUICKSORT4VECTIIII%LOW (   ë)  @   a   QUICKSORT4VECTIIII%HIGH #   +*         QUICKSORT4VECTIIIR )   ˛*     a   QUICKSORT4VECTIIIR%VECT1 )   >+     a   QUICKSORT4VECTIIIR%VECT2 )   Ę+     a   QUICKSORT4VECTIIIR%VECT3 )   V,     a   QUICKSORT4VECTIIIR%VECT4 '   â,  @   a   QUICKSORT4VECTIIIR%LOW (   "-  @   a   QUICKSORT4VECTIIIR%HIGH #   b-         QUICKSORT4VECTIIRI )   é-     a   QUICKSORT4VECTIIRI%VECT1 )   u.     a   QUICKSORT4VECTIIRI%VECT2 )   /     a   QUICKSORT4VECTIIRI%VECT3 )   /     a   QUICKSORT4VECTIIRI%VECT4 '   0  @   a   QUICKSORT4VECTIIRI%LOW (   Y0  @   a   QUICKSORT4VECTIIRI%HIGH #   0         QUICKSORT4VECTIIRR )    1     a   QUICKSORT4VECTIIRR%VECT1 )   Ź1     a   QUICKSORT4VECTIIRR%VECT2 )   82     a   QUICKSORT4VECTIIRR%VECT3 )   Ä2     a   QUICKSORT4VECTIIRR%VECT4 '   P3  @   a   QUICKSORT4VECTIIRR%LOW (   3  @   a   QUICKSORT4VECTIIRR%HIGH #   Đ3         QUICKSORT4VECTIRII )   W4     a   QUICKSORT4VECTIRII%VECT1 )   ă4     a   QUICKSORT4VECTIRII%VECT2 )   o5     a   QUICKSORT4VECTIRII%VECT3 )   ű5     a   QUICKSORT4VECTIRII%VECT4 '   6  @   a   QUICKSORT4VECTIRII%LOW (   Ç6  @   a   QUICKSORT4VECTIRII%HIGH #   7         QUICKSORT4VECTIRIR )   7     a   QUICKSORT4VECTIRIR%VECT1 )   8     a   QUICKSORT4VECTIRIR%VECT2 )   Ś8     a   QUICKSORT4VECTIRIR%VECT3 )   29     a   QUICKSORT4VECTIRIR%VECT4 '   ž9  @   a   QUICKSORT4VECTIRIR%LOW (   ţ9  @   a   QUICKSORT4VECTIRIR%HIGH #   >:         QUICKSORT4VECTIRRI )   Ĺ:     a   QUICKSORT4VECTIRRI%VECT1 )   Q;     a   QUICKSORT4VECTIRRI%VECT2 )   Ý;     a   QUICKSORT4VECTIRRI%VECT3 )   i<     a   QUICKSORT4VECTIRRI%VECT4 '   ő<  @   a   QUICKSORT4VECTIRRI%LOW (   5=  @   a   QUICKSORT4VECTIRRI%HIGH #   u=         QUICKSORT4VECTIRRR )   ü=     a   QUICKSORT4VECTIRRR%VECT1 )   >     a   QUICKSORT4VECTIRRR%VECT2 )   ?     a   QUICKSORT4VECTIRRR%VECT3 )    ?     a   QUICKSORT4VECTIRRR%VECT4 '   ,@  @   a   QUICKSORT4VECTIRRR%LOW (   l@  @   a   QUICKSORT4VECTIRRR%HIGH #   Ź@         QUICKSORT4VECTRIII )   3A     a   QUICKSORT4VECTRIII%VECT1 )   żA     a   QUICKSORT4VECTRIII%VECT2 )   KB     a   QUICKSORT4VECTRIII%VECT3 )   ×B     a   QUICKSORT4VECTRIII%VECT4 '   cC  @   a   QUICKSORT4VECTRIII%LOW (   ŁC  @   a   QUICKSORT4VECTRIII%HIGH #   ăC         QUICKSORT4VECTRIIR )   jD     a   QUICKSORT4VECTRIIR%VECT1 )   öD     a   QUICKSORT4VECTRIIR%VECT2 )   E     a   QUICKSORT4VECTRIIR%VECT3 )   F     a   QUICKSORT4VECTRIIR%VECT4 '   F  @   a   QUICKSORT4VECTRIIR%LOW (   ÚF  @   a   QUICKSORT4VECTRIIR%HIGH #   G         QUICKSORT4VECTRIRI )   ĄG     a   QUICKSORT4VECTRIRI%VECT1 )   -H     a   QUICKSORT4VECTRIRI%VECT2 )   šH     a   QUICKSORT4VECTRIRI%VECT3 )   EI     a   QUICKSORT4VECTRIRI%VECT4 '   ŃI  @   a   QUICKSORT4VECTRIRI%LOW (   J  @   a   QUICKSORT4VECTRIRI%HIGH #   QJ         QUICKSORT4VECTRIRR )   ŘJ     a   QUICKSORT4VECTRIRR%VECT1 )   dK     a   QUICKSORT4VECTRIRR%VECT2 )   đK     a   QUICKSORT4VECTRIRR%VECT3 )   |L     a   QUICKSORT4VECTRIRR%VECT4 '   M  @   a   QUICKSORT4VECTRIRR%LOW (   HM  @   a   QUICKSORT4VECTRIRR%HIGH #   M         QUICKSORT4VECTRRII )   N     a   QUICKSORT4VECTRRII%VECT1 )   N     a   QUICKSORT4VECTRRII%VECT2 )   'O     a   QUICKSORT4VECTRRII%VECT3 )   łO     a   QUICKSORT4VECTRRII%VECT4 '   ?P  @   a   QUICKSORT4VECTRRII%LOW (   P  @   a   QUICKSORT4VECTRRII%HIGH #   żP         QUICKSORT4VECTRRIR )   FQ     a   QUICKSORT4VECTRRIR%VECT1 )   ŇQ     a   QUICKSORT4VECTRRIR%VECT2 )   ^R     a   QUICKSORT4VECTRRIR%VECT3 )   ęR     a   QUICKSORT4VECTRRIR%VECT4 '   vS  @   a   QUICKSORT4VECTRRIR%LOW (   śS  @   a   QUICKSORT4VECTRRIR%HIGH #   öS         QUICKSORT4VECTRRRI )   }T     a   QUICKSORT4VECTRRRI%VECT1 )   	U     a   QUICKSORT4VECTRRRI%VECT2 )   U     a   QUICKSORT4VECTRRRI%VECT3 )   !V     a   QUICKSORT4VECTRRRI%VECT4 '   ­V  @   a   QUICKSORT4VECTRRRI%LOW (   íV  @   a   QUICKSORT4VECTRRRI%HIGH #   -W         QUICKSORT4VECTRRRR )   ´W     a   QUICKSORT4VECTRRRR%VECT1 )   @X     a   QUICKSORT4VECTRRRR%VECT2 )   ĚX     a   QUICKSORT4VECTRRRR%VECT3 )   XY     a   QUICKSORT4VECTRRRR%VECT4 '   äY  @   a   QUICKSORT4VECTRRRR%LOW (   $Z  @   a   QUICKSORT4VECTRRRR%HIGH 