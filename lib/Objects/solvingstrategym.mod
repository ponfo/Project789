  f  B   k820309              19.1        6��^                                                                                                          
       src/SolvingStrategy/SolvingStrategy.f90 SOLVINGSTRATEGYM              NEWSOLVINGSTRATEGYDT SOLVINGSTRATEGYDT gen@SETSOLVINGSTRATEGY                      @                              
                            @                              
                            @                              
                                                              u #CONSTRUCTOR    &         @   @                                �                      #STRATEGY    #SOLVINGSTRATEGYDT              
  @                                   �             #NEWSOLVINGSTRATEGYDT                      @               �             '�                   #PROCESSDT    #BUILDERANDSOLVER    #SCHEME &                � $                                   �                      #PROCESSDT 	                 @  @               @          	     '�                    #PROCESS 
   #INIT    #CHANGE    #RUN                � $                             
                            #NEWPROCESSDT                   @  @                               '                      #USEPROCESS    1         �   � $                      �                       #NEWPROCESS_PROCEDURE    #         @     @                                	               #THIS              
                                                    #NEWPROCESSDT    1         �   � $                      �                        #INIT    #         @     @                                                #THIS    #PROCESS              
                                     �               #PROCESSDT 	             
                                                     #NEWPROCESSDT    1         �   � $                      �                        #CHANGE    #         @     @                                                #THIS    #NEWPROCESS              
                                     �               #PROCESSDT 	             
                                                     #NEWPROCESSDT    1         �   � $                      �                        #RUN    #         @     @                                                #THIS              
                                     �               #PROCESSDT 	               � $                                  �       �              #BUILDERANDSOLVERDT                   @  @               @               '�                    #BUILDERANDSOLVER    #INIT    #CHANGE "               � $                                                         #NEWBUILDERANDSOLVERDT                   @  @                               '                        1         �   � $                      �                        #INIT    #         @     @                                                #THIS     #BUILDERANDSOLVER !             
                                      �               #BUILDERANDSOLVERDT              
                                 !                    #NEWBUILDERANDSOLVERDT    1         �   � $                      �      "                  #CHANGE #   #         @     @                            #                    #THIS $   #NEWBUILDERANDSOLVER %             
                                $     �               #BUILDERANDSOLVERDT              
                                 %                    #NEWBUILDERANDSOLVERDT                � $                             &     �                     #SCHEMEDT '                  @  @               @          '     '�                    #SCHEME (   #INIT /   #CHANGE 3               � $                             (                            #NEWSCHEMEDT )                  @  @                          )     '                      #INTEGRATE *   1         �   � $                      $�     *                   #INTEGRATOR_INTERFACE +   #         @     @                           +     	               #THIS ,   #DT -   #MULTI_STEP .             
                               ,                     #NEWPROCESSDT              
                                -     
                
                                 .           1         �   � $                      �      /                  #INIT 0   #         @     @                            0                    #THIS 1   #SCHEME 2             
                                1     �               #SCHEMEDT '             
                                 2                    #NEWSCHEMEDT )   1         �   � $                      �      3                  #CHANGE 4   #         @     @                            4                    #THIS 5   #NEWSCHEME 6             
                                5     �               #SCHEMEDT '             
                                 6                    #NEWSCHEMEDT )                     @               �               '�                    #STRATEGY 7   #INIT 8   #CHANGE <               � $                             7     �                     #NEWSOLVINGSTRATEGYDT    1         �   � $                     �      8                  #INIT 9   #         @     @                            9                    #THIS :   #STRATEGY ;             
D                                :     �               #SOLVINGSTRATEGYDT              
                                 ;     �             #NEWSOLVINGSTRATEGYDT    1         �   � $                      �      <                  #CHANGE =   #         @     @                             =                    #THIS >   #NEWSTRATEGY ?             
D                                >     �               #SOLVINGSTRATEGYDT              
                                 ?     �             #NEWSOLVINGSTRATEGYDT       �   A      fn#fn &   �   N   b   uapp(SOLVINGSTRATEGYM    /  @   J  PROCESSM "   o  @   J  BUILDERANDSOLVERM    �  @   J  SCHEMEM '   �  Q       gen@SETSOLVINGSTRATEGY    @  u      CONSTRUCTOR %   �  b   a   CONSTRUCTOR%STRATEGY %     �       NEWSOLVINGSTRATEGYDT /   �  _   a   NEWSOLVINGSTRATEGYDT%PROCESSDT #   �  |      PROCESSDT+PROCESSM +   s  b   a   PROCESSDT%PROCESS+PROCESSM &   �  `      NEWPROCESSDT+PROCESSM 1   5  b   a   NEWPROCESSDT%USEPROCESS+PROCESSM .   �  R      NEWPROCESS_PROCEDURE+PROCESSM 3   �  Z   a   NEWPROCESS_PROCEDURE%THIS+PROCESSM (   C  R   a   PROCESSDT%INIT+PROCESSM    �  _      INIT+PROCESSM #   �  W   a   INIT%THIS+PROCESSM &   K  Z   a   INIT%PROCESS+PROCESSM *   �  T   a   PROCESSDT%CHANGE+PROCESSM     �  b      CHANGE+PROCESSM %   [  W   a   CHANGE%THIS+PROCESSM +   �  Z   a   CHANGE%NEWPROCESS+PROCESSM '   	  Q   a   PROCESSDT%RUN+PROCESSM    ]	  R      RUN+PROCESSM "   �	  W   a   RUN%THIS+PROCESSM 6   
  h   a   NEWSOLVINGSTRATEGYDT%BUILDERANDSOLVER 5   n
  |      BUILDERANDSOLVERDT+BUILDERANDSOLVERM F   �
  k   a   BUILDERANDSOLVERDT%BUILDERANDSOLVER+BUILDERANDSOLVERM 8   U  P      NEWBUILDERANDSOLVERDT+BUILDERANDSOLVERM :   �  R   a   BUILDERANDSOLVERDT%INIT+BUILDERANDSOLVERM '   �  h      INIT+BUILDERANDSOLVERM ,   _  `   a   INIT%THIS+BUILDERANDSOLVERM 8   �  c   a   INIT%BUILDERANDSOLVER+BUILDERANDSOLVERM <   "  T   a   BUILDERANDSOLVERDT%CHANGE+BUILDERANDSOLVERM )   v  k      CHANGE+BUILDERANDSOLVERM .   �  `   a   CHANGE%THIS+BUILDERANDSOLVERM =   A  c   a   CHANGE%NEWBUILDERANDSOLVER+BUILDERANDSOLVERM ,   �  ^   a   NEWSOLVINGSTRATEGYDT%SCHEME !     r      SCHEMEDT+SCHEMEM (   t  a   a   SCHEMEDT%SCHEME+SCHEMEM $   �  _      NEWSCHEMEDT+SCHEMEM .   4  b   a   NEWSCHEMEDT%INTEGRATE+SCHEMEM -   �  j      INTEGRATOR_INTERFACE+SCHEMEM 2      Z   a   INTEGRATOR_INTERFACE%THIS+SCHEMEM 0   Z  @   a   INTEGRATOR_INTERFACE%DT+SCHEMEM 8   �  @   a   INTEGRATOR_INTERFACE%MULTI_STEP+SCHEMEM &   �  R   a   SCHEMEDT%INIT+SCHEMEM    ,  ^      INIT+SCHEMEM "   �  V   a   INIT%THIS+SCHEMEM $   �  Y   a   INIT%SCHEME+SCHEMEM (   9  T   a   SCHEMEDT%CHANGE+SCHEMEM    �  a      CHANGE+SCHEMEM $   �  V   a   CHANGE%THIS+SCHEMEM )   D  Y   a   CHANGE%NEWSCHEME+SCHEMEM "   �  t       SOLVINGSTRATEGYDT +     j   a   SOLVINGSTRATEGYDT%STRATEGY '   {  R   a   SOLVINGSTRATEGYDT%INIT    �  `      INIT    -  _   a   INIT%THIS    �  b   a   INIT%STRATEGY )   �  T   a   SOLVINGSTRATEGYDT%CHANGE    B  c      CHANGE    �  _   a   CHANGE%THIS #     b   a   CHANGE%NEWSTRATEGY 