  �^  �   k820309              19.1        ���^                                                                                                          
       src/condition/Condition.f90 CONDITIONM              CONDITIONDT                                                     
                       @�                                  
                            @                              
                            @                              
                            @                              
                      @  @                              '�                    #PTR                 �$                                  �                      #NODEDT                      @               �               '�                    #POINTDT 	   #DOF >   #SOURCE O   #INITNODE1D �   #INITNODE2D �   #INITNODE3D �   #ASSIGNSOURCE �   #ASSIGNDOF �   #FIXDOF �   #FREEDOF �   #GETNDOF �                � $                              	     P                      #POINTDT 
                 @  @               D           
     'P                    #ID    #COORD    #INITPOINT1D    #INITPOINT2D    #INITPOINT3D    #SETID    #SETX #   #SETY '   #SETZ +   #GETID /   #GETX 2   #GETY 5   #GETZ 8   #GETDIMENSION ;                � $                                                          � $                                                          
            &                                           1         �   � $                      �                        #INITPOINT1D    #         @     @                                                #THIS    #ID    #X              
                                     P               #POINTDT 
             
                                                      
                                      
      1         �   � $                      �                        #INITPOINT2D    #         @     @                                                #THIS    #ID    #X    #Y              
                                     P               #POINTDT 
             
                                                      
                                      
                
                                      
      1         �   � $                      �                        #INITPOINT3D    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                     P               #POINTDT 
             
                                                      
                                      
                
                                      
                
                                      
      1         �   � $                      �                        #SETID     #         @     @                                                 #THIS !   #ID "             
                                !     P               #POINTDT 
             
                                 "           1         �   � $                      �      #                  #SETX $   #         @     @                            $                    #THIS %   #X &             
                                %     P               #POINTDT 
             
                                 &     
      1         �   � $                      �      '                  #SETY (   #         @     @                            (                    #THIS )   #Y *             
                                )     P               #POINTDT 
             
                                 *     
      1         �   � $                      �      +             	     #SETZ ,   #         @     @                            ,                    #THIS -   #Z .             
                                -     P               #POINTDT 
             
                                 .     
      1         �   � $                     �      /             
     #GETID 0   %         @   @                           0                           #THIS 1             
                                1     P               #POINTDT 
   1         �   � $                     �      2              	    #GETX 3   %         @   @                           3                    
       #THIS 4             
                                4     P               #POINTDT 
   1         �   � $                     �      5              
    #GETY 6   %         @   @                           6                    
       #THIS 7             
                                7     P               #POINTDT 
   1         �   � $                     �      8                  #GETZ 9   %         @   @                           9                    
       #THIS :             
                                :     P               #POINTDT 
   1         �   � $                     �      ;                  #GETDIMENSION <   %         @   @                           <                           #THIS =             
                                =     P               #POINTDT 
             � $                              >            P                    #DOFDT ?             &                                                          @  @              D           ?     '                    #VAL @   #FIXEDVAL A   #ISFIXED B   #INIT C   #FIXDOF H   #FREEDOF L                �$                             @                
               � $                             A               
                � $                              B                  1         �   � $                      �      C                  #INIT D   #         @     @                            D                    #THIS E   #DOF F   #ISFIXED G             
                                E                    #DOFDT ?             
                                 F     
                
                                  G           1         �   � $                      �      H                  #FIXDOF I   #         @     @                            I                    #THIS J   #FIXEDVAL K             
                                J                    #DOFDT ?             
                                 K     
      1         �   � $                      �      L                  #FREEDOF M   #         @     @                            M                    #THIS N             
                                N                    #DOFDT ?                �$                              O     P       �              #SOURCEDT P                  @  @              �           P     'P                    #NDIM Q   #FUNC R   #INIT |                � $                             Q                              � $                              R                   8	            #EQUATIONPARSER S             &                                                          @  @              E         S     '8	                   #BYTECODE T   #BYTECODESIZE U   #IMMED V   #IMMEDSIZE W   #STACK X   #STACKSIZE Y   #STACKPTR Z   #FUNCSTRING [   #FUNCSTRINGORIG \   #VARIABLENAMES ]   #EVALUATE ^   #PARSE b   #COMPILE e   #ADDCOMPILEDBYTE h   #COMPILESUBSTR l   #MATHITEMINDEX q   #CHECKSYNTAX v   #FINALIZE y              �$                             T                                         &                                                                   l              y                                                           � $                              U     H                                    l                                                      0               �$                             V            P                
            &                                                                   l              y
                                                           � $                              W     �                                    l                                                      0               �$                             X            �                
            &                                                                   l              y
                                                           � $                              Y     �                                    l                                                      0                � $                              Z     �                                    l                                                      0                � $                             [            �                                     l                                      }              C                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                � $                             \            �      	                              l                                      }              C                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    .           � $                             ]            �             
               &                                                   1         �   � $                     �      ^                  #EVALUATE _   %         @    @                           _                    
       #THIS `   #VAL a                                             `     8	              #EQUATIONPARSER S             
                                 a                   
              &                                           1         �   � D                     �      b                  #PARSE c   #         @     @                            c                    #THIS d                                             d     8	              #EQUATIONPARSER S   1         �   � D                     �      e                  #COMPILE f   #         @     @                            f                    #THIS g                                             g     8	              #EQUATIONPARSER S   1         �   � D                     �      h                  #ADDCOMPILEDBYTE i   #         @     @                            i                    #THIS j   #B k                                             j     8	              #EQUATIONPARSER S             
                                 k           1         �   � D                     �      l                  #COMPILESUBSTR m   #         @ �   @                            m                    #THIS n   #B o   #E p                                             n     8	              #EQUATIONPARSER S             
                                  o                     
                                  p           1         �   � D                    �      q                  #MATHITEMINDEX r   %         @    @                           r                           #THIS s   #B t   #E u                                             s     8	              #EQUATIONPARSER S             
                                  t                     
                                  u           1         �   � D                     �      v                  #CHECKSYNTAX w   #         @     @                            w                    #THIS x                                             x     8	              #EQUATIONPARSER S   2         �   �                              y             #FINALIZE z   #         @     @                           z                    #THIS {                                              {     8	              #EQUATIONPARSER S   1         �   � $                      �      |                  #INIT }   #         @     @                            }                    #THIS ~   #NVAR    #NDIM �   #VAR �   #FUNC �             
                                ~     P               #SOURCEDT P             
                                                      
                                 �           ,         
                                �                         p          5 O p            5 O p                          1 ,         
                                �                         p          5 O p            5 O p                          1 1         �   � $                      �      �                  #INITNODE1D �   #         @     @                            �                    #THIS �   #ID �   #NDOF �   #X �             
                                �     �               #NODEDT              
                                 �                     
                                 �                     
                                 �     
      1         �   � $                      �      �                  #INITNODE2D �   #         @     @                            �                    #THIS �   #ID �   #NDOF �   #X �   #Y �             
                                �     �               #NODEDT              
                                 �                     
                                 �                     
                                 �     
                
                                 �     
      1         �   � $                      �      �                  #INITNODE3D �   #         @     @                            �                    #THIS �   #ID �   #NDOF �   #X �   #Y �   #Z �             
                                �     �               #NODEDT              
                                 �                     
                                 �                     
                                 �     
                
                                 �     
                
                                 �     
      1         �   � $                      �      �                  #ASSIGNSOURCE �   #         @     @                            �                    #THIS �   #SOURCE �             
                                �     �               #NODEDT              
                                 �     P              #SOURCEDT P   1         �   � $                      �      �                  #ASSIGNDOF �   #         @     @                            �                    #THIS �   #IDOF �   #DOF �             
                                �     �               #NODEDT              
                                 �                     
                                 �     
      1         �   � $                      �      �             	     #FIXDOF �   #         @     @                            �                    #THIS �   #IDOF �   #FIXEDVAL �             
                                �     �               #NODEDT              
                                 �                     
                                 �     
      1         �   � $                      �      �             
     #FREEDOF �   #         @     @                            �                    #THIS �   #IDOF �             
                                �     �               #NODEDT              
                                 �           1         �   � $                     �      �                  #GETNDOF �   %         @   @                           �                           #THIS �             
                                �     �               #NODEDT                  @  @               �          �     '                    #NNODE �   #INTEGRATOR �                � $                             �                                � $                              �                         #INTEGRATORDT �                  @  @              D           �     '                   #GAUSSORDER �   #INTEGTERMS �   #WEIGHT �   #GPOINT �   #SHAPEFUNC �   #DSHAPEFUNC �   #DDSHAPEFUNC �   #INIT �   #VALUEGPOINTS �   #GETG1D �   #GETGTRIANGLE �   #GETGSQUARE �   #GETGTETRAHEDRON �                � $                             �                                � $                             �                            � $                             �                             
            &                                                     � $                             �            P                 
            &                   &                                                      � $                             �            �                 
            &                   &                                                      � $                             �                            
            &                   &                   &                                                      � $                             �            �                
            &                   &                   &                   &                                           1         �   � $                      �      �                  #INIT �   #         @     @                            �                    #THIS �   #GAUSSORDER �   #TYPE �             
                                �                   #INTEGRATORDT �             
                                 �                     
                                �                    1 1         �   � $                      �      �             	     #VALUEGPOINTS �   #         @     @                            �                    #THIS �   #TYPE �             
                                �                   #INTEGRATORDT �             
                                �                    1 1         �   � D                     �      �             
     #GETG1D �   #         @     @                            �                    #THIS �             
                                �                   #INTEGRATORDT �   1         �   � D                     �      �                  #GETGTRIANGLE �   #         @     @                            �                    #THIS �             
                                �                   #INTEGRATORDT �   1         �   � D                     �      �                  #GETGSQUARE �   #         @     @                            �                    #THIS �             
                                �                   #INTEGRATORDT �   1         �   � D                     �      �                  #GETGTETRAHEDRON �   #         @     @                            �                    #THIS �             
                                �                   #INTEGRATORDT �                     @               A          �     '�                    #ID �   #NODE �   #GEOMETRY �   #ASSIGNGEOMETRY �   #ASSIGNNODE �   #CALCULATELHS �   #CALCULATERHS �   #CALCULATERESULTS �                � $                             �                              � $                              �                   �             #NODEPTRDT              &                                                        �$                             �            P              #GEOMETRYDT �   1         �   � $                      �      �                  #ASSIGNGEOMETRY �   #         @     @                             �                    #THIS �   #GEOMETRY �             
D                                �     �               #CONDITIONDT �             
                                 �                   #GEOMETRYDT �   1         �   � $                      �      �                  #ASSIGNNODE �   #         @     @                             �                    #THIS �   #INDEX �   #NODE �             
D                                �     �               #CONDITIONDT �             
                                 �                     
                                  �     �              #NODEDT    1         �   � $                     �      �                  #CALCULATELHS �   (        D   @                            �                                   
    #THIS �             &                   &                                                     
                                �     �               #CONDITIONDT �   1         �   � $                     �      �                  #CALCULATERHS �   (        D   @                            �                                   
    #THIS �             &                                                     
                                �     �               #CONDITIONDT �   1         �   � $                      �      �                  #CALCULATERESULTS �   #         @     @                             �                    #THIS �             
                                �     �               #CONDITIONDT �      �   /      fn#fn     �      b   uapp(CONDITIONM    �   @   J  UTILITIESM    +  @   J  DEBUGGERM    k  @   J  NODEM    �  @   J  NODEPTRM    �  @   J  GEOMETRYM #   +  Y      NODEPTRDT+NODEPTRM '   �  \   a   NODEPTRDT%PTR+NODEPTRM    �  �       NODEDT+NODEM %   �  ]   a   NODEDT%POINTDT+NODEM    &  �      POINTDT+POINTM "      H   a   POINTDT%ID+POINTM %   h  �   a   POINTDT%COORD+POINTM +   �  Y   a   POINTDT%INITPOINT1D+POINTM #   U  a      INITPOINT1D+POINTM (   �  U   a   INITPOINT1D%THIS+POINTM &     @   a   INITPOINT1D%ID+POINTM %   K  @   a   INITPOINT1D%X+POINTM +   �  Y   a   POINTDT%INITPOINT2D+POINTM #   �  h      INITPOINT2D+POINTM (   L  U   a   INITPOINT2D%THIS+POINTM &   �  @   a   INITPOINT2D%ID+POINTM %   �  @   a   INITPOINT2D%X+POINTM %   !	  @   a   INITPOINT2D%Y+POINTM +   a	  Y   a   POINTDT%INITPOINT3D+POINTM #   �	  o      INITPOINT3D+POINTM (   )
  U   a   INITPOINT3D%THIS+POINTM &   ~
  @   a   INITPOINT3D%ID+POINTM %   �
  @   a   INITPOINT3D%X+POINTM %   �
  @   a   INITPOINT3D%Y+POINTM %   >  @   a   INITPOINT3D%Z+POINTM %   ~  S   a   POINTDT%SETID+POINTM    �  Z      SETID+POINTM "   +  U   a   SETID%THIS+POINTM     �  @   a   SETID%ID+POINTM $   �  R   a   POINTDT%SETX+POINTM      Y      SETX+POINTM !   k  U   a   SETX%THIS+POINTM    �  @   a   SETX%X+POINTM $      R   a   POINTDT%SETY+POINTM    R  Y      SETY+POINTM !   �  U   a   SETY%THIS+POINTM       @   a   SETY%Y+POINTM $   @  R   a   POINTDT%SETZ+POINTM    �  Y      SETZ+POINTM !   �  U   a   SETZ%THIS+POINTM    @  @   a   SETZ%Z+POINTM %   �  S   a   POINTDT%GETID+POINTM    �  Z      GETID+POINTM "   -  U   a   GETID%THIS+POINTM $   �  R   a   POINTDT%GETX+POINTM    �  Z      GETX+POINTM !   .  U   a   GETX%THIS+POINTM $   �  R   a   POINTDT%GETY+POINTM    �  Z      GETY+POINTM !   /  U   a   GETY%THIS+POINTM $   �  R   a   POINTDT%GETZ+POINTM    �  Z      GETZ+POINTM !   0  U   a   GETZ%THIS+POINTM ,   �  Z   a   POINTDT%GETDIMENSION+POINTM $   �  Z      GETDIMENSION+POINTM )   9  U   a   GETDIMENSION%THIS+POINTM !   �  �   a   NODEDT%DOF+NODEM    -  �      DOFDT+DOFM    �  H   a   DOFDT%VAL+DOFM $     H   a   DOFDT%FIXEDVAL+DOFM #   T  H   a   DOFDT%ISFIXED+DOFM     �  R   a   DOFDT%INIT+DOFM    �  h      INIT+DOFM    V  S   a   INIT%THIS+DOFM    �  @   a   INIT%DOF+DOFM "   �  @   a   INIT%ISFIXED+DOFM "   )  T   a   DOFDT%FIXDOF+DOFM    }  `      FIXDOF+DOFM !   �  S   a   FIXDOF%THIS+DOFM %   0  @   a   FIXDOF%FIXEDVAL+DOFM #   p  U   a   DOFDT%FREEDOF+DOFM    �  R      FREEDOF+DOFM "     S   a   FREEDOF%THIS+DOFM $   j  ^   a   NODEDT%SOURCE+NODEM !   �  n      SOURCEDT+SOURCEM &   6  H   a   SOURCEDT%NDIM+SOURCEM &   ~  �   a   SOURCEDT%FUNC+SOURCEM -   &  i     EQUATIONPARSER+FORTRANPARSER 6   �  �   a   EQUATIONPARSER%BYTECODE+FORTRANPARSER :   �  �   a   EQUATIONPARSER%BYTECODESIZE+FORTRANPARSER 3   (   �   a   EQUATIONPARSER%IMMED+FORTRANPARSER 7   !  �   a   EQUATIONPARSER%IMMEDSIZE+FORTRANPARSER 3   �!  �   a   EQUATIONPARSER%STACK+FORTRANPARSER 7   �"  �   a   EQUATIONPARSER%STACKSIZE+FORTRANPARSER 6   Z#  �   a   EQUATIONPARSER%STACKPTR+FORTRANPARSER 8   �#  �  a   EQUATIONPARSER%FUNCSTRING+FORTRANPARSER <   �(  �  a   EQUATIONPARSER%FUNCSTRINGORIG+FORTRANPARSER ;   y-  �   a   EQUATIONPARSER%VARIABLENAMES+FORTRANPARSER 6   .  V   a   EQUATIONPARSER%EVALUATE+FORTRANPARSER '   k.  c      EVALUATE+FORTRANPARSER ,   �.  \   a   EVALUATE%THIS+FORTRANPARSER +   */  �   a   EVALUATE%VAL+FORTRANPARSER 9   �/  S   %   EQUATIONPARSER%PARSE+FORTRANPARSER=PARSE $   	0  R      PARSE+FORTRANPARSER )   [0  \   a   PARSE%THIS+FORTRANPARSER =   �0  U   %   EQUATIONPARSER%COMPILE+FORTRANPARSER=COMPILE &   1  R      COMPILE+FORTRANPARSER +   ^1  \   a   COMPILE%THIS+FORTRANPARSER M   �1  ]   %   EQUATIONPARSER%ADDCOMPILEDBYTE+FORTRANPARSER=ADDCOMPILEDBYTE .   2  Y      ADDCOMPILEDBYTE+FORTRANPARSER 3   p2  \   a   ADDCOMPILEDBYTE%THIS+FORTRANPARSER 0   �2  @   a   ADDCOMPILEDBYTE%B+FORTRANPARSER I   3  [   %   EQUATIONPARSER%COMPILESUBSTR+FORTRANPARSER=COMPILESUBSTR ,   g3  `      COMPILESUBSTR+FORTRANPARSER 1   �3  \   a   COMPILESUBSTR%THIS+FORTRANPARSER .   #4  @   a   COMPILESUBSTR%B+FORTRANPARSER .   c4  @   a   COMPILESUBSTR%E+FORTRANPARSER I   �4  [   %   EQUATIONPARSER%MATHITEMINDEX+FORTRANPARSER=MATHITEMINDEX ,   �4  h      MATHITEMINDEX+FORTRANPARSER 1   f5  \   a   MATHITEMINDEX%THIS+FORTRANPARSER .   �5  @   a   MATHITEMINDEX%B+FORTRANPARSER .   6  @   a   MATHITEMINDEX%E+FORTRANPARSER E   B6  Y   %   EQUATIONPARSER%CHECKSYNTAX+FORTRANPARSER=CHECKSYNTAX *   �6  R      CHECKSYNTAX+FORTRANPARSER /   �6  \   a   CHECKSYNTAX%THIS+FORTRANPARSER 6   I7  N   a   EQUATIONPARSER%FINALIZE+FORTRANPARSER '   �7  R      FINALIZE+FORTRANPARSER ,   �7  \   a   FINALIZE%THIS+FORTRANPARSER &   E8  R   a   SOURCEDT%INIT+SOURCEM    �8  y      INIT+SOURCEM "   9  V   a   INIT%THIS+SOURCEM "   f9  @   a   INIT%NVAR+SOURCEM "   �9  @   a   INIT%NDIM+SOURCEM !   �9  �   a   INIT%VAR+SOURCEM "   �:  �   a   INIT%FUNC+SOURCEM (   6;  X   a   NODEDT%INITNODE1D+NODEM !   �;  k      INITNODE1D+NODEM &   �;  T   a   INITNODE1D%THIS+NODEM $   M<  @   a   INITNODE1D%ID+NODEM &   �<  @   a   INITNODE1D%NDOF+NODEM #   �<  @   a   INITNODE1D%X+NODEM (   =  X   a   NODEDT%INITNODE2D+NODEM !   e=  r      INITNODE2D+NODEM &   �=  T   a   INITNODE2D%THIS+NODEM $   +>  @   a   INITNODE2D%ID+NODEM &   k>  @   a   INITNODE2D%NDOF+NODEM #   �>  @   a   INITNODE2D%X+NODEM #   �>  @   a   INITNODE2D%Y+NODEM (   +?  X   a   NODEDT%INITNODE3D+NODEM !   �?  y      INITNODE3D+NODEM &   �?  T   a   INITNODE3D%THIS+NODEM $   P@  @   a   INITNODE3D%ID+NODEM &   �@  @   a   INITNODE3D%NDOF+NODEM #   �@  @   a   INITNODE3D%X+NODEM #   A  @   a   INITNODE3D%Y+NODEM #   PA  @   a   INITNODE3D%Z+NODEM *   �A  Z   a   NODEDT%ASSIGNSOURCE+NODEM #   �A  ^      ASSIGNSOURCE+NODEM (   HB  T   a   ASSIGNSOURCE%THIS+NODEM *   �B  V   a   ASSIGNSOURCE%SOURCE+NODEM '   �B  W   a   NODEDT%ASSIGNDOF+NODEM     IC  e      ASSIGNDOF+NODEM %   �C  T   a   ASSIGNDOF%THIS+NODEM %   D  @   a   ASSIGNDOF%IDOF+NODEM $   BD  @   a   ASSIGNDOF%DOF+NODEM $   �D  T   a   NODEDT%FIXDOF+NODEM    �D  j      FIXDOF+NODEM "   @E  T   a   FIXDOF%THIS+NODEM "   �E  @   a   FIXDOF%IDOF+NODEM &   �E  @   a   FIXDOF%FIXEDVAL+NODEM %   F  U   a   NODEDT%FREEDOF+NODEM    iF  \      FREEDOF+NODEM #   �F  T   a   FREEDOF%THIS+NODEM #   G  @   a   FREEDOF%IDOF+NODEM %   YG  U   a   NODEDT%GETNDOF+NODEM    �G  Z      GETNDOF+NODEM #   H  T   a   GETNDOF%THIS+NODEM %   \H  k      GEOMETRYDT+GEOMETRYM +   �H  H   a   GEOMETRYDT%NNODE+GEOMETRYM 0   I  b   a   GEOMETRYDT%INTEGRATOR+GEOMETRYM )   qI       INTEGRATORDT+INTEGRATORM 4   �J  H   a   INTEGRATORDT%GAUSSORDER+INTEGRATORM 4   �J  H   a   INTEGRATORDT%INTEGTERMS+INTEGRATORM 0   K  �   a   INTEGRATORDT%WEIGHT+INTEGRATORM 0   �K  �   a   INTEGRATORDT%GPOINT+INTEGRATORM 3   XL  �   a   INTEGRATORDT%SHAPEFUNC+INTEGRATORM 4   M  �   a   INTEGRATORDT%DSHAPEFUNC+INTEGRATORM 5   �M  �   a   INTEGRATORDT%DDSHAPEFUNC+INTEGRATORM .   �N  R   a   INTEGRATORDT%INIT+INTEGRATORM !   �N  l      INIT+INTEGRATORM &   bO  Z   a   INIT%THIS+INTEGRATORM ,   �O  @   a   INIT%GAUSSORDER+INTEGRATORM &   �O  L   a   INIT%TYPE+INTEGRATORM 6   HP  Z   a   INTEGRATORDT%VALUEGPOINTS+INTEGRATORM )   �P  \      VALUEGPOINTS+INTEGRATORM .   �P  Z   a   VALUEGPOINTS%THIS+INTEGRATORM .   XQ  L   a   VALUEGPOINTS%TYPE+INTEGRATORM 7   �Q  T   %   INTEGRATORDT%GETG1D+INTEGRATORM=GETG1D #   �Q  R      GETG1D+INTEGRATORM (   JR  Z   a   GETG1D%THIS+INTEGRATORM C   �R  Z   %   INTEGRATORDT%GETGTRIANGLE+INTEGRATORM=GETGTRIANGLE )   �R  R      GETGTRIANGLE+INTEGRATORM .   PS  Z   a   GETGTRIANGLE%THIS+INTEGRATORM ?   �S  X   %   INTEGRATORDT%GETGSQUARE+INTEGRATORM=GETGSQUARE '   T  R      GETGSQUARE+INTEGRATORM ,   TT  Z   a   GETGSQUARE%THIS+INTEGRATORM I   �T  ]   %   INTEGRATORDT%GETGTETRAHEDRON+INTEGRATORM=GETGTETRAHEDRON ,   U  R      GETGTETRAHEDRON+INTEGRATORM 1   ]U  Z   a   GETGTETRAHEDRON%THIS+INTEGRATORM    �U  �       CONDITIONDT    �V  H   a   CONDITIONDT%ID !   �V  �   a   CONDITIONDT%NODE %   pW  `   a   CONDITIONDT%GEOMETRY +   �W  \   a   CONDITIONDT%ASSIGNGEOMETRY    ,X  `      ASSIGNGEOMETRY $   �X  Y   a   ASSIGNGEOMETRY%THIS (   �X  X   a   ASSIGNGEOMETRY%GEOMETRY '   =Y  X   a   CONDITIONDT%ASSIGNNODE    �Y  g      ASSIGNNODE     �Y  Y   a   ASSIGNNODE%THIS !   UZ  @   a   ASSIGNNODE%INDEX     �Z  T   a   ASSIGNNODE%NODE )   �Z  Z   a   CONDITIONDT%CALCULATELHS    C[  �      CALCULATELHS "   \  Y   a   CALCULATELHS%THIS )   Z\  Z   a   CONDITIONDT%CALCULATERHS    �\  �      CALCULATERHS "   Z]  Y   a   CALCULATERHS%THIS -   �]  ^   a   CONDITIONDT%CALCULATERESULTS !   ^  R      CALCULATERESULTS &   c^  Y   a   CALCULATERESULTS%THIS 