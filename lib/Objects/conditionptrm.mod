  ^  �   k820309              19.1        ���^                                                                                                          
       src/condition/ConditionPtr.f90 CONDITIONPTRM              CONDITIONPTRDT                      @                              
                      @  @               A               '�                    #ID    #NODE    #GEOMETRY �   #ASSIGNGEOMETRY �   #ASSIGNNODE �   #CALCULATELHS �   #CALCULATERHS �   #CALCULATERESULTS �                � $                                                           � $                                                 �             #NODEPTRDT              &                                                          @  @                              '�                    #PTR                 �$                                  �                      #NODEDT                   @  @               �               '�                    #POINTDT    #DOF =   #SOURCE N   #INITNODE1D �   #INITNODE2D �   #INITNODE3D �   #ASSIGNSOURCE �   #ASSIGNDOF �   #FIXDOF �   #FREEDOF �   #GETNDOF �                � $                                   P                      #POINTDT 	                 @  @               D           	     'P                    #ID 
   #COORD    #INITPOINT1D    #INITPOINT2D    #INITPOINT3D    #SETID    #SETX "   #SETY &   #SETZ *   #GETID .   #GETX 1   #GETY 4   #GETZ 7   #GETDIMENSION :                � $                             
                             � $                                                          
            &                                           1         �   � $                      �                        #INITPOINT1D    #         @     @                                                #THIS    #ID    #X              
                                     P               #POINTDT 	             
                                                      
                                      
      1         �   � $                      �                        #INITPOINT2D    #         @     @                                                #THIS    #ID    #X    #Y              
                                     P               #POINTDT 	             
                                                      
                                      
                
                                      
      1         �   � $                      �                        #INITPOINT3D    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                     P               #POINTDT 	             
                                                      
                                      
                
                                      
                
                                      
      1         �   � $                      �                        #SETID    #         @     @                                                #THIS     #ID !             
                                      P               #POINTDT 	             
                                 !           1         �   � $                      �      "                  #SETX #   #         @     @                            #                    #THIS $   #X %             
                                $     P               #POINTDT 	             
                                 %     
      1         �   � $                      �      &                  #SETY '   #         @     @                            '                    #THIS (   #Y )             
                                (     P               #POINTDT 	             
                                 )     
      1         �   � $                      �      *             	     #SETZ +   #         @     @                            +                    #THIS ,   #Z -             
                                ,     P               #POINTDT 	             
                                 -     
      1         �   � $                     �      .             
     #GETID /   %         @   @                           /                           #THIS 0             
                                0     P               #POINTDT 	   1         �   � $                     �      1              	    #GETX 2   %         @   @                           2                    
       #THIS 3             
                                3     P               #POINTDT 	   1         �   � $                     �      4              
    #GETY 5   %         @   @                           5                    
       #THIS 6             
                                6     P               #POINTDT 	   1         �   � $                     �      7                  #GETZ 8   %         @   @                           8                    
       #THIS 9             
                                9     P               #POINTDT 	   1         �   � $                     �      :                  #GETDIMENSION ;   %         @   @                           ;                           #THIS <             
                                <     P               #POINTDT 	             � $                              =            P                    #DOFDT >             &                                                          @  @              D           >     '                    #VAL ?   #FIXEDVAL @   #ISFIXED A   #INIT B   #FIXDOF G   #FREEDOF K                �$                             ?                
               � $                             @               
                � $                              A                  1         �   � $                      �      B                  #INIT C   #         @     @                            C                    #THIS D   #DOF E   #ISFIXED F             
                                D                    #DOFDT >             
                                 E     
                
                                  F           1         �   � $                      �      G                  #FIXDOF H   #         @     @                            H                    #THIS I   #FIXEDVAL J             
                                I                    #DOFDT >             
                                 J     
      1         �   � $                      �      K                  #FREEDOF L   #         @     @                            L                    #THIS M             
                                M                    #DOFDT >                �$                              N     P       �              #SOURCEDT O                  @  @              �           O     'P                    #NDIM P   #FUNC Q   #INIT {                � $                             P                              � $                              Q                   8	            #EQUATIONPARSER R             &                                                          @  @              E         R     '8	                   #BYTECODE S   #BYTECODESIZE T   #IMMED U   #IMMEDSIZE V   #STACK W   #STACKSIZE X   #STACKPTR Y   #FUNCSTRING Z   #FUNCSTRINGORIG [   #VARIABLENAMES \   #EVALUATE ]   #PARSE a   #COMPILE d   #ADDCOMPILEDBYTE g   #COMPILESUBSTR k   #MATHITEMINDEX p   #CHECKSYNTAX u   #FINALIZE x              �$                             S                                         &                                                                   �              y                                                           � $                              T     H                                    �                                                      0               �$                             U            P                
            &                                                                   �              y
                                                           � $                              V     �                                    �                                                      0               �$                             W            �                
            &                                                                   �              y
                                                           � $                              X     �                                    �                                                      0                � $                              Y     �                                    �                                                      0                � $                             Z            �                                     �                                      e              C                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                � $                             [            �      	                              �                                      e              C                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    .           � $                             \            �             
               &                                                   1         �   � $                     �      ]                  #EVALUATE ^   %         @    @                           ^                    
       #THIS _   #VAL `                                             _     8	              #EQUATIONPARSER R             
                                 `                   
              &                                           1         �   � D                     �      a                  #PARSE b   #         @     @                            b                    #THIS c                                             c     8	              #EQUATIONPARSER R   1         �   � D                     �      d                  #COMPILE e   #         @     @                            e                    #THIS f                                             f     8	              #EQUATIONPARSER R   1         �   � D                     �      g                  #ADDCOMPILEDBYTE h   #         @     @                            h                    #THIS i   #B j                                             i     8	              #EQUATIONPARSER R             
                                 j           1         �   � D                     �      k                  #COMPILESUBSTR l   #         @ �   @                            l                    #THIS m   #B n   #E o                                             m     8	              #EQUATIONPARSER R             
                                  n                     
                                  o           1         �   � D                    �      p                  #MATHITEMINDEX q   %         @    @                           q                           #THIS r   #B s   #E t                                             r     8	              #EQUATIONPARSER R             
                                  s                     
                                  t           1         �   � D                     �      u                  #CHECKSYNTAX v   #         @     @                            v                    #THIS w                                             w     8	              #EQUATIONPARSER R   2         �   �                              x             #FINALIZE y   #         @     @                           y                    #THIS z                                              z     8	              #EQUATIONPARSER R   1         �   � $                      �      {                  #INIT |   #         @     @                            |                    #THIS }   #NVAR ~   #NDIM    #VAR �   #FUNC �             
                                }     P               #SOURCEDT O             
                                 ~                     
                                            ,         
                                �                         p          5 O p            5 O p                          1 ,         
                                �                         p          5 O p            5 O p                          1 1         �   � $                      �      �                  #INITNODE1D �   #         @     @                            �                    #THIS �   #ID �   #NDOF �   #X �             
                                �     �               #NODEDT              
                                 �                     
                                 �                     
                                 �     
      1         �   � $                      �      �                  #INITNODE2D �   #         @     @                            �                    #THIS �   #ID �   #NDOF �   #X �   #Y �             
                                �     �               #NODEDT              
                                 �                     
                                 �                     
                                 �     
                
                                 �     
      1         �   � $                      �      �                  #INITNODE3D �   #         @     @                            �                    #THIS �   #ID �   #NDOF �   #X �   #Y �   #Z �             
                                �     �               #NODEDT              
                                 �                     
                                 �                     
                                 �     
                
                                 �     
                
                                 �     
      1         �   � $                      �      �                  #ASSIGNSOURCE �   #         @     @                            �                    #THIS �   #SOURCE �             
                                �     �               #NODEDT              
                                 �     P              #SOURCEDT O   1         �   � $                      �      �                  #ASSIGNDOF �   #         @     @                            �                    #THIS �   #IDOF �   #DOF �             
                                �     �               #NODEDT              
                                 �                     
                                 �     
      1         �   � $                      �      �             	     #FIXDOF �   #         @     @                            �                    #THIS �   #IDOF �   #FIXEDVAL �             
                                �     �               #NODEDT              
                                 �                     
                                 �     
      1         �   � $                      �      �             
     #FREEDOF �   #         @     @                            �                    #THIS �   #IDOF �             
                                �     �               #NODEDT              
                                 �           1         �   � $                     �      �                  #GETNDOF �   %         @   @                           �                           #THIS �             
                                �     �               #NODEDT                 �$                             �            P              #GEOMETRYDT �                 @  @               �          �     '                    #NNODE �   #INTEGRATOR �                � $                             �                                � $                              �                         #INTEGRATORDT �                  @  @              D           �     '                   #GAUSSORDER �   #INTEGTERMS �   #WEIGHT �   #GPOINT �   #SHAPEFUNC �   #DSHAPEFUNC �   #DDSHAPEFUNC �   #INIT �   #VALUEGPOINTS �   #GETG1D �   #GETGTRIANGLE �   #GETGSQUARE �   #GETGTETRAHEDRON �                � $                             �                                � $                             �                            � $                             �                             
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
                                �                   #INTEGRATORDT �   1         �   � $                      �      �                  #ASSIGNGEOMETRY �   #         @     @                            �                    #THIS �   #GEOMETRY �             
                                �     �               #CONDITIONDT              
                                 �                   #GEOMETRYDT �   1         �   � $                      �      �                  #ASSIGNNODE �   #         @     @                            �                    #THIS �   #INDEX �   #NODE �             
                                �     �               #CONDITIONDT              
                                 �                     
                                  �     �              #NODEDT    1         �   � $                     �      �                  #CALCULATELHS �   (        D   @                           �                                   
    #THIS �             &                   &                                                     
                                �     �               #CONDITIONDT    1         �   � $                     �      �                  #CALCULATERHS �   (        D   @                           �                                   
    #THIS �             &                                                     
                                �     �               #CONDITIONDT    1         �   � $                      �      �                  #CALCULATERESULTS �   #         @     @                            �                    #THIS �             
                                �     �               #CONDITIONDT                      @                           �     '�                    #PTR �                �$                             �     �                      #CONDITIONDT       �   5      fn#fn #   �      b   uapp(CONDITIONPTRM    �   @   J  CONDITIONM '   4  �      CONDITIONDT+CONDITIONM *     H   a   CONDITIONDT%ID+CONDITIONM ,   J  �   a   CONDITIONDT%NODE+CONDITIONM #   �  Y      NODEPTRDT+NODEPTRM '   F  \   a   NODEPTRDT%PTR+NODEPTRM    �  �      NODEDT+NODEM %   �  ]   a   NODEDT%POINTDT+NODEM    �  �      POINTDT+POINTM "   �  H   a   POINTDT%ID+POINTM %   *  �   a   POINTDT%COORD+POINTM +   �  Y   a   POINTDT%INITPOINT1D+POINTM #     a      INITPOINT1D+POINTM (   x  U   a   INITPOINT1D%THIS+POINTM &   �  @   a   INITPOINT1D%ID+POINTM %     @   a   INITPOINT1D%X+POINTM +   M  Y   a   POINTDT%INITPOINT2D+POINTM #   �  h      INITPOINT2D+POINTM (   	  U   a   INITPOINT2D%THIS+POINTM &   c	  @   a   INITPOINT2D%ID+POINTM %   �	  @   a   INITPOINT2D%X+POINTM %   �	  @   a   INITPOINT2D%Y+POINTM +   #
  Y   a   POINTDT%INITPOINT3D+POINTM #   |
  o      INITPOINT3D+POINTM (   �
  U   a   INITPOINT3D%THIS+POINTM &   @  @   a   INITPOINT3D%ID+POINTM %   �  @   a   INITPOINT3D%X+POINTM %   �  @   a   INITPOINT3D%Y+POINTM %      @   a   INITPOINT3D%Z+POINTM %   @  S   a   POINTDT%SETID+POINTM    �  Z      SETID+POINTM "   �  U   a   SETID%THIS+POINTM     B  @   a   SETID%ID+POINTM $   �  R   a   POINTDT%SETX+POINTM    �  Y      SETX+POINTM !   -  U   a   SETX%THIS+POINTM    �  @   a   SETX%X+POINTM $   �  R   a   POINTDT%SETY+POINTM      Y      SETY+POINTM !   m  U   a   SETY%THIS+POINTM    �  @   a   SETY%Y+POINTM $     R   a   POINTDT%SETZ+POINTM    T  Y      SETZ+POINTM !   �  U   a   SETZ%THIS+POINTM      @   a   SETZ%Z+POINTM %   B  S   a   POINTDT%GETID+POINTM    �  Z      GETID+POINTM "   �  U   a   GETID%THIS+POINTM $   D  R   a   POINTDT%GETX+POINTM    �  Z      GETX+POINTM !   �  U   a   GETX%THIS+POINTM $   E  R   a   POINTDT%GETY+POINTM    �  Z      GETY+POINTM !   �  U   a   GETY%THIS+POINTM $   F  R   a   POINTDT%GETZ+POINTM    �  Z      GETZ+POINTM !   �  U   a   GETZ%THIS+POINTM ,   G  Z   a   POINTDT%GETDIMENSION+POINTM $   �  Z      GETDIMENSION+POINTM )   �  U   a   GETDIMENSION%THIS+POINTM !   P  �   a   NODEDT%DOF+NODEM    �  �      DOFDT+DOFM    �  H   a   DOFDT%VAL+DOFM $   �  H   a   DOFDT%FIXEDVAL+DOFM #     H   a   DOFDT%ISFIXED+DOFM     ^  R   a   DOFDT%INIT+DOFM    �  h      INIT+DOFM      S   a   INIT%THIS+DOFM    k  @   a   INIT%DOF+DOFM "   �  @   a   INIT%ISFIXED+DOFM "   �  T   a   DOFDT%FIXDOF+DOFM    ?  `      FIXDOF+DOFM !   �  S   a   FIXDOF%THIS+DOFM %   �  @   a   FIXDOF%FIXEDVAL+DOFM #   2  U   a   DOFDT%FREEDOF+DOFM    �  R      FREEDOF+DOFM "   �  S   a   FREEDOF%THIS+DOFM $   ,  ^   a   NODEDT%SOURCE+NODEM !   �  n      SOURCEDT+SOURCEM &   �  H   a   SOURCEDT%NDIM+SOURCEM &   @  �   a   SOURCEDT%FUNC+SOURCEM -   �  i     EQUATIONPARSER+FORTRANPARSER 6   Q  �   a   EQUATIONPARSER%BYTECODE+FORTRANPARSER :   E   �   a   EQUATIONPARSER%BYTECODESIZE+FORTRANPARSER 3   �   �   a   EQUATIONPARSER%IMMED+FORTRANPARSER 7   �!  �   a   EQUATIONPARSER%IMMEDSIZE+FORTRANPARSER 3   �"  �   a   EQUATIONPARSER%STACK+FORTRANPARSER 7   w#  �   a   EQUATIONPARSER%STACKSIZE+FORTRANPARSER 6   $  �   a   EQUATIONPARSER%STACKPTR+FORTRANPARSER 8   �$  �  a   EQUATIONPARSER%FUNCSTRING+FORTRANPARSER <   ~)  �  a   EQUATIONPARSER%FUNCSTRINGORIG+FORTRANPARSER ;   ;.  �   a   EQUATIONPARSER%VARIABLENAMES+FORTRANPARSER 6   �.  V   a   EQUATIONPARSER%EVALUATE+FORTRANPARSER '   -/  c      EVALUATE+FORTRANPARSER ,   �/  \   a   EVALUATE%THIS+FORTRANPARSER +   �/  �   a   EVALUATE%VAL+FORTRANPARSER 9   x0  S   %   EQUATIONPARSER%PARSE+FORTRANPARSER=PARSE $   �0  R      PARSE+FORTRANPARSER )   1  \   a   PARSE%THIS+FORTRANPARSER =   y1  U   %   EQUATIONPARSER%COMPILE+FORTRANPARSER=COMPILE &   �1  R      COMPILE+FORTRANPARSER +    2  \   a   COMPILE%THIS+FORTRANPARSER M   |2  ]   %   EQUATIONPARSER%ADDCOMPILEDBYTE+FORTRANPARSER=ADDCOMPILEDBYTE .   �2  Y      ADDCOMPILEDBYTE+FORTRANPARSER 3   23  \   a   ADDCOMPILEDBYTE%THIS+FORTRANPARSER 0   �3  @   a   ADDCOMPILEDBYTE%B+FORTRANPARSER I   �3  [   %   EQUATIONPARSER%COMPILESUBSTR+FORTRANPARSER=COMPILESUBSTR ,   )4  `      COMPILESUBSTR+FORTRANPARSER 1   �4  \   a   COMPILESUBSTR%THIS+FORTRANPARSER .   �4  @   a   COMPILESUBSTR%B+FORTRANPARSER .   %5  @   a   COMPILESUBSTR%E+FORTRANPARSER I   e5  [   %   EQUATIONPARSER%MATHITEMINDEX+FORTRANPARSER=MATHITEMINDEX ,   �5  h      MATHITEMINDEX+FORTRANPARSER 1   (6  \   a   MATHITEMINDEX%THIS+FORTRANPARSER .   �6  @   a   MATHITEMINDEX%B+FORTRANPARSER .   �6  @   a   MATHITEMINDEX%E+FORTRANPARSER E   7  Y   %   EQUATIONPARSER%CHECKSYNTAX+FORTRANPARSER=CHECKSYNTAX *   ]7  R      CHECKSYNTAX+FORTRANPARSER /   �7  \   a   CHECKSYNTAX%THIS+FORTRANPARSER 6   8  N   a   EQUATIONPARSER%FINALIZE+FORTRANPARSER '   Y8  R      FINALIZE+FORTRANPARSER ,   �8  \   a   FINALIZE%THIS+FORTRANPARSER &   9  R   a   SOURCEDT%INIT+SOURCEM    Y9  y      INIT+SOURCEM "   �9  V   a   INIT%THIS+SOURCEM "   (:  @   a   INIT%NVAR+SOURCEM "   h:  @   a   INIT%NDIM+SOURCEM !   �:  �   a   INIT%VAR+SOURCEM "   P;  �   a   INIT%FUNC+SOURCEM (   �;  X   a   NODEDT%INITNODE1D+NODEM !   P<  k      INITNODE1D+NODEM &   �<  T   a   INITNODE1D%THIS+NODEM $   =  @   a   INITNODE1D%ID+NODEM &   O=  @   a   INITNODE1D%NDOF+NODEM #   �=  @   a   INITNODE1D%X+NODEM (   �=  X   a   NODEDT%INITNODE2D+NODEM !   '>  r      INITNODE2D+NODEM &   �>  T   a   INITNODE2D%THIS+NODEM $   �>  @   a   INITNODE2D%ID+NODEM &   -?  @   a   INITNODE2D%NDOF+NODEM #   m?  @   a   INITNODE2D%X+NODEM #   �?  @   a   INITNODE2D%Y+NODEM (   �?  X   a   NODEDT%INITNODE3D+NODEM !   E@  y      INITNODE3D+NODEM &   �@  T   a   INITNODE3D%THIS+NODEM $   A  @   a   INITNODE3D%ID+NODEM &   RA  @   a   INITNODE3D%NDOF+NODEM #   �A  @   a   INITNODE3D%X+NODEM #   �A  @   a   INITNODE3D%Y+NODEM #   B  @   a   INITNODE3D%Z+NODEM *   RB  Z   a   NODEDT%ASSIGNSOURCE+NODEM #   �B  ^      ASSIGNSOURCE+NODEM (   
C  T   a   ASSIGNSOURCE%THIS+NODEM *   ^C  V   a   ASSIGNSOURCE%SOURCE+NODEM '   �C  W   a   NODEDT%ASSIGNDOF+NODEM     D  e      ASSIGNDOF+NODEM %   pD  T   a   ASSIGNDOF%THIS+NODEM %   �D  @   a   ASSIGNDOF%IDOF+NODEM $   E  @   a   ASSIGNDOF%DOF+NODEM $   DE  T   a   NODEDT%FIXDOF+NODEM    �E  j      FIXDOF+NODEM "   F  T   a   FIXDOF%THIS+NODEM "   VF  @   a   FIXDOF%IDOF+NODEM &   �F  @   a   FIXDOF%FIXEDVAL+NODEM %   �F  U   a   NODEDT%FREEDOF+NODEM    +G  \      FREEDOF+NODEM #   �G  T   a   FREEDOF%THIS+NODEM #   �G  @   a   FREEDOF%IDOF+NODEM %   H  U   a   NODEDT%GETNDOF+NODEM    pH  Z      GETNDOF+NODEM #   �H  T   a   GETNDOF%THIS+NODEM 0   I  `   a   CONDITIONDT%GEOMETRY+CONDITIONM %   ~I  k      GEOMETRYDT+GEOMETRYM +   �I  H   a   GEOMETRYDT%NNODE+GEOMETRYM 0   1J  b   a   GEOMETRYDT%INTEGRATOR+GEOMETRYM )   �J       INTEGRATORDT+INTEGRATORM 4   �K  H   a   INTEGRATORDT%GAUSSORDER+INTEGRATORM 4   �K  H   a   INTEGRATORDT%INTEGTERMS+INTEGRATORM 0   :L  �   a   INTEGRATORDT%WEIGHT+INTEGRATORM 0   �L  �   a   INTEGRATORDT%GPOINT+INTEGRATORM 3   zM  �   a   INTEGRATORDT%SHAPEFUNC+INTEGRATORM 4   &N  �   a   INTEGRATORDT%DSHAPEFUNC+INTEGRATORM 5   �N  �   a   INTEGRATORDT%DDSHAPEFUNC+INTEGRATORM .   �O  R   a   INTEGRATORDT%INIT+INTEGRATORM !   P  l      INIT+INTEGRATORM &   �P  Z   a   INIT%THIS+INTEGRATORM ,   �P  @   a   INIT%GAUSSORDER+INTEGRATORM &   Q  L   a   INIT%TYPE+INTEGRATORM 6   jQ  Z   a   INTEGRATORDT%VALUEGPOINTS+INTEGRATORM )   �Q  \      VALUEGPOINTS+INTEGRATORM .    R  Z   a   VALUEGPOINTS%THIS+INTEGRATORM .   zR  L   a   VALUEGPOINTS%TYPE+INTEGRATORM 7   �R  T   %   INTEGRATORDT%GETG1D+INTEGRATORM=GETG1D #   S  R      GETG1D+INTEGRATORM (   lS  Z   a   GETG1D%THIS+INTEGRATORM C   �S  Z   %   INTEGRATORDT%GETGTRIANGLE+INTEGRATORM=GETGTRIANGLE )    T  R      GETGTRIANGLE+INTEGRATORM .   rT  Z   a   GETGTRIANGLE%THIS+INTEGRATORM ?   �T  X   %   INTEGRATORDT%GETGSQUARE+INTEGRATORM=GETGSQUARE '   $U  R      GETGSQUARE+INTEGRATORM ,   vU  Z   a   GETGSQUARE%THIS+INTEGRATORM I   �U  ]   %   INTEGRATORDT%GETGTETRAHEDRON+INTEGRATORM=GETGTETRAHEDRON ,   -V  R      GETGTETRAHEDRON+INTEGRATORM 1   V  Z   a   GETGTETRAHEDRON%THIS+INTEGRATORM 6   �V  \   a   CONDITIONDT%ASSIGNGEOMETRY+CONDITIONM *   5W  `      ASSIGNGEOMETRY+CONDITIONM /   �W  Y   a   ASSIGNGEOMETRY%THIS+CONDITIONM 3   �W  X   a   ASSIGNGEOMETRY%GEOMETRY+CONDITIONM 2   FX  X   a   CONDITIONDT%ASSIGNNODE+CONDITIONM &   �X  g      ASSIGNNODE+CONDITIONM +   Y  Y   a   ASSIGNNODE%THIS+CONDITIONM ,   ^Y  @   a   ASSIGNNODE%INDEX+CONDITIONM +   �Y  T   a   ASSIGNNODE%NODE+CONDITIONM 4   �Y  Z   a   CONDITIONDT%CALCULATELHS+CONDITIONM (   LZ  �      CALCULATELHS+CONDITIONM -   
[  Y   a   CALCULATELHS%THIS+CONDITIONM 4   c[  Z   a   CONDITIONDT%CALCULATERHS+CONDITIONM (   �[  �      CALCULATERHS+CONDITIONM -   c\  Y   a   CALCULATERHS%THIS+CONDITIONM 8   �\  ^   a   CONDITIONDT%CALCULATERESULTS+CONDITIONM ,   ]  R      CALCULATERESULTS+CONDITIONM 1   l]  Y   a   CALCULATERESULTS%THIS+CONDITIONM    �]  Y       CONDITIONPTRDT #   ^  a   a   CONDITIONPTRDT%PTR 