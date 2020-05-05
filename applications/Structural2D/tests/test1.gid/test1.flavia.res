GiD Post Result File 1.0
 Result "Displacement" "test1"            1  Vector onNodes
 Values
           1  1.067243939085251E-007  1.440004209070755E-007
           2  6.318098928363887E-008  1.156215759330755E-007
           3  7.834554893452516E-008  8.181900278291410E-008
           4  5.059274315131458E-008  6.713621669354126E-008
           5  5.851567883602891E-008  1.087685865151422E-007
           6  7.149255951659182E-008  5.851567883602896E-008
           7  5.166268941809558E-008  7.916894987308563E-008
           8  6.053093637381042E-008  5.166268941809564E-008
           9  0.000000000000000E+000  0.000000000000000E+000
End Values

GaussPoints "PointsNormalStressOnTriangs" ElemType Triangle
Number of GaussPoints: 4
Natural Coordinates: Given
        0.3333333333333333          0.3333333333333333
        0.2000000000000000          0.2000000000000000
        0.6000000000000000          0.2000000000000000
        0.2000000000000000          0.6000000000000000
End gausspoints
Result "NormalStressOnTriangs" "test1" 1 Vector onGaussPoints "PointsNormalStressOnTriangs"
Values
1      22921.3747866923076799      22921.3747866923222318
          22921.3747866923076799      22921.3747866923222318
          22921.3747866923076799      22921.3747866923222318
          22921.3747866923076799      22921.3747866923222318
2       5535.1068375615204786       5535.1068375615077457
           5535.1068375615204786       5535.1068375615077457
           5535.1068375615204786       5535.1068375615077457
           5535.1068375615204786       5535.1068375615077457
3      10520.9128077964796830      28353.2808248799883586
          10520.9128077964796830      28353.2808248799883586
          10520.9128077964796830      28353.2808248799883586
          10520.9128077964796830      28353.2808248799883586
4       5535.1068375615377590       5535.1068375615332116
           5535.1068375615377590       5535.1068375615332116
           5535.1068375615377590       5535.1068375615332116
           5535.1068375615377590       5535.1068375615332116
5       7177.7048719721369707       7177.7048719721369707
           7177.7048719721369707       7177.7048719721369707
           7177.7048719721369707       7177.7048719721369707
           7177.7048719721369707       7177.7048719721369707
6      -1750.4187071664032374       4462.3526741182504338
          -1750.4187071664032374       4462.3526741182504338
          -1750.4187071664032374       4462.3526741182504338
          -1750.4187071664032374       4462.3526741182504338
7      28353.2808248799665307      10520.9128077964815020
          28353.2808248799665307      10520.9128077964815020
          28353.2808248799665307      10520.9128077964815020
          28353.2808248799665307      10520.9128077964815020
8       4462.3526741182540718      -1750.4187071663973256
           4462.3526741182540718      -1750.4187071663973256
           4462.3526741182540718      -1750.4187071663973256
           4462.3526741182540718      -1750.4187071663973256
End Values

GaussPoints "PointsShearStressOnTriangs" ElemType Triangle
Number of GaussPoints: 4
Natural Coordinates: Given
        0.3333333333333333          0.3333333333333333
        0.2000000000000000          0.2000000000000000
        0.6000000000000000          0.2000000000000000
        0.2000000000000000          0.6000000000000000
End gausspoints
Result "ShearStressOnTriangs" "test1" 1 Scalar onGaussPoints "PointsShearStressOnTriangs"
Values
1      0.1707862521330770E+05
          0.1707862521330770E+05
          0.1707862521330770E+05
          0.1707862521330770E+05
2      0.5535106837561517E+04
          0.5535106837561517E+04
          0.5535106837561517E+04
          0.5535106837561517E+04
3      0.1365072069287991E+05
          0.1365072069287991E+05
          0.1365072069287991E+05
          0.1365072069287991E+05
4      0.5535106837561517E+04
          0.5535106837561517E+04
          0.5535106837561517E+04
          0.5535106837561517E+04
5      0.1078115724741963E+05
          0.1078115724741963E+05
          0.1078115724741963E+05
          0.1078115724741963E+05
6      0.7206762779195761E+04
          0.7206762779195761E+04
          0.7206762779195761E+04
          0.7206762779195761E+04
7      0.1365072069287992E+05
          0.1365072069287992E+05
          0.1365072069287992E+05
          0.1365072069287992E+05
8      0.7206762779195752E+04
          0.7206762779195752E+04
          0.7206762779195752E+04
          0.7206762779195752E+04
End Values

GaussPoints "PointsStrainStressOnTriangs" ElemType Triangle
Number of GaussPoints: 4
Natural Coordinates: Given
        0.3333333333333333          0.3333333333333333
        0.2000000000000000          0.2000000000000000
        0.6000000000000000          0.2000000000000000
        0.2000000000000000          0.6000000000000000
End gausspoints
Result "StrainStressOnTriangs" "test1" 1 Vector onGaussPoints "PointsStrainStressOnTriangs"
Values
1          0.0000001135153799          0.0000000000000000
              0.0000001135153799          0.0000000000000000
              0.0000001135153799          0.0000000000000000
              0.0000001135153799          0.0000000000000000
2          0.0000000274119577          0.0000000000000000
              0.0000000274119577          0.0000000000000000
              0.0000000274119577          0.0000000000000000
              0.0000000274119577          0.0000000000000000
3          0.0000000962599080          0.0000000000000000
              0.0000000962599080          0.0000000000000000
              0.0000000962599080          0.0000000000000000
              0.0000000962599080          0.0000000000000000
4          0.0000000274119577          0.0000000000000000
              0.0000000274119577          0.0000000000000000
              0.0000000274119577          0.0000000000000000
              0.0000000274119577          0.0000000000000000
5          0.0000000355467289          0.0000000000000000
              0.0000000355467289          0.0000000000000000
              0.0000000355467289          0.0000000000000000
              0.0000000355467289          0.0000000000000000
6          0.0000000067152651          0.0000000000000000
              0.0000000067152651          0.0000000000000000
              0.0000000067152651          0.0000000000000000
              0.0000000067152651          0.0000000000000000
7          0.0000000962599080          0.0000000000000000
              0.0000000962599080          0.0000000000000000
              0.0000000962599080          0.0000000000000000
              0.0000000962599080          0.0000000000000000
8          0.0000000067152651          0.0000000000000000
              0.0000000067152651          0.0000000000000000
              0.0000000067152651          0.0000000000000000
              0.0000000067152651          0.0000000000000000
End Values
