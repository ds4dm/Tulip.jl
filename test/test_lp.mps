NAME          TESTLP                                                             
ROWS     
 E  ROW1
 L  ROW2
 G  ROW3
 N  COST
COLUMNS
    X1        ROW1                11   ROW2                21
    X1        ROW3                31   COST                01
    X2        ROW1                12   ROW2                22
    X2        ROW3                32   COST                02
    X3        ROW1                13   ROW2                23
    X3        ROW3                33   COST                03
    X4        ROW1                14   ROW2                24
    X4        ROW3                34   COST                04
    X5        ROW1                15   ROW2                25
    X5        ROW3                35   COST                05
    X6        ROW1                16   ROW2                26
    X6        ROW3                36   COST                06
    X7        ROW1                17   ROW2                27
    X7        ROW3                37   COST                07
    X8        ROW1                18   ROW2                28
    X8        ROW3                38   COST                08
    X9        ROW1                19   ROW2                29
    X9        ROW3                39   COST                09
RHS
    B         ROW1               1.0   ROW2                2.0
    B         ROW3               3.0    
* Following lines should be ignored
    B2        ROW1                1.   ROW2                0.
BOUNDS
 LO BND1      X1                  1.
 UP BND1      X2                  1.
 FX BND1      X3                  1.
 FR BND1      X4
 MI BND1      X5
 UP BND1      X5                  5.
 PL BND1      X6
 BV BND1      X7
 LI BND1      X8                  1.
 UI BND1      X9                  1.
* Following lines should be ignored
 UP BND2      X1                  2.
ENDATA