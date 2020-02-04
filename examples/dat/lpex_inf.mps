NAME          LP1                                                             

*   Problem:
*   min     x1 + x2
*   s.t.    x1 + x2 =  1
*           x1 - x2 =  0
*                x2 =  1
*           x1,  x2 >= 0

ROWS     
 E  ROW1
 E  ROW2
 E  ROW3
 N  COST
COLUMNS
    X1        ROW1                1.   ROW2                1.
    X1        ROW3                0.   COST                1.
    X2        ROW1                1.   ROW2               -1.
    X2        ROW3                1.   COST                1.
RHS
    B         ROW1                1.   ROW3                1.
ENDATA