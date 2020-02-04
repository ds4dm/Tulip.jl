NAME          LP2                                                          

*   Problem:
*   min     -x1 - x2
*   s.t.     x1 - x2 =  1
*            x1,  x2 >= 0

ROWS     
 E  ROW1
 N  COST
COLUMNS
    X1        ROW1                1.   COST               -1.
    X2        ROW1               -1.   COST               -1.
RHS
    B         ROW1                1.
ENDATA