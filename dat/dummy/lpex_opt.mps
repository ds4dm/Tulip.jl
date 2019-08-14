NAME          LP1                                                             

*   Problem:
*   min     x1 + 2*x2
*   s.t.    x1 +   x2 = 1
*           x1 -   x2 = 0
*           0 <= x1, x2, <= 1

ROWS     
 E  ROW1
 E  ROW2
 N  COST
COLUMNS
    X1        ROW1                1.   ROW2                1.
    X1        COST                1.
    X2        ROW1                1.   ROW2               -1.
    X2        COST                2.
RHS
    B         ROW1                1.   ROW2                0.
BOUNDS
 UP BND1      X1                  1.
 UP BND1      X2                  1.
ENDATA