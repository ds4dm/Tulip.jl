NAME          LP_FREE_VARS

*   Problem:
*   min     x1   + x2 + x3
*   s.t.    2 x1 + x2 >= 2
*           x1 + 2 x2 >= 2
*           x1 + x2 + x3 >= 0

ROWS     
 G  ROW1
 G  ROW2
 G  ROW3
 N  COST
COLUMNS
    X1        ROW1                2.   ROW2                1.
    X1        COST                1.   ROW3                1.
    X2        ROW1                1.   ROW2                2.
    X2        COST                1.   ROW3                1.
    X3        ROW3                1.   COST                1.
RHS
    B         ROW1                2.   ROW2                2.
BOUNDS
 FR BND1      X1
 FR BND1      X2
 FR BND1      X3
ENDATA