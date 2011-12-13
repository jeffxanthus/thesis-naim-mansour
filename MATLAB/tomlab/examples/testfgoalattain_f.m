function f = testfgoalattain_f(x,Prob)

A =  [ -0.5  0  0;  0  -2  10;  0  1  -2 ];
B =  [ 1  0;  -2  2;  0  1 ];
C =  [ 1  0  0;  0  0  1 ];
f = sort(eig(A+B*reshape(x,2,2)*C));