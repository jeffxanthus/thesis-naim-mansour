% TomSym version of lpQG.

c     = [-7 -5]';
A     = [ 1  2; 4  1 ];
b_U   = [ 6 12 ]';
x_L   = [ 0  0 ]';

toms 2x1 x       % Symbolic variable

solution = ezsolve(c'*x, {A*x<=b_U, x_L<=x});

disp('x =');
disp(solution.x);