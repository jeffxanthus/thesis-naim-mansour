% TomSym version of qpQG.

Name = 'QP Example';
toms x y
guess = struct('x',0,'y',1);

solution = ezsolve(3*x-4*y+4*(x^2+y^2)+x*y, ...
    {x+y<=5, x==y, x>=0, y>=0}, guess, Name);