% funchandleQG is a small example for defining and solving problems
% using nested functions, anonymous functions and subfunctions.

function [R1, R2, R3] = funchandleQG

Name='funchandleQG';

x_L = [0 0 0 0 0 ]';
x_U = [inf inf 1 1 1 ]';

A = [1 0     1  0 0 ; ...
    0 1.333 0  1 0 ; ...
    0 0    -1 -1 1 ];

b_L = [];
b_U = [1.6 ; 3 ; 0];

c_L = [1.25;3];
c_U = c_L;

x_0 = ones(5,1);

HessPattern = spalloc(5,5,0);
ConsPattern = [ 1 0 1 0 0; ...
    0 1 0 1 0 ];

% Example with local sub functions
Prob1 = minlpAssign(@my_f, @my_g, @my_H, HessPattern, ...
    x_L, x_U, Name, x_0, [], [], [], [], ...
    A, b_L, b_U, @my_c, @my_dc, @my_d2c, ...
    ConsPattern, c_L, c_U);

R1 = tomRun('knitro', Prob1, 1);

constr = @(x) [ x(1)^2+x(3) ; sqrt(x(2)^3)+1.5*x(4)];

% Example with local sub functions and anonymous function
% One directly into assign routine and one as variable input
Prob2 = minlpAssign(@(x) [2 3 1.5 2 -0.5]*x, @my_g, @my_H, HessPattern, ...
    x_L, x_U, Name, x_0, [], [], [], [], ...
    A, b_L, b_U, constr, @my_dc, @my_d2c, ...
    ConsPattern, c_L, c_U);

R2 = tomRun('knitro', Prob2, 1);

% Example with local sub functions, anonymous function
% and nested function
Prob3 = minlpAssign(@mynested_f, @my_g, @my_H, HessPattern, ...
    x_L, x_U, Name, x_0, [], [], [], [], ...
    A, b_L, b_U, constr, @my_dc, @my_d2c, ...
    ConsPattern, c_L, c_U);

    function f = mynested_f(x)
        f = [2 3 1.5 2 -0.5]*x;
    end

R3 = tomRun('knitro', Prob3, 1);

end

function f = my_f(x)
f = [2 3 1.5 2 -0.5]*x;
end

function g = my_g(x)
g = [2 3 1.5 2 -0.5]';
end

function H = my_H(x)
H = spalloc(5,5,0);
end

function c = my_c(x)
c = [ x(1)^2+x(3) ; sqrt(x(2)^3)+1.5*x(4)];
end

function dc = my_dc(x)
dc = [2*x(1)     0.0       1.0  0.0  0.0 ; ...
    0.0   1.5*sqrt(x(2)) 0.0  1.5  0.0];
end

function d2c = my_d2c(x,lam)
d2c      = spalloc(5,5,2);
d2c(1,1) = lam(1)*2;
d2c(2,2) = lam(2)*3/(4*sqrt(x(2)));
end