%% Simple Trigonometric Example
% TomSym implementation of GAMS Example (TRIG,SEQ=261)
%
% Simple trigonometric problem from the LGO library
%
% Janos Pinter, LGO - Users Guide, Pinter Consulting Services, Halifax,
% Canada, 2003.

toms x1

obj = sin(11*x1) + cos(13*x1) - sin(17*x1) - cos(19*x1);

eq1 = {-2 <= x1 <= 5
    -x1+5*sin(x1) <= 0};
x0 ={x1 == 1};

options = struct;
options.solver = 'conopt';
solution = ezsolve(obj,eq1,x0,options);