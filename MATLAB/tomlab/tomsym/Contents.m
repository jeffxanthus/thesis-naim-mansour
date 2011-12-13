% TomSym - A symbolic library for optimization with TOMLAB
% Version 7.7 2011-05-09
%
% TomSym provides symbolic scalars and matrices that can be used to define
% optimization problems within the TOMLAB environment. Unlike normal Matlab
% matrices, which contain numbers, tomSym variables contain symbolic
% expressions. For example:
%
% a = tom('a') % Creates a symbolic variable 'a'
% 
% a = tomSym(1x1):
% 
%   a
% 
% >> b = 0.2*a + cos(a) % Creates an expression involving 'a'
% 
% b = tomSym(1x1):
% 
%   0.2*a+cos(a)
%
% The variables 'a' and 'b' now contain symbolic expressions. These can be
% used to define and solve an optimization problem:
%
% >> ezsolve(b,-5<=a<=5,a==0) % Find a to minimize b, subject to -5 <= a <= 5, 
%                             % with starting guess a=0.
%
% It is also possible to compute derivatives, and to generate matlab code
% from the expressions:
%
% >> mcode(derivative(b,a))
