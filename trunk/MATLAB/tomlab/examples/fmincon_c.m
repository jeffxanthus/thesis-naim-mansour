% fmincon_c - test problem for fmincon
%
% Nonlinear constraints from file glc4_c: 
% Nonlinear constraints for Floudas-Pardalos 3.3
%
% function [c,cEQ] = fmincon_c(x, Prob)

function [c,cEQ] = fmincon_c(x, Prob)

% Two nonlinear constraints (quadratic)
% Must be written as c(x) <= 0
%
% Was 4 <= c(x) as glc4_c.m

c   = [4-(x(3)-3)^2+x(4); 4-(x(5)-3)^2+x(6)];

cEQ = [];