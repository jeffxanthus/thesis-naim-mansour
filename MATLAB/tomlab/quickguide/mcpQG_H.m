% mcpQG_H - Hessian matrix for mpec quick guide example
%
% function H = mcpQG_H(x, Prob)

function H = mcpQG_H(x,Prob)

% H = [ 2 0 0 0 0 ; 0 8 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 0 ];

H = zeros(5,5);
H(1,1) = 2.0;
H(2,2) = 8.0;