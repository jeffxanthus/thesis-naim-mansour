% mcpQG_g - gradient vector for mpec quick guide example
%
% function g = mcpQG_g(x, Prob)

function g = mcpQG_g(x,Prob)

g = [ 2*(x(1)-5), 4*(2*x(2)+1), 0, 0, 0 ]';