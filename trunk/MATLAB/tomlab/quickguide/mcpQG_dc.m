% mcpQG_dc - nonlinear constraint gradient matrix 
%          for mpec quick guide example
%
% function dc = mcpQG_dc(x, Prob)

function dc = mcpQG_dc(x,Prob)

dc = zeros(4,5);

dc(1,:) = [ -1.5  2.0  1.0  -0.5  1.0 ];
dc(2,:) = [  3.0 -1.0  0.0   0.0  0.0 ];
dc(3,:) = [ -1.0  0.5  0.0   0.0  0.0 ];
dc(4,:) = [ -1.0 -1.0  0.0   0.0  0.0 ];