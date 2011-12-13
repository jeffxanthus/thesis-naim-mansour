%-----------------------------------------------------
% function f = un_f(y, Prob);
%
% Compute un_f using tomsol dll
%-----------------------------------------------------

function f = un_f(y, Prob)

f = abs(tomsol(26,y));