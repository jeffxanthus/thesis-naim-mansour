% TOMLAB gateway routine 
% Callback from KNITRO
%
% ktr_grad computes the objective gradient g and nonlinear constraint
% Jacobian dc
%
% function [g,dc]=ktr_grad(x, Prob)

% Anders Göran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2007-2007 by Tomlab Optimization Inc., $Release: 5.9.0$
% Written Aug 29, 2007.   Last modified Aug 30, 2007.

function [g,dc] = ktr_grad(x,Prob)

g  = nlp_g(x,Prob);
dc = sparse(nlp_dc(x,Prob));

% MODIFICATION LOG:
% 070829 ango Wrote file
% 070830 ango Add help