% cpt_dc.m
%
% function [dc,nn]=cpt_dc(x, Prob)
%
% cpt_dc is called by CONOPT to evaluate the Jacobian of the
% nonlinear constraints and the objective function at the point x.
%
% Automatically handles the extended constraint Jacobian when
% both upper and lower bounded constraints are present.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written June 2, 2003.    Last modified Jun 13, 2007.

function [dc,nn] = cpt_dc(x,Prob)

cexidx = Prob.cexidx;
m2     = Prob.m2;

g  = feval('nlp_g',x,Prob);
if m2>0
    dc = feval('nlp_dc',x,Prob);
    if isempty(dc)
        dc = sparse([],[],[],Prob.mNonLin,Prob.N);
    end
    dc = sparse( [ dc ; dc(cexidx(m2+1:end),:) ; g' ] )';
else
    dc = sparse(g);
end

% Number of nonzeros is also returned
nn = nnz(dc);

% MODIFICATION LOG
%
% 030602 ango Wrote file
% 050902 ango Avoid nlp_dc call if m2 is zero
% 070613 ango Handle empty dc from nlp_dc