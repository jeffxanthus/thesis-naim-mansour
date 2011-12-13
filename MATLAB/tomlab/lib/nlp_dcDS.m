% function dc = nlp_dcDS(x, Prob)
%
% nlp_dcDS calls the TOMLAB gateway function nlp_dc
% that computes the constraint gradient at x
%
% nlp_dcDS is used when calling MINLPBB and FILTERSQP in sparse version
%
% nlp_cdcS returns the DYNAMIC SPARSE (DS) non-zero values of dc

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Dec 29, 2002.    Last modified Dec 30, 2002.

function dc = nlp_dcDS(x, Prob)
dc = sparse(feval('nlp_dc',x, Prob)');

% MODIFICATION LOG:
%
% 021229 hkh Written function
% 021230 hkh Changed to only compute dc as sparse