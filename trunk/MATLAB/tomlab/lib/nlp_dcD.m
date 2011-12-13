% function z = nlp_dcD(x, Prob)
%
% nlp_dcD calls the TOMLAB gateway function nlp_dc, 
% which computes the gradient for all constraints at x, dc, for test problem P
%
% nlp_dcD is used when calling MINLPBB and FILTERSQP in dense version
%
% dc could be either dense or sparse, and converted to a full (dense) matrix
% dc is always first transposed because in the Tomlab format,
% every row corresponds to a constraint whereas in the DUNDEE format
% every row corresponds to a variable

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Dec 23, 2002.   Last modified Dec 23, 2002.

function z = nlp_dcD(x, Prob)

% Call nlp_dc to get the constraint gradient
dc = feval('nlp_dc',x, Prob)';

if isempty(dc)
   z=[];
else
   if issparse(dc)
      z=full(dc);
      z=z(:);
   else
      z=dc(:);
   end
end

% MODIFICATION LOG:
%
% 021223 hkh  Written