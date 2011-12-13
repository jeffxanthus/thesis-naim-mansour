% function z = ktr_dcS(x, Prob)
%
% ktr_dcS calls the TOMLAB gateway function nlp_dc,
% that computes the gradient dc for all constraints at x
%
% nlp_dcS is used by KNITRO for the sparse nonlinear Jacobian. 
%
% nlp_dcS returns the static non-zero values of dcS in the vector z
%
% dcS could be either dense or sparse.
% The non-zero pattern is described in Prob.ConsPattern (should be sparse)
% If isempty(Prob.ConsPattern), then a dense matrix is assumed.
%

% Anders Goran, Tomlab Optimization Inc, E-mail: anders@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Nov 11, 2002.    Last modified Jun 07, 2007

function z = ktr_dcS(x, Prob)

% Call nlp_dc to get constraint gradient
dcS = feval('nlp_dc',x, Prob);

if isempty(dcS)
   if isempty(Prob.ConsPattern)
      z = zeros(Prob.N*Prob.mNonLin,1);
   else
      z = zeros(nnz(Prob.ConsPattern),1);
   end
   return
end

if isempty(Prob.ConsPattern)
   if issparse(dcS)
      z=full(dcS);
      z=z(:);
   else
      z=dcS(:);
   end
else
   if ~isempty(Prob.ConsIdx)
      z = full(dcS(Prob.ConsIdx));
   else
      z = full(dcS(:));
   end
end

% MODIFICATION LOG:
%
% 021111 ango Wrote function, based on nlp_cdcS.m
% 021125 ango Complete rewrite, to use nlp_dc for value.
% 021223 hkh  Change sparse handling similar to nlp_cdcS
% 020206 hkh  Use name filterSQP, filterSQPTL
% 030429 ango This file was adopted by OQNLP, as it is no longer used by the DUNDEE solvers
% 030429 ango OQNLP needs untransposed constraint gradient
% 040602 ango Returns correct shape z regardless of Prob.ConsIdx present or not. 
% 070607 ango Better handling if nlp_dc returns empty