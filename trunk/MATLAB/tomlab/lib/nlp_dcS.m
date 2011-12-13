% function z = nlp_dcS(x, Prob)
%
% nlp_dcS calls the TOMLAB gateway function nlp_dc,
% that computes the gradient dc for all constraints at x
%
% nlp_dcS is used when calling MINLPBB and FILTERSQP in sparse version
%
% nlp_dcS returns the static non-zero values of dcS in the vector z
%
% dcS could be either dense or sparse.
% The non-zero pattern is described in Prob.ConsPattern (should be sparse)
% If isempty(Prob.ConsPattern), then a dense matrix is assumed.
%
% To convert from dynamic sparse Matlab format to static sparse Fortran 0-1
% format determined by Prob.ConsPattern, a linear index vector is computed
% in minlpbbTL.m and filterSQPTL.m, and sent in Prob.ConsIdx.
% The 2 lines of code needed are commented below (find/sub2ind), for users
% that call MINLPBB and filterSQP directly without using minlpbbTL.m
% and filterSQPTL.m

% Anders Goran, Tomlab Optimization Inc, E-mail: anders@tomopt.com
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 11, 2002.    Last modified Jun 2, 2004

function z = nlp_dcS(x, Prob)

% Call nlp_dc to get constraint gradient
dcS = feval('nlp_dc',x, Prob);

if isempty(dcS)
   z=[];
else
   if isempty(Prob.ConsPattern)
      if issparse(dcS)
         z=full(dcS);
         z=z(:);
      else
         z=dcS(:);
      end
   else
      if issparse(dcS) | all(size(dcS) > 1) 
         if ~isempty(Prob.ConsIdx)
            z = full(dcS(Prob.ConsIdx));
         else
            z = full(dcS(:));
         end
      else
         z = dcS(:);
      end
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