% nlp_cdcS.m
%
% function [Mode, c, z]=nlp_cdcS(x, Prob, Mode, nState)
%
% Sparse version of nlp_cdc.m
%
% nlp_cdcS calls the TOMLAB gateway function nlp_c,
% which evaluates the constraints c at x for the test problem P (Prob.P).
%
% It also calls the TOMLAB gateway function nlp_dc,
% which computes the gradient for all constraints at x, dc, for test problem P
%
% Mode and nState is sent as Prob.Mode, Prob.nState to nlp_c and nlp_dc.
%
% Mode = 0 Assign function values
% Mode = 1 Assign known derivatives, unknown set as -11111 (=missing value)
% Mode = 2 Assign function and known derivatives, unknown derivatives -11111
%
% nState = 1         First call
% nState = 0         Other calls
% nState = 2+Inform  Last call, see Inform parameter for the solver
%
% nlp_cdcS is used when calling MINOS and SNOPT in sparse version
%
% nlp_cdcS returns c and the static non-zero values of dcS in the vector z
%
% dcS could be either dense or sparse.
% The non-zero pattern is described in Prob.ConsPattern (should be sparse)
% If isempty(Prob.ConsPattern), then a dense matrix is assumed.
%
% To convert from dynamic sparse Matlab format to static sparse Fortran 0-1
% format determined by Prob.ConsPattern, a linear index vector is computed
% in minosTL.m and snoptTL.m, and sent in Prob.ConsIdx.
% The 2 lines of code needed are commented below (find/sub2ind), for users
% that call SNOPT and MINOS directly without using snoptTL.m and minosTL.m


% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Mar 5, 1997.    Last modified July 24, 2011.

function [Mode, c, z]=nlp_cdcS(x, Prob, Mode, nState)

if ~isempty(Prob.FUNCS.cdc)
   if Mode == 0
      [Mode,c]=feval(Prob.FUNCS.cdc,x,Prob,Mode,nState);
   else
      [Mode,c,dcS]=feval(Prob.FUNCS.cdc,x,Prob,Mode,nState);
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
            if issparse(dcS) | all(size(dcS) > 1 )
               z = full(dcS(Prob.ConsIdx));
               % This is self contained, when Prob.ConsIdx not computed:
               % [i,j]=find(Prob.ConsPattern);
               % ind = sub2ind(size(Prob.ConsPattern),i,j);
               % z = full(dcS(ind));
               % The old way of adding an epsilon to get the static picture
               %if nnz(Prob.ConsPattern) ~= nnz(dcS)
               %   dcS = dcS + eps*(Prob.ConsPattern - (dcS~=0));
               %end
               %z = dcS(:);
            else
               % If dcS is a vector, assume the user has put the non-zero 
               % elements in the vector properly compared to ConsPattern
               z = dcS(:);
            end
         end
      end
   end
elseif length(x) < Prob.N
   Prob.Mode   = Mode;
   Prob.nState = nState;
   c = nlp_c( [x;zeros(Prob.N-length(x),1)], Prob);
   if Mode > 0
      if Prob.ConsDiff < 6
         dcS=nlp_dc([x;zeros(Prob.N-length(x),1)], Prob); 
      else
         dcS = [];
      end
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
            if issparse(dcS) | all(size(dcS) > 1 ) | Prob.mNonLin == 1
               z = full(dcS(Prob.ConsIdx));
            else
               z = dcS(:);
            end
         end
         return
      end
   end
else
   Prob.Mode   = Mode;
   Prob.nState = nState;
   c = nlp_c( x, Prob);
   if Mode > 0
      if Prob.ConsDiff ~= 6
         dcS=nlp_dc(x, Prob); 
      else
         dcS = [];
      end
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
            if issparse(dcS) | all(size(dcS) > 1 ) | Prob.mNonLin == 1
               z = full(dcS(Prob.ConsIdx));
            else
               z = dcS(:);
            end
         end
         return
      end
   end
end

% MODIFICATION LOG:
%
% 010626 hkh Changing dynamic/static strategy to epsilon strategy 
% 020212 hkh Changing dynamic/static strategy to find/sub2ind/full
% 020301 hkh Safeguard against empty dcS returned
% 020409 hkh Adding Mode to Prob.Mode, nState to Prob.nState
% 020919 hkh No call to nlp_dc if ConsDiff == 6, because no check in nlp_dc
% 021223 hkh Assume user has defined non-zero elements according to pattern
%            in ConsPattern ONLY if dcS is a vector.
% 050616 hkh Avoid isfield, assume Prob.FUNCS.cdc is defined
% 060814 med FUNCS used for callbacks instead
% 091014 hkh Avoid crash when mNonLin==1 and dcS full, test mNonLin == 1
% 100928 hkh Add special code to handle length(x) < Prob.N
% 110724  hkh  Test on ConsDiff ~= 6, not < 6, to enable parfor options
