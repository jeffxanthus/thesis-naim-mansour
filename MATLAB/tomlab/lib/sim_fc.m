% function [f,c] = sim_fc(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of function values f(x)
% together with constraint values c(x)
%
% sim_fc calls the routine Prob.FUNCS.fc
% either as [f,c]=feval(Prob.FUNCS.fc, x) or
%           [f,c]=feval(Prob.FUNCS.fc, x, Prob or
%           [f,c]=feval(Prob.FUNCS.fc, x, Prob,varargin)
% depending on the number of inputs
%
% The global counters n_f and n_c are incremented

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 24, 2003.   Last modified Aug 14, 2006.

function [f,c] = sim_fc(x, Prob, varargin)

global n_f n_c BUILDP NARG 

% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f  
global NLP_xc NLP_c

x=x(:);

if ~isempty(NLP_x) 
   % No test on NLP_xc, because it should be == NLP_x
   if length(x)~=length(NLP_x) 
      NLP_x  = [];
      NLP_xc = [];
   elseif all(x == NLP_x) 
      if ~isempty(NLP_f)
         f      = NLP_f;
         c      = NLP_c;
         return
      end
      %NLP_xc = x;
   end
end
n_f=n_f+1;
n_c=n_c+1;

Func = Prob.FUNCS.fc;

if isempty(Func)
   NLP_x=[];
   NLP_f=[];
   NLP_xc=[];
   NLP_c=[];
   f=NaN;
   c=NaN;
   return
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(10);
   end
   BUILDP0=BUILDP;
   BUILDP=[];
   if p > 2
      [f,c]=feval(Func, x, Prob, varargin{:} );
   elseif p > 1
      [f,c]=feval(Func, x, Prob);
   else
      [f,c]=feval(Func, x);
   end
   BUILDP=BUILDP0;
end
if isnan(f), f= 1E5; end
c(isnan(c))=1E5;
c=c(:);

NLP_x=x;
NLP_f=f;
NLP_xc=x;
NLP_c=c;

%if isnan(f)
%   error('sim_f: The user function returned NaN. Safeguard code!');
%end

% MODIFICATION LOG
%
% 030524  hkh   Written
% 040304  hkh   Safeguard for NaN
% 040318  hkh   Check for empty f if NLP_x==x, but most likely bug fixed now
% 040506  hkh   Safeguard c, changing to a column vector
% 060814  med   FUNCS used for callbacks instead