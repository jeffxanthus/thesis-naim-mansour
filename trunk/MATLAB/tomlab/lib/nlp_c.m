% TOMLAB gateway routine for computation of the constraint vector
%
% nlp_c calls the routine Prob.FUNCS.c
% either as c=feval(Prob.FUNCS.c, x) or
%           c=feval(Prob.FUNCS.c, x, Prob,varargin{:})
% depending on the number of inputs
%
% The global counter variable n_c is incremented
%
% function c = nlp_c(x, Prob, varargin)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written Oct 10, 1998. Last modified Dec 12, 2006.

function c = nlp_c(x, Prob, varargin)

global n_c NARG 
global mad_c
global NLP_xc NLP_c  % Communication nlp_c/dc

x=x(:);

if Prob.simType > 0
   [f,c] = sim_fc(x, Prob, varargin{:} );
   return
end

N = min(length(x),Prob.N);

if ~isempty(NLP_xc)
   if length(x)~=length(NLP_xc)
      NLP_xc=[];
   elseif all(x==NLP_xc)
      c=NLP_c;
      return
   end
end

n_c=n_c+1;

Func = Prob.FUNCS.c;

if isempty(Func)
   NLP_xc=[];
   NLP_c=[];
   c=zeros(0,1);
   return
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(4);
   end
   if Prob.ADCons == 1
      if p > 2
         %mad_c=feval(Func,fmad(x(1:N),speye(length(x(1:N)))),Prob,varargin{:});
         mad_c=feval(Func, fmad(x(1:N),speye(N)),Prob,varargin{:});
      elseif p > 1
         mad_c=feval(Func, fmad(x(1:N),speye(N)), Prob);
      else
         mad_c=feval(Func, fmad(x(1:N),speye(N)));
      end
      c=getvalue(mad_c);
   else
      if p > 2
         c=feval(Func, x(1:N), Prob, varargin{:});
      elseif p > 1
         c=feval(Func, x(1:N), Prob);
      else
         c=feval(Func, x(1:N));
      end
      c = full(c(:));
   end
end

NLP_xc=x;
NLP_c=c;

% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation.
% 981018   hkh   Delete ctrl vector for scaling
% 981029   hkh   Communicating NLP_x,f,c,dc to be used in numerical diffs
% 981102   hkh   Just using NLP_c and NLP_xc.
% 981120   hkh   Check if x and NLP_xc has same length
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
% 011109   hkh   Safeguard, always make c a column vector
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020413   hkh   Just send Prob.N first x variables
% 020429   hkh   length(x) < Prob.N if nnJac < nnObj for SOL solvers, use min
% 021230   hkh   Safe guard c=full(c), if user gives c as sparse
% 030524   hkh   Add handling of simulation problems
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 040318   hkh   Move sim_fc first, before check on NLP_xc
% 040409   hkh   More efficient call to MAD
% 040526   hkh   Use 1:N in Func calls, e.g.  c=feval(Func, x(1:N), Prob);
% 040901   med   getvalue lower case
% 060814   med   FUNCS used for callbacks instead
% 061212   med   ADMAT removed