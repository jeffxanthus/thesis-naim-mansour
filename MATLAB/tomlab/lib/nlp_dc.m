% TOMLAB gateway routine for computation of the constraint Jacobian matrix
%
% nlp_dc calls the routine Prob.FUNCS.dc either as 
%           dc=feval(Prob.FUNCS.dc, x) or
%           dc=feval(Prob.FUNCS.dc, x, Prob) or
%           dc=feval(Prob.FUNCS.dc, x, Prob, varargin{:})
% depending on the number of inputs in Prob.FUNCS.dc
%
% function dc = nlp_dc(x, Prob, varargin)
%
% dc is a mL x n matrix, where mL = number of nonlinear constraints
%                               n = number of variables
%
% The global counter variable n_dc is incremented
% If Prob.ConsDiff > 0, dc is estimated numerically
% If Prob.CheckNaN > 0, NaN elements in dc are estimated numerically
%
% Numerical method is dependent on Prob.ConsDiff and Prob.CheckNaN

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 10, 1998.   Last modified July 22, 2011.

function dc = nlp_dc(x, Prob, varargin)

global n_dc NARG
global mad_c mad_dc
global NLP_xc NLP_xdc NLP_c NLP_dc  % Communication nlp_c/dc

if Prob.simType > 0
   [g,dc] = sim_gdc(x, Prob, varargin{:} );
   return
end

x=x(:);

N = min(length(x),Prob.N);

if ~isempty(NLP_xdc)
   if length(x)~=length(NLP_xdc)
      NLP_xdc=[];
   elseif all(x==NLP_xdc)
      dc=NLP_dc;
      return
   end
end

n_dc=n_dc+1;

Func = Prob.FUNCS.dc;

if Prob.ADCons == 1
   if ~isempty(NLP_xc)
      if all(x==NLP_xc)
         cx=NLP_c;
      else
         cx=[];
      end
   else
      cx=[];
   end
   if isempty(cx)
      if isempty(NARG)
         p = xnargin(Func);
      else
         p = NARG(4);
      end
      Func  = Prob.FUNCS.c;
      fdvar = Prob.FDVar;
      if fdvar == 0 | length(fdvar) == N
         Z = speye(N);
      else
         z = zeros(N,1);
         z(fdvar) = 1;
         Z = spdiags(z,0,N,N);
      end
      if p > 2
         madc=feval(Func, fmad(x(1:N),Z), Prob,varargin{:});
      elseif p > 1
         madc=feval(Func, fmad(x(1:N),Z), Prob);
      else
         madc=feval(Func, fmad(x(1:N),Z));
      end
      dc=getinternalderivs(madc);
   else
      dc=getinternalderivs(mad_c);
   end
elseif Prob.ADCons == -1  
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(5);
   end
   if p > 2
      mad_dc=feval(Func, fmad(x(1:N),speye(N)), Prob, varargin{:});
   elseif p > 1
      mad_dc=feval(Func, fmad(x(1:N),speye(N)), Prob);
   else
      mad_dc=feval(Func, fmad(x(1:N),speye(N)));
   end
   dc=getvalue(mad_dc);
elseif isempty(Func) | Prob.ConsDiff > 0
   if isempty(Prob.FUNCS.c)
      dc=zeros(0,N);
      NLP_xdc = [];
      NLP_dc  = [];
      return
   else
      % Numerical difference computation of dc
      if ~isempty(NLP_xc)
         if all(x==NLP_xc)
            cx=NLP_c;
         else
            cx=[];
         end
      else
         cx=[];
      end
      ConsDiff = Prob.ConsDiff;
      if any(ConsDiff == [2 3 4])
         [dc,cx]=FDJac2(x(1:N), Prob, 'nlp_c', cx, varargin{:}); % Using splines
      elseif ConsDiff == 5
         % Using complex number technique
         [dc,cx]=FDJac3(x(1:N), Prob, 'nlp_c', cx, varargin{:}); 
      elseif ConsDiff == 7
         [dc,cx]=FDJac_par(x(1:N), Prob, 'nlp_c', cx, varargin{:}); % FD parfor Algorithm
      else
         [dc,cx]=FDJac(x(1:N), Prob, 'nlp_c', cx, varargin{:}); % FD Algorithm
      end   
      NLP_xc = x;
      NLP_c  = cx;
   end
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(5);
   end

   if p>2
      dc = feval(Func, x(1:N), Prob, varargin{:});
   elseif p==2
      dc = feval(Func, x(1:N), Prob);
   else
      dc = feval(Func, x(1:N));
   end
   if Prob.CheckNaN > 0
      [iN,jN,dcN] = find(isnan(dc));
      if ~isempty(iN)    % There are elements set to NaN, to be estimated
         Prob.ConsPattern = sparse(iN,jN,dcN,size(dc,1),N);
         if ~isempty(NLP_xc)
            if all(x==NLP_xc)
               cx=NLP_c;
            else
               cx=[];
            end
         else
            cx=[];
         end
         if any(Prob.CheckNaN == [2 3 4])
            [dcN,cx]=FDJac2(x(1:N), Prob, 'nlp_c', cx, varargin{:}); 
         elseif Prob.CheckNaN == 5
            % Using complex number technique
            [dcN,cx]=FDJac3(x(1:N), Prob, 'nlp_c', cx, varargin{:}); 
         else
            [dcN,cx]=FDJac(x(1:N), Prob, 'nlp_c', cx, varargin{:});%FD Algorithm
         end   
         % Merge analytic and numerical dc
         dc(isnan(dc)) = dcN(isnan(dc));
         NLP_xc = x;
         NLP_c  = cx;
      end
   end
end

NLP_xdc=x;
NLP_dc=dc;
% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation 
% 981028   hkh   Add code for numerical differences using splines
% 981102   hkh   Just using NLP_c and NLP_xc. Change empty setting of dc
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
% 001114   hkh   Use Prob.ConsDiff instead of NumDiff for explicit control
% 001212   hkh   Use Prob.ConsDiff > 5 for solver estimate of gradient
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020413   hkh   Just send the first Prob.N x variables to Func
% 020416   hkh   Do not set dc empty if ConsDiff == 6, will not be called by SOL
% 020429   hkh   length(x) < Prob.N if nnJac < nnObj for SOL solvers, use min
% 030127   hkh   Check for NaN elements, estimate numerically
% 030524   hkh   Add handling of simulation problems
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 040331   hkh   Use 2nd parameter cx from FDJac-routines, dc could get []
% 040331   hkh   Save NLP_dc, and use if x==NLP_xdc
% 040407   hkh   Always call FDJac2, no check on spline TB routines
% 040409   hkh   More efficent handling of repeated nlp_d2c calls for AD-num.der
% 040526   hkh   Use x(1:N) in all function calls
% 040901   med   getvalue lower case
% 060814   med   FUNCS used for callbacks instead
% 061212   med   ADMAT removed
% 110724   hkh   Added option ConsDiff=7, using parfor in FDJac_par
