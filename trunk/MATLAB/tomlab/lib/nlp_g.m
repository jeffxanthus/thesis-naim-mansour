% function g = nlp_g(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the gradient vector
%
% nlp_g calls the routine Prob.FUNCS.g
% either as g=feval(Prob.FUNCS.g, x) or
%           g=feval(Prob.FUNCS.g, x, Prob,varargin)
% depending on the number of inputs
%
% The global counter variable n_g is incremented

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 10, 1998.   Last modified July 24, 2011.

function g = nlp_g(x, Prob, varargin)

global n_g BUILDP NARG PartSep
global mad_f mad_g
% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f NLP_pSepIndex
% Communication nlp_g/H
global NLP_g NLP_xg

if Prob.GATEF == 2
   g = feval(Prob.FUNCS.g, x, Prob);
   return
elseif Prob.GATEF == 3
   g = fdng(x(:), Prob, [], []); % FD Algorithm
   return
elseif Prob.simType > 0
   g = sim_gdc(x, Prob, varargin{:} );
   return
end

n_g=n_g+1;

Func = Prob.FUNCS.g;
N    = min(length(x),Prob.N);

if Prob.ADObj == 1 & ~strcmp(Func,'ls_g')
   if all(NLP_x == x)
      g=squeeze(getderivs(mad_f));
      g=g(:);
      if isempty(g)
         g = zeros(N,1);
      end
   else
      if isempty(NARG)
         p = xnargin(Func);
      else
         p = NARG(1);
      end
      Func  = Prob.FUNCS.f;
      fdvar = Prob.cols;
      if isempty(fdvar) | length(fdvar) == N
         Z = speye(N);
      else
         z = zeros(N,1);
         z(Prob.cols) = 1;
         Z = spdiags(z,0,N,N);
      end
      if p > 2
         madf=feval(Func, fmad(x(1:N),Z), Prob,varargin{:});
      elseif p > 1
         madf=feval(Func, fmad(x(1:N),Z), Prob);
      else
         madf=feval(Func, fmad(x(1:N),Z));
      end
      g=squeeze(getderivs(madf));
      g=g(:);
   end
elseif Prob.ADObj == -1
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(2);
   end
   if p > 2
      mad_g=feval(Func, fmad(x(1:N),speye(N)), Prob, varargin{:});
   elseif p > 1
      mad_g=feval(Func, fmad(x(1:N),speye(N)), Prob);
   else
      mad_g=feval(Func, fmad(x(1:N),speye(N)));
   end
   g=getvalue(mad_g);
   g=g(:);
elseif isempty(Func) | (Prob.NumDiff > 0 & isempty(Prob.FUNCS.r))
   % Here we should call a numerical difference routine
   % Check if we have computed f(x) already (the normal case)
   fx=[];
   if PartSep
      if ~isempty(NLP_x)
         if NLP_pSepIndex == Prob.PartSep.index
            if all(x == NLP_x)
               fx=NLP_f;
            end
         end
      end
   elseif ~isempty(NLP_x)
      if all(x == NLP_x)
         fx=NLP_f;
      end
   end
   BUILDP0=BUILDP;
   BUILDP=[];

   % if Prob.g_k defined Only NaN values will be estimated
   switch Prob.NumDiff
    case 1
      g=fdng(x(1:N), Prob, Prob.g_k, fx, varargin{:}); % FD Algorithm
    case {2,3,4}
      g=fdng2(x(1:N), Prob, Prob.g_k, fx, varargin{:}); % Using splines

      %if exist('csapi') & exist('csaps') & exist('spaps')
      %   g=fdng2(x, Prob, Prob.g_k, fx, varargin{:}); % Using splines
      %else
      %   fprintf('Can not find directory SPLINES\n')
      %   fprintf('Have you got a license for Spline Toolbox?\n')
      %   fprintf('Using finite differece routine fdng.m\n')
      %   g=fdng(x, Prob, Prob.g_k, fx, varargin{:}); % FD Algorithm
      %end
    case 5
      % Using complex number technique
      g=fdng3(x(1:N), Prob, Prob.g_k, fx, varargin{:});
    case 7
      g=fdng_par(x(1:N), Prob, Prob.g_k, fx, varargin{:}); % FD parfor Algorithm
    otherwise
      g=fdng(x(1:N), Prob, Prob.g_k, fx, varargin{:}); % FD Algorithm
   end
   NLP_x = x;
   NLP_f = fx;
   BUILDP=BUILDP0;
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(2);
   end
   if p > 2
      g=feval(Func, x(1:N), Prob, varargin{:});
   elseif p > 1
      g=feval(Func, x(1:N), Prob);
   else
      g=feval(Func, x(1:N));
   end
   g = full(g(:));
   if Prob.CheckNaN > 0 | isempty(g)
      if isempty(g), g=NaN*ones(N,1); end
      if any(isnan(g)) % There are elements set to NaN, to be estimated
         fx=[];
         if PartSep
            if ~isempty(NLP_x)
               if NLP_pSepIndex == Prob.PartSep.index
                  if all(x == NLP_x)
                     fx=NLP_f;
                  end
               end
            end
         elseif ~isempty(NLP_x)
            if all(x == NLP_x)
               fx=NLP_f;
            end
         end
         BUILDP0=BUILDP;
         BUILDP=[];
         switch abs(Prob.CheckNaN)
          case 1
            g=fdng(x(1:N), Prob, g, fx, varargin{:}); % FD Algorithm
          case {2,3,4}
            g=fdng2(x(1:N), Prob, g, fx, varargin{:}); % Using splines
          case 5
            g=fdng3(x(1:N), Prob, g, fx, varargin{:});% Using complex numbers
          otherwise
            g=fdng(x(1:N), Prob, g, fx, varargin{:}); % FD Algorithm
         end
         NLP_x = x;
         NLP_f = fx;
         BUILDP=BUILDP0;
      end
   end
end
NLP_xg = x;
NLP_g  = g;

% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation
% 981026   hkh   Added numerical differentiation using fdng
% 981027   mbk   Added numerical differentiation using fdng2
% 981027   hkh   Unnecessary return. No Prob.ctrl used
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
% 981029   hkh   Communicating NLP_x,f,c,dc to be used in numerical diffs
%                If NLLS problem, numerical differences should be done for
%                the Jacobian, not the gradient. Check if isempty(residFunc).
% 981102   hkh   Must check if partially separable function. Add to global.
% 981110   mbk   Changed "elseif Prob.NumDiff == 2" to
%                "elseif (Prob.NumDiff >= 1) & (Prob.NumDiff <= 4)".
% 981123   mbk   Changed "elseif (Prob.NumDiff >= 1) & (Prob.NumDiff <= 4)"
%                to "elseif (Prob.NumDiff > 1) & (Prob.NumDiff <= 4)".
% 981127   hkh   Change to safeguarded check on PartSep
% 990307   hkh   Add the new complex number diff technique in fdng3.m
% 001212   hkh   Use Prob.NumDiff > 5 for solver estimate of gradient
% 011201   hkh   Change BUILDP to avoid saving of numerical gradient when
%                determining search direction automatically
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 020416   hkh   Do not set g empty if NumDiff == 6, will not be called by SOL
% 020528   hkh   More efficient handling of partially separable functions
% 021223   hkh   Safe guard g=full(g), if user gives g as sparse.
% 030127   hkh   Check NaN elements if Prob.CheckNaN > 0, estimate numerically
% 030210   hkh   Safe guard: Change g=full(g); to g=full(g(:));
% 030522   hkh   Safe guard: Change g=full(g); to g=full(g(:));
% 030524   hkh   Add handling of simulation problems
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 040407   hkh   Always call fdng2, no check on spline TB routines
% 040407   hkh   Add g from Prob.g_k in call to FD, only estimate NaN values
% 040409   hkh   More efficient call to MAD
% 040526   hkh   Use x(1:N) in all function calls
% 040901   med   getvalue lower case
% 050616   hkh   Avoid isfield on Prob.g_k, always defined in Prob
% 050819   hkh   For LS, AD must be applied in nlp_J, otherwise weights wrong
% 060218   hkh   Check if isempty and estimate gradient numerically if so
% 060814   med   FUNCS used for callbacks instead
% 061212   med   ADMAT removed
% 080608   hkh   Add fast calls if GATEF > 1
% 110724   hkh   Add call to fdng_par (parfor) if NumDiff = 7
