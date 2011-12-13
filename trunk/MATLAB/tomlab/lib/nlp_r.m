% function [rW, rWeight] = nlp_r(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the residual vector
%
% INPUT:
% x       Current iterate x, where the residual is vector function r=r(x).
% Prob    TOMLAB problem structure
%
% To compute r(x) nlp_r calls the routine Prob.FUNCS.r either as
%            r=feval(Prob.FUNCS.r, x) or
%            r=feval(Prob.FUNCS.r, x, Prob, varargin{:})
% depending on the number of inputs
%
% nlp_r computes a weighted residual, using Prob.LS.weightType/weightY/y.
%
% The global counter variable n_r is incremented
%
% Weighting is computed and saved in the global variable wLS,
% which is used by gateway routine nlp_J to compute Jacobian weighting
%
% r and x are stored in globals LS_r and LS_x, to avoid recomputing.
%
% OUTPUT:
% rW      Weighted residual vector, where rW=rWeight*r; r is residual.
% rWeight Computed and returned as a matrix if user calls with two arguments.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 10, 1998.  Last modified July 22, 2011.

function [rW, rWeight] = nlp_r(x, Prob , varargin)

global wLS
global n_r NARG
global LS_x LS_r LS_xJ LS_J
global SEP_z SEP_Jz
global mad_r
x=x(:);
N=Prob.N; % N is the number of parameters to be sent to user routines

if ~isempty(LS_x)
   if length(x)~=length(LS_x)
      LS_x=[];
   elseif all(x==LS_x)
      rW=LS_r;
      rWeight=[];
      return
   end
end

n_r=n_r+1;
Func = Prob.FUNCS.r;

if isempty(Func)
   LS_x=x;
   LS_r=[];
   LS_J=[];
   rW=[];
   rWeight=[];
   return
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(7);
   end
   if Prob.ADObj == 1
      if p > 2
         mad_r=feval(Func, fmad(x(1:N),speye(length(x(1:N)))),Prob,varargin{:});
      elseif p > 1
         mad_r=feval(Func, fmad(x(1:N),speye(length(x(1:N)))), Prob);
      else
         mad_r=feval(Func, fmad(x(1:N),speye(length(x(1:N)))));
      end
      r=getvalue(mad_r);
   elseif Prob.ADObj == 0
      if p > 2
         if Prob.LS.SepAlg
            parr=nargout(Func);
            if parr < 4 
               r = feval(Func, x(1:N), Prob);
               fprintf('nlp_r: ERROR!!! Flag Prob.LS.SepAlg=1, but ');
               fprintf('%s',Func);
               fprintf(' has only %d output parameters\n',parr);
               fprintf('Must return r, J, z and Jz, at least 4 parameters!\n');
               fprintf('Try to run ordinary nonlinear least squares.\n');
               z=[]; J=[]; Jz=[];
               Prob.LS.SepAlg=0;
            else
               [r J z Jz] = feval(Func, x(1:N), Prob, varargin{:} );
            end
         else
            r=feval(Func, x(1:N), Prob, varargin{:} );
         end
      elseif p == 2
         if Prob.LS.SepAlg
            parr=nargout(Func);
            if parr < 4 
               r = feval(Func, x(1:N), Prob);
               fprintf('nlp_r: ERROR!!! Flag Prob.LS.SepAlg=1, but ');
               fprintf('%s has only %d output parameters\n',Func,parr);
               fprintf('Must return r, J, z and Jz, at least 4 parameters!\n');
               fprintf('Try to run ordinary nonlinear least squares.\n');
               z=[]; J=[]; Jz=[];
               Prob.LS.SepAlg=0;
            else
               [r J z Jz] = feval(Func, x(1:N), Prob);
            end
         else
            r=feval(Func, x(1:N), Prob);
         end
      else
         if Prob.LS.SepAlg
            [r J z Jz] = feval(Func, x(1:N));
         else
            r=feval(Func, x(1:N));
         end
      end
   end
end

weightType=Prob.LS.weightType;

if weightType==0
   wLS=[];
elseif weightType==1
   y=Prob.LS.y;
   if isempty(y) || sum(y)==0
      wLS=[];
   else
      ix=find(y(:,1)~=0);
      wLS=zeros(size(y,1),1);
      wLS(ix)=1./abs(y(ix,1));
   end
elseif weightType==2
   wLS=Prob.LS.weightY;
elseif weightType==3
   % The weights are given by a call to a function with params (x,r,Prob)
   wFunc=Prob.LS.weightY;
   if exist(wFunc,'file') & ischar(wFunc)
      if xnargin(wFunc) > 2
         wLS=feval(wFunc,x(1:N),r,Prob);
      elseif xnargin(wFunc) > 1
         wLS=feval(wFunc,x(1:N),r);
      else
         wLS=feval(wFunc,x(1:N));
      end
   else
      fprintf('nlp_r: Error!!! Given weight routine in wFunc not found!\n'); 
      wLS=[];
   end
else
   wLS=[];
end
if ~isempty(wLS)
   if any(size(wLS)==1) & any(size(wLS)==length(r))
      wLS=wLS(:);
      rW=full(wLS.*r(:));
   elseif all(size(wLS)==length(r))  % Must be square weighting matrix
      rW=full(wLS*r(:));
   else  % Trouble!!!, ignore weights
      fprintf('nlp_r: Error!!! Given weight vector/matrix has wrong size!\n');
      wLS=[];
      rW=r(:);
   end
else
   rW=full(r(:));
end
if nargout > 1
   if isempty(wLS)
      rWeight=eye(length(r),length(r));
   elseif any(size(wLS)==1)
      rWeight=diag(wLS);
   else
      rWeight=wLS;
   end
end

if Prob.LS.SepAlg
   % Jz and J is already weighted in exp_r and expLS
   SEP_z=z; SEP_Jz=Jz;

   [rW Js Q R pRank] = Sep_rJ(x,rW,J,z,Jz,Prob);

   % Projected and weighted Jacobian now computed
   LS_xJ=x;
   LS_J=Js;
end
LS_x=x;
LS_r=rW;

pbuild(x,rW,1);

% MODIFICATION LOG
%
% 981011   hkh   Added automatic differentiation.
% 981018   hkh   Added weighted least squares handling
% 981023   hkh   Make the routine more general. Use globals.
%                Return rW if already computed (used by g and H routines).
% 981024   hkh   Add LS_xJ as global
% 981108   hkh   Add extra parameter z with the separable solution
% 981118   hkh   Call pbuild with rW, not r.
% 981120   hkh   Check if x and LS_x has same length
% 981126   hkh   Use xnargin as filter, to avoid bug in Matlab5.1
%                Add printing about the use of automatic differentiation
% 981210   hkh   weightType==1, Use 1/y, not sqrt(1/y).
% 990523   hkh   Dangerous change for minimax. Use N=Prob.N to determine how
%                many variables to send to user routine (to allow for extra
%                variables in x vector).
% 011109   hkh   Safeguard, always make residual column vector
% 020409   hkh   Use global NARG instead of calling xnargin every time
% 030808   ango  Correct comments
% 031201   hkh   Revising AD handling, new for MAD, changes for ADMAT
% 040427   hkh   Do full(rW) in case of sparse weight vector/matrix
% 040901   med   getvalue lower case
% 050414   ango  Fixed an incorrect fprintf statement
% 050801   med   isstr replaced by ischar
% 060814   med   FUNCS used for callbacks instead
% 061212   med   ADMAT removed
% 110722   hkh   Corrected comments about weighting
