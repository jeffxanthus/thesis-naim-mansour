% function pbuild(x,f,LS)
%
% Build matrix with solver search directions from iteration points x
% or just the points tried (non-linesearch algorithms)
%
% f is the function value f(x), and LS is a flag if a nonlinear least
% squares problem is solved
%
% IF global BUILDP == 1 then
%    determine line search steps
%    determine limits for the iterations, used in default plotting
% IF global BUILDP == 2 then
%    save all points in p_dx, alphaV is then empty
%    determine limits for the iterations, used in default plotting
%
% Global variables used:
%   p_dx:	Matrix with all the search directions
%   alphaV:     Vector with the alfa steps (line search lengths)
%   X_min:	Smallest x_k(i)-value
%   X_max:	Biggest x_k(i)-value
%   F_X:        Matrix with function values. Rows: [Iteration_no f(x)]
%
% Internal global variables:
%   X_OLD       Last known base point x_k
%   X_NEW       Last x point in line search. Possible new x_k
%   pLen        Number of iterations so far

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Nov 6, 1997.    Last modified Jul 17, 2009.

function pbuild(x,f,LS)

if nargin < 3
   LS=0;
end

global BUILDP

if isempty(BUILDP), return, end
if BUILDP==0, return, end

global p_dx alphaV X_min X_max F_X
global X_OLD X_NEW pLen

n=length(x);
% Special case, when putting pbuild call in nlp_r for LS problems
if LS
   f=0.5*(f(:)'*f(:));
end

if isempty(p_dx)
   p_dx=zeros(n,1);
   pLen=0;
   %alphaV=ones(300,1);
   alphaV=zeros(0,2);
   X_OLD=[];
   X_NEW=[];
   F_X=[0 f];
end

if isempty(X_OLD)
   X_OLD=x;
   X_min = x;
   X_max = x;
   if BUILDP > 1
      pLen=1;
      p_dx(:,pLen)=x;
   end
else
   if pLen==0
      %if all(x==X_OLD)
      %   disp('all x equal')
      %   pause
      %end
      if all(x==X_OLD), return; end  % Equal e.g. if estimating gradients
      pLen=1;
      p_dx(:,pLen)=x-X_OLD;
      X_NEW=x;
   elseif BUILDP > 1
      pLen=pLen+1;
      p_dx(:,pLen)=x;
      %X_OLD=X_NEW;
      %X_NEW=x;
   else
      % Check if line search or new step
      pNew=x-X_OLD;
      pOld=p_dx(:,pLen);
      i=find(pOld~=0);
      alpha=pNew(i(1))/pOld(i(1));
      if all(abs(alpha*pOld - pNew) < 1E-11) & alpha < 200
         % Assume the same direction, i.e. a line search step
         % disp([pLen alpha]);
         m=size(alphaV,2);
         alphaV=[alphaV;[pLen,alpha,zeros(1,m-2)]];
         X_NEW=x;
      else
         % Assume new direction
         pLen=pLen+1;
         p_dx=[p_dx,x-X_NEW];
         X_OLD=X_NEW;
         X_NEW=x;
      end
   end
   F_X=[F_X;[pLen f]];
   X_min=min(X_min,x);   % Determine plotting region
   X_max=max(X_max,x);
end

% MODIFICATION LOG
%
% 981022  hkh  Added 3rd parameter LS for least squares.
% 030325  hkh  Change tolerance from 1E-14 to 1E-11 and do all(...) in test
% 040728  med  MIDEVA pragma removed
