% function [X,norm,nCon] = constrainedLHD(norm,nSample,Prob,nTrial,nCon,PriLev)
%
% Use available maximin LHDs to find nSample feasible points.
%
% INPUT:
% norm      Choose from Maximin LHDs with specified norm
% nSample   Try to find nSample feasible points
% Prob      Problem structure
% nTrial    Maximum number of iterations
% PriLev    PriLev = 0, no output
%
% OUTPUT:
% X         Sample point matrix, d x M
% norm      Norm used in the final design

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Feb 2, 2008.    Last modified Sep 18, 2009.

function [X,norm,nCon] = constrainedLHD(norm,nSample,Prob,nTrial,nCon,PriLev)

if nargin < 6
   PriLev = 1;
end
if nargin < 5
   nCon = 0;
end
if nargin < 4
   MAXITER = 100;
elseif isempty(nTrial)  || ~isfinite(nTrial)
   MAXITER = 100;
else
   MAXITER = max(100,nTrial);
end

dim = Prob.N;
x_L = Prob.x_L;
x_U = Prob.x_U;
x_D = x_U - x_L;

dLin = Prob.mLin;
dCon = Prob.mNonLin;

M = dLin+dCon;
N = nSample + M;

iter = 1;
add1 = 0;
sub1 = 0;
plus = 0;
minus= 0;
absN = inf;

plusX = [];
XX = [];
nF = 0;
while iter <= MAXITER
   X  = getLHD(norm,dim,N);
   if isempty(X)
      % Use best found design
      if PriLev > 0
         fprintf(' Big enough design not available.')
      end
      if dim == 2
         norm = 1;
         if PriLev > 0
            fprintf(' Switch to 1-norm.\n')
         end
         [X,norm,nCon] = constrainedLHD(norm,nSample,Prob,MAXITER,nCon,PriLev);
         return
      elseif norm ~= 4
         norm = 4;
         if PriLev > 0
            fprintf(' Switch to AE-norm.\n')
         end
         [X,norm,nCon] = constrainedLHD(norm,nSample,Prob,MAXITER,nCon,PriLev);
         return
      else
         X = XX;
         if PriLev > 0
            fprintf(['\n Break instead with ' int2str(nF) ' feasible points.\n'])
         end
         return
      end
   else
      X = repmat(x_L,1,N) + repmat(x_D,1,N).*(X/(N-1));
   end
   [X,nF,nCon] = getFeasiblePoints(X,Prob,nCon,x_D,x_U);
   XX = X;
   n  = nSample - nF;
   
   if n > 0
      plus = 1;
   else
      minus = 1;
   end
   if abs(n) < absN
      absN = abs(n);
   elseif plus  &&  minus
      break
   end
   if n == 1
      if sub1 == 1
         break
      end
      add1 = 1;
   end
   if n == -1
      if add1 == 1
         break
      end
      sub1 = 1;
   end
   
   if n == 0
      return
   elseif n > 0
      % Haven't found enough feasible points, increase N.
      N = N + max(floor(n/dim),1); %max(floor(log2(n)),1);  %n;
   elseif n < 0
      % Have found too many feasible points, decrease N.
      N = N - max(floor(log2(-n)),1);
      plusX = X;
   end
   iter = iter+1;
end

if ~isempty(plusX)
   X = plusX;
   nF = size(X,2);
end

if PriLev > 0
   if iter > MAXITER
      disp([' MAXITER reached, cannot find exactly ' int2str(nSample) ' points.'])
   else
      disp([' Alternating between two designs, cannot find exactly ' int2str(nSample) ' points.'])
   end
   disp([' Break instead with ' int2str(nF) ' feasible points.'])
end




function X = getLHD(norm,dim,nSample)
%    function X = getLHD(norm,dim,nSample)
%
%     Norm     Type               dim    nSample
%    ----------------------------------------------
%       1   Maximin L1-norm         2     2-inf
%                                   3     2-16
%       2   Maximin L2-norm       2-4     2-300
%                                5-10     2-100
%       3   Maximin Inf-norm        2     2-inf
%                                   3     2-17,27-28,64-65
%                                 4-8     2-10
%       4   Audze-Eglais (AE)    2-10     2-100
%       5   Minimax Designs         2     2-27

X = [];
if ( dim < 2   ||   nSample < 2   ||   ~ismember(norm,1:5) )                                       || ...
   ( norm == 1    &&  ( dim >  3  ||  ( dim == 3  &&  nSample > 16 ) ) )                           || ...
   ( norm == 2    &&  ( dim > 10  ||  ( ismember(dim,2:4)   &&  nSample > 300 ) ...
                                  ||  ( ismember(dim,5:10)  &&  nSample > 100 ) ) )                || ...
   ( norm == 3    &&  ( dim >  8  ||  ( dim == 3  &&  ~ismember(nSample,[2:17,27:28,64:65]) ) ...
                                  ||  ( ismember(dim,4:8)   &&  nSample > 10 ) ) )                 || ...
   ( norm == 4    &&  ( dim > 10  ||  nSample > 100 ) )                                            || ...
	( norm == 5    && ~( dim == 2  &&  nSample <= 27 ) )
   %disp('No such design available. Returns X = [];')
   return
end

if norm == 2
   X = getMaxiMinLHDnorm2(dim,nSample);
elseif norm == 4
   X = getLHDnormAE(dim,nSample);
else
   X = getMaxiMinLHDnorm135(norm,dim,nSample);
end



function [XX,nF,nCon] = getFeasiblePoints(XX,Prob,nCon,x_D,x_U)

dLin = Prob.mLin;
dCon = Prob.mNonLin;
bTol = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol = Prob.optParam.cTol;      % Constraint feasibility tolerance
IntVars = Prob.MIP.IntVars;

% If any Integer Constraints, round LHD.
if ~isempty(IntVars)
   XX = daceInts(XX,x_D,x_U,IntVars);
end

% Check if any linear constraints, find feasible points.
if dLin > 0
   fIdx = all( [-Prob.A;Prob.A]*XX <= repmat([-Prob.b_L;Prob.b_U],1,size(XX,2))+bTol );
   XX = XX(:,fIdx);
end
nF = size(XX,2);

% Check if any nonlinear constraints.
fIdx = true(1,nF);
if dCon > 0
   % loop over linear feasible points and check which fulfill all constrints
   for i = 1:nF
      %C = nlp_c(O(:,i), Prob, varargin{:});
      %nCon = nCon + 1;
      C = nlp_c(XX(:,i), Prob);
      nCon = nCon + 1;
      cErr = sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
      if cErr > 0
         fIdx(i) = 0;
      end
   end
   XX = XX(:,fIdx);
   nF = size(XX,2);
end

% MODIFICATION LOG:
%
% 080707  nhq  Written
% 080716  nhq  Added PriLev and nCon
% 090918  hkh  In daceInts, generate with equal probability the integer values
% 090918  hkh  Change call to external daceInts
% 090918  hkh  Change call to getFeasiblePoints, change comments
