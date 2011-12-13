% function = sepMIP_f(IntVars)
%
% Function computing MINLP/NLP/NLLS subproblem in separable MINLP approach
% when using genIP for IP part, called from solver sepMINLP.
%
% genIP ensures that the logical constraint
% lower number <= sum(integer selection variables) <= upper number
%
% Global variables used: allIP feasIP ProbIP fBest ResBest IntBest
%
% sepMIP_f checks that simple bounds on the integer variables are fulfilled
% (some variables might be fixed on bounds)
% sepMIP_f checks that any logical constraints are fulfilled
% sepMIP_f is setting up the shrinked NLP problem and calls nlpblend to
% solve it

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 2004-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Apr 10, 2004. Last modified July 25, 2011.

function sepMIP_f(IntVars)

ix               = IntVars(:);    % Indices for current selection of variables
% Print level for bounds and logical constraints on integer variables
PriLev = 1; 

global allIP feasIP ProbIP fBest ResBest IntBest
allIP    = allIP + 1;  % Global counter including skipped points
x_k      = zeros(ProbIP.N,1);
x_k(ix)  = 1;
if any(ProbIP.x_L > x_k | ProbIP.x_U < x_k)
   % Simple bounds violated
   if PriLev > 2
      fprintf('Skip. Simple bounds violated!')
      fprintf(' %d',ix)
      fprintf('\n')
   end
   return
end
AA = ProbIP.MIP.A;
if ~isempty(AA)
   % b_L = ProbIP.MIP.b_L;
   % b_U = ProbIP.MIP.b_U;
   % Ax  = sum(AA(:,IntVars)')';
   % x_k   = zeros(size(AA,2),1);
   % x_k(IntVars) = 1;
   Ax   = AA*x_k;
   if any(Ax > ProbIP.MIP.b_U | Ax < ProbIP.MIP.b_L) 
      if PriLev > 1
         fprintf('Skip. Linear logical constraints violated')
         fprintf(' %d',ix)
         fprintf('\n')
      end
      return
   end
end

feasIP = feasIP + 1;   % Global counter which NLP sub problem
if PriLev > 0
   fprintf('Problem %d out of ',feasIP);
   fprintf('%d ',allIP);
   xprinti(IntVars,'IV:')
end

if isempty(ProbIP.MIP.FUNCS)
   % Run full problem with fixed variables
   P  = ProbIP.Prob;    % MINLP/NLP subproblem in Prob sub structure
   IP = find(P.A(1,:) == 1);
   nI = length(IP);
   % Set IP variables as fixed
   P.x_L(IP) = x_k;
   P.x_U(IP) = x_k;
   P.x_0(IP) = x_k;
else
   % Run shrinked problem 
   P  = ProbIP.Prob;    % MINLP/NLP subproblem in Prob sub structure
   ix = ones(P.N,1);
   IV = zeros(P.N,1);
   IV(P.MIP.IntVars) = 1;
   IP = ProbIP.MIP.IP;
   ix(IP) = 0;
   IV(IP) = 0;
   % The variables that are 1 and not 0 in the deleted variable set IP
   P.MIP.IP = IP(IntVars);
   ixV = find(ix);
   IVV = find(IV);
   % The integer variables not part of the set IP
   P.MIP.IntVars = IVV;
   N     = length(ixV);
   P.N   = N;
   P.x_L = P.x_L(ixV);
   P.x_U = P.x_U(ixV);
   P.x_0 = P.x_0(ixV);
   if P.mNonLin > 0
      % Shrink ConsPattern matrix
      P.ConsPattern    = P.ConsPattern(:,ixV);
      P.ConIx          = findpatt(P.ConsPattern);
   end
   P.Pnr = P.P;
   P.P   = feasIP;
   P.P2  = allIP;
   i1 = ProbIP.MIP.i1(:,1);
   if isempty(i1)
      i1v = [];
   else
      i1v =ProbIP.MIP.i1(:,2);
   end
   i2 = ProbIP.MIP.i2;
   if ~isempty(i1)
      v = full(sum(P.A(i1,IP)'));
      %v
      %[P.x_L P.x_U]
      for i = 1:length(v)
          j = i1v(i);
          k = i1(i);
          P.x_L(j) = max(P.x_L(j),P.b_L(k)-v(i)/P.A(k,j));
          P.x_U(j) = min(P.x_U(j),P.b_U(k)-v(i)/P.A(k,j));
      end
      %[P.x_L P.x_U]
   end
   if ~isempty(i2)
      v = full(sum(P.A(i2,IP)'));
   end
   % must compute A*x for IP part, sum A(IP)
   if ~isempty(P.LINCON)
      % Linear constraints  for sub problem
      P.A     = P.LINCON.A;
      P.b_L   = P.LINCON.b_L-v;
      P.b_U   = P.LINCON.b_U-v;
      P.mLin  = size(P.A,1);
   else
      P.A     = [];
      P.b_L   = [];
      P.b_U   = [];
      P.mLin  = 0;
   end
   if isfield(ProbIP.MIP.FUNCS,'Init')
      P = feval(ProbIP.MIP.FUNCS.Init,P);
   end
   P.FUNCS = ProbIP.MIP.FUNCS;
end

R = tomRun(ProbIP.Solver,P,0);
f = R.f_k;
if R.ExitFlag == 0
   % PrintResult(R,2);
   if f < fBest
      fBest   = f;
      ResBest = R;
      IntBest = IntVars;
   end
end
if PriLev > 0
   fprintf('f_k %f fBest %f ',f,fBest);
   xprinti(IntVars,'IV:')
end
return
%================ END OF CODE; the rest is old

P                = ProbIP.Prob;    % NLP subproblem in Prob sub structure
%BLENDS           = P.user.BLENDS;
%RM               = P.user.RM;
% Create a smaller problem
n                = length(ix);
N                = n*BLENDS;      % Number of variables in the shrinked prob
ixF = [];
v   = RM;
j   = 0;
x_L = zeros(N,1);
x_U = ones(N,1);
for i = 1:BLENDS
    ixF = [ixF;v+ix];             % Indices in full problem
    v = v + RM;
    % Set lower and upper bounds for the selected variables
    x_L(j+1:j+n) = ProbIP.user.XX_L(ix,i);
    x_U(j+1:j+n) = ProbIP.user.XX_U(ix,i);
    j = j + n;
end
P.user.c = ProbIP.user.c(ixF);
%if ~isinf(fBest)                  % Simple check that cost is too high
%   c1 = sum(P.user.c*0.03);
%   c2 = min(P.user.c)*(1-0.03*n);
%   fNew = c1 + c2;
%   if fNew > fBest
%      if PriLev > 0
%         fprintf('SKIP Point %d, best f(x) possible: %f. ',fNew);
%         fprintf('Current best f(x):  %f\n',fNew);
%      end
%      return
%   end
%end

P.N              = N;
%P.x_L            = 0.03*ones(N,1);
%P.x_U            = ones(N,1);
P.x_L            = x_L;
P.x_U            = x_U;
if P.mNonLin > 0
   P.ConsPattern    = P.ConsPattern(:,ixF);
end
P.user.blx       = ProbIP.user.blx(ix);
P.user.tgNx      = ProbIP.user.tgNx(ix);
P.user.tgGx      = ProbIP.user.tgGx(ix);
P.user.safax     = ProbIP.user.safax(ix);
P.user.pufax     = ProbIP.user.pufax(ix);
P.user.h2mx      = ProbIP.user.h2mx(ix);
P.user.Nmodels   = ProbIP.user.Nmodels;
P.user.RM        = n;
P.user.RMIP      = 0;

P.P              = feasIP;
P.P2             = allIP;
P.Pnr            = inf;
% Test presolve analysis
%PriLevOpt        = P.PriLevOpt;
%P.PriLevOpt      = 1;
P.MIP.IntVars    = [];
% Comment this row if to avoid presolve analysis
%'presolve'
%P.PriLevOpt =1;
%P.A
%P.x_L
%P.x_U
if 0
   AA = P.A;
   BL = P.b_L;
   BU = P.b_U;
   P                = preSolve(P);
   P.A   = AA;
   P.b_L = BL;
   P.b_U = BU;
else
   %P                = preSolve(P);
end
%P.PriLevOpt =0;
%P.A
%P.b_L
%P.b_U
%P.x_L
%P.x_U
%keyboard

%P.PriLevOpt      = PriLevOpt;
%Generate a general x_0, with sum(x_0) = 1
x_L              = P.x_L;
x_U              = P.x_U;
%fprintf('x_U: ');
%fprintf('%4.2f ',x_U);
%fprintf('\n');
x_U(x_U < 0.03)  = 0;
j                = 0;
x_0              = zeros(N,1);
for i = 1:BLENDS
    i1               = find(x_L(j+1:j+n) ~= 0);
    i0               = find(x_L(j+1:j+n) == 0 & x_U(j+1:j+n) > 0);
    x_0(j+i1)        = 0.5*(x_U(j+i1)+x_L(j+i1));
    sumX             = sum(x_0(j+i1));
    z                = (1-sumX)/length(i0);
    if z < 0.03
       z             = 0.03;
       % Scale down the initial x_0
       %'Scale down'
       x_0(j+i1)     = (1-0.03*length(i0))/sumX*x_0(j+i1);
    end
    x_0(j+i0)        = z;
    j                = j + n;
end
P.x_0            = x_0;
%fprintf('x_0: ');
%fprintf('%4.2f ',x_0);
%fprintf('\n');
[f,Result] = nlpblend(ix, P);
if ~isinf(f)
   global xFeas fFeas cFeas nFeas iFeas
   nFeas = nFeas+1;
   xFeas(:,nFeas) = Result.x_k;
   iFeas(:,nFeas) = ix;
   fFeas(nFeas) = f;
   cFeas(nFeas) = Result.h_k;
end

if f < fBest
   fBest   = f;
   ResBest = Result;
   IntBest = IntVars;
end

% MODIFICATION LOG:
%
% 070519  hkh  Written
% 110725  hkh  ConsIx changed to ConIx
