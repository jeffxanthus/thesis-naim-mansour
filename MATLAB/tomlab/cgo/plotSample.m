%function fig = plotSample(Result,Action)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified April 10, 2008.

function fig = plotSample(Result,Action)

d = Result.Prob.N;

if d ~= 2
    return
end

if nargin == 2
    if Action == 0
        pause off
    end
end

% Info from Result
%x_k         = Result.x_k;
%f_k         = Result.f_k;

Digits = [];
if isfield('Result','Digits')
    if ~isempty(Result.DIGIT)
        Digits  = unique(Result.DIGIT.FuncEv);
    end
end

% Info from Result.Prob
%prob_type   = Result.Prob.probFile;
%prob_number = Result.Prob.P;
x_opt       = Result.Prob.x_opt;
f_opt       = Result.Prob.f_opt;
%x_L         = Result.Prob.x_L;
%x_U         = Result.Prob.x_U;
%Func_f      = Result.Prob.FUNCS.f;

% Info from Result.CGO
X           = Result.CGO.WarmStartInfo.O;
F           = Result.CGO.WarmStartInfo.F;
nInit       = Result.CGO.WarmStartInfo.nInit;
%fMinIdx     = Result.CGO.WarmStartInfo.fMinIdx;
%Scale       = Result.Prob.CGO;

fig = plotProblem(Result.Prob);

plot3(x_opt(1),x_opt(2),f_opt,'r*','MarkerSize',18)
plot3(X(1,1:nInit),X(2,1:nInit),F(1:nInit),'b.','MarkerSize',15)

for k = nInit+1:length(F)
    pause
    if ismember(k,Digits)
        plot3(X(1,k),X(2,k),F(k),'g.','MarkerSize',15)
    else
        plot3(X(1,k),X(2,k),F(k),'k.','MarkerSize',15)
    end
end

if nargin == 2
    if Action == 0
        pause on
    end
end


% MODIFICATION LOG:
%
% 080410 hkh Written 
