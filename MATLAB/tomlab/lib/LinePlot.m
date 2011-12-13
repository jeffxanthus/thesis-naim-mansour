% LinePlot.m
%
% function LinePlot(f, x, p, f_0, g_0, LineParam, alphaMax, ...
%          pType, alpha, varargin);
%
% See LineSearch for variable description
%
% Plot along line (search direction)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written June 11, 1995.  Last modified Sep 24, 2000.

function LinePlot(f, x, p, f_0, g_0, LineParam, alphaMax, ...
    pType, alpha, varargin)

global  n_f

fType = CheckFunc( f, pType);

rho     = LineParam.rho;
fLowBnd = LineParam.fLowBnd;

fp_0 = g_0'*p;			% derivative in x in search direction,
% fp(x) =  g(x)' * p

if fp_0 == 0
    my = alphaMax; % maximal value for alpha
else
    my = min(alphaMax,(fLowBnd - f_0) / (rho*fp_0)); % maximal value for alpha
end
if my < 0
    my=alphaMax;
end

alpha_max = min([my 3 ]);

fprintf('alphaMax = %20.5f\n',alpha_max);
m=300;

step=alpha_max/m;
X=0:step:alpha_max;
Y=zeros(m+1,1);
Y(1)=f_0;
n_f_0=n_f;

for k=1:m
    alpha_k=X(k+1);
    xEval=x+alpha_k*p;
    f_alpha_k = evalFunc(f, xEval, fType, varargin{:});
    Y(k+1)=f_alpha_k;
end
ymin=min(Y);
log=0;
if log
    if ymin <= 0
        add=1;
        Y=log10(Y-ymin+1E-15);
    else
        add=1;
        Y=log10(Y-ymin+1E-15);
    end
else
end
plot(X,Y);
hold on;
title('One dimensional line search problem');
xlabel('alpha');
j=min(m+1,1+round(alpha/step));
plot(X(j),Y(j),'x');
if log
    if add==0
        ylabel('log10(f(x))');
    else
        ylabel(['log10(f(x)+min(f)=' num2str(ymin) ' +1E-15)']);
    end
else
    ylabel('f(x)');
end
hold off;

n_f=n_f_0; % Reset the count of function evaluations. Using m func.eval/plot.
pause;

% ============================================================
function [fType] = CheckFunc( f, pType)
% ============================================================
if strcmp(f,'nlp_f')
    fType=1;
elseif strcmp(f,'nlp_r')
    fType=2;
elseif strcmp(f,'con_fm')
    fType=3;
else
    fType=4+pType;
end

% ============================================================
function f_k = evalFunc( f, x, fType, varargin)
% ============================================================

nargin;
switch fType
    case 1
        f_k = nlp_f(x,varargin{:});
    case 2
        r_k = nlp_r(x, varargin{:} );     % Function value
        f_k = 0.5 * (r_k' * r_k);
    case 3
        f_k = con_fm(x, [], varargin{:});
    case 4
        f_k = feval(f,x,varargin{:});
    case {5,6}
        r_k = feval(f,x,varargin{:});
        f_k = 0.5 * (r_k' * r_k);
    case {7,8}
        f_k = feval(f, x, [], varargin{:});
end

% MODIFICATION LOG:
%
% 981005  hkh  Revised the code for the new NLPLIB structure approach
% 981017  hkh  Remove g as input argument
% 981026  hkh  Get fLow from optParam.fLow
% 000928  hkh  Revise the parameter handling