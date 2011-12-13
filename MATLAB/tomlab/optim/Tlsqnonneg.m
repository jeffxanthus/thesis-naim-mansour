% Tlsqnonneg is the TOMLAB equivalent to LSQNONNEG in Optimization TB 3.x
%
% Tlsqnonneg is present for the users running Matlab 5.2 or older.
% In Matlab 5.3 lsqnonneg is available in \matlabr11\toolbox\matlab\matfun
%
% But this version is much faster, doing a fast MEX-file solution with Tnnls
%
% TOMLAB lsqnonneg solves:
%        linear least squares problems with nonnegativity constraints
%
%	min              ||C*x -d||^2 = f(x)
%    x
%                               x  >= 0
%
% x is a n-dimensional unknown parameter vector
%
% function [x, f_k, r_k, ExitFlag, Output, Lambda] = Tlsqnonneg(...
%           C, d, x_0, options)
%
% INPUT: ( 2 arguments always needed )
%
% C        Matrix C in C*x-d.
% d        Data vector in C*x-d
% x_0      Starting value for the x variables
% options  Replaces the default optimization parameters
%          Fields used: Display, TolX
%             Increase TolX is routine fails to find solution.
%
% OUTPUT:
%
% x        Optimal x parameters
% f_k      Optimal residual sum of squares, sum {r(x).^2} (Note! no 0.5)
% r_k      The optimal residual vector
% ExitFlag exit condition of lsqlin.
%      > 0 lsqlin converged to a solution X.
%        0 Reached the maximum number of iterations without convergence
%      < 0 Errors
%
% Output   Structure. Fields:
%   Output.iterations    Number of iterations
%   Output.algorithm     Type of algorithm used
%
% Lambda   Lagrange multipliers (dual parameters) at the solution
%          Lambda(i) <=0 when x(i) close to 0
%          Lambda(i) close to 0 when x(i) > 0

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sep 9, 1999.    Last modified Sept 10, 2009.

function [x, f_k, r_k, ExitFlag, Output, Lambda] = Tlsqnonneg(...
    C, d, x_0, options)

if nargin < 4, options = [];
    if nargin < 2
        error('lsqnonneg requires two input arguments C and d');
    end
end

TolXX = [];

if ~isempty(options)
    if isfield(options,'TolX')
        TolXX = options.TolX;
    end
end

[x,rNorm,mode,Iter,Lambda] = Tnnls(C, d, [],[],0,1,[],[],TolXX);

if mode == 0
    ExitFlag = 1;
    Output.message   = 'Convergence';
else
    ExitFlag = -1;
    Output.message   = 'Maximal number of iterations reached';
end

Output.algorithm   = 'TOMLAB Tnnls MEX-solver';
Output.iterations  = Iter;

% Use same sign convention as Matlab 5.3 lsqnonneg

r_k=d-C*x;
f_k=r_k'*r_k;


% MODIFICATIONS
%
% 001031 hkh Use wnnls MEX-file interface
% 011204 hkh Tiny changes
% 030116 hkh Change MEX call to Tnnls
% 030128 hkh Change name to Tlsqnonneg
% 040414 hkh Remove everything with MIDEVA
% 090910 hkh Iter now returned from Tnnls as 4th parameter
% 090910 hkh Output.message defined
