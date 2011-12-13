% mcpAssign is a direct way of setting up a Mixed Complementarity Problem (MCP)
% in the TOMLAB format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = mcpAssign(....)
%
% It is then possible to solve the MCP using the TOMLAB KNITRO solver
% knitroTL with the call:  Result = tomRun('knitro', Prob, PriLev);
%
% See the file tomlab\examples\mcptest1.m for an example
%
% Syntax of mcpAssign:
%
% INPUT (Call with at least seven parameters)
%
% function Prob = mcpAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
%                           MPEC, fLowBnd, ...
%                           A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
%                           x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least seven parameters)
%
% f           Name of the function that computes the function value f(x)
% g           Name of the function that computes the n x 1 gradient vector
% H           Name of the function that computes the n x n Hessian matrix
% HessPattern n x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the Hessian and ones indicate values that might
%             be non-zero. If empty indicates estimation of all elements
%             HessPattern is used when estimating the Hessian numerically.
%             Estimated before solve, if Prob.LargeScale==1, HessPattern==[]
%
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as a nx1  Inf vector.
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
% Note:       The number n of the unknown variables x are taken as
%             max(length(x_L),length(x_U),length(x_0))
%             You must specify at least one of these with correct length,
%             then the others are given default values
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
%
% MPEC        Each row of mpec is one pair. The example below says
%             x(3) _|_ c(2) should be complementary.
%             Exactly two nonzeros per row is allowed.
%
%               mpec = [ ...
%                  3,0,2,0,0,0; ...
%                  4,0,3,0,0,0; ...
%                  5,0,4,0,0,0; ...
%                  ];
%
%             mpec = [   var1,var2 , lin1,lin2 , non1,non2  ; ... ];
%
% fLowBnd     A lower bound on the function value at optimum. Default -1E300
%             A good estimate is not critical. Use [] if not known at all.
%
% L I N E A R   C O N S T R A I N T S
% A           mA x n matrix A, linear constraints b_L <= A*x <= b_U. Dense or sparse
% b_L         Lower bound vector in linear constraints b_L <= A*x <= b_U.
% b_U         Upper bound vector in linear constraints b_L <= A*x <= b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints
% dc          Name of function that computes the constraint Jacobian mN x n
% d2c         Name of function that computes the second part of the
%             Lagrangian function (only needed for some solvers)
%             See the help gateway routine nlp_d2c for an explanation of d2c
%
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
%
% c_L         Lower bound vector in nonlinear constraints c_L <= c(x) <= c_U.
% c_U         Upper bound vector in nonlinear constraints c_L <= c(x) <= c_U.
%
% A D D I T I O N A L   P A R A M E T E R S
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written Apr 14, 2004.    Last modified Dec 12, 2006.

function Prob = mcpAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
                          MPEC, fLowBnd, ...
                          A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ... 
                          x_min, x_max, f_opt, x_opt)

if nargin < 23
   x_opt=[];
   if nargin < 22
      f_opt=[];
      if nargin < 21
         x_max=[];
         if nargin < 20
            x_min=[];
            if nargin < 19
               c_U=[];
               if nargin < 18
                  c_L=[];
                  if nargin < 17
                     ConsPattern=[];
                     if nargin < 16
                        d2c=[];
                        if nargin < 15
                           dc=[];
                           if nargin < 14
                              c=[];
                              if nargin < 13
                                 b_U=[];
                                 if nargin < 12
                                    b_L=[];
                                    if nargin < 11
                                       A=[];
                                       if nargin < 10
                                          fLowBnd=[];
                                          if nargin < 9
                                             MPEC=[];
                                             if nargin < 8
                                                x_0=[];
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end
                      
if (~issparse(MPEC))
    MPEC = sparse(MPEC);    % Make sure MPEC is sparse
end

if size(MPEC,2) ~= 6
    error('Input parameter MPEC should have 6 columns');
end

if isempty(f)
    f = 'mcpdummy_f';
    g = 'mcpdummy_g';
    H = 'mcpdummy_H';
    N = max([length(x_L),length(x_U),length(x_0)]);
    HessPattern = sparse(N,N);
end

Prob = conAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
                          [], fLowBnd, ...
                          A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ... 
                          x_min, x_max, f_opt, x_opt);

Prob.probType = checkType('mcp');

Prob = BuildMPEC(Prob,MPEC);

% MODIFICATION LOG
%
% 060524  med  Written, based on conAssign
% 060705  med  Updated help
% 060818  hkh  Updated help regarding fLowBnd
% 061212  med  Help updated (mpec has 6 columns)