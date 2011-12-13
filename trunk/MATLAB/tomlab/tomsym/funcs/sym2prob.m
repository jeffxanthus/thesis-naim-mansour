% sym2prob - Compile symbolic function/constraints into a Prob struct.
%
% Prob = sym2prob(type,f,c) creates a Prob structure, suitable for tomRun,
% with the objective of minimizing f, subject to c, with respect to all
% symbols that they contain.
%
% The problem type can be (among others):
%
%   'lp'    - Linear programming
%   'qp'    - Quadratic programming
%   'con'   - Nonlinear programming
%   'qpcon' - Quadratic problem with nonlinear constraints
%   'minlp' - Mixed integer nonlinear programming
%   'sdp'   - Linear semidefinite programming
%   'bmi'   - Bilinear semidefinite programming
%
% Prob = sym2prob(f,c) calls tomDiagnose to attempt to guess the problem
% type, and prints a warning.
%
% The objective f must be a tomSym symbolic object, while the constraints
% list c should be a cell array.
%
% If the objective f is a vector or matrix, then a least-square problem is
% solved, minimizing 0.5*sum(vec(f)) (Half the sum-of-squares of the 
% elements of f). If f is an equality on the form lhs == rhs, then the
% sum-of-squares of the difference (rhs-lhs) is minimized.
%
% sym2prob(type,f,c,x0) supplies an initial guess for one or more of the
% unknowns. The guess x0 should be a struct, where each field is named
% as a symbol and contains a numeric array of the correct size.
%
% sym2prob(type,f,c,x0,OPTIONS) where OPTIONS is a structure sets options.
% The options structure can have the following fields.
%
%   OPTIONS.name    - The name of the problem
%   OPTIONS.use_d2c - (boolean) true = compute symbolic d2c
%   OPTIONS.use_H   - (boolean) true = compute symbolic H
%   OPTIONS.prilev  - 0 = quiet, 1 = display some information
%   OPTIONS.Prob    - Further options that will be copied to the Prob
%                     structure.
%
% Linear and box constraints will be automatically detected, if they are
% formulated using simple addition and multiplication. For example,
%   -3*(x+2) <= 4+x
% is automatically converted to
%          x >= -2.5
%
% For nonlinear problems, sym2prob will generate temporary m-files
% representing the nonlinear functions and their derivatives (These files
% are usually small and harmless, and many modern operating systems
% automatically clean up old temporary files). In order to remove these
% files when they are no longer needed, it is recommended to run
% tomCleanup(Prob) after the problem has been solved.
