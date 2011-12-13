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
%   'lpcon' - Linear objective with nonlinear constraints
%   'qpcon' - Quadratic objective with nonlinear constraints
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
% If the objective f is a vector or matrix, or has the form of an equation,
% then the norm of f (or the residual of the equation) is minimized.
% Allowed values for options.norm are:
% - 'L2' (default) minimizes 0.5*sum(vec(f)) (A least-squares problem.) 
% - 'L1' minimizes sum(abs(f))
% - 'LInf' minimizes max(abs(f))
% The objective will be transformed, introducing extra unknowns in order to
% avoid the nonlinearities that the functions "max" and "abs" represent.
%
% sym2prob(type,f,c,x0) supplies an initial guess for one or more of the
% unknowns. The guess x0 can be a struct, where each field is named
% as a symbol and contains a numeric array of the correct size. x0 can also
% be a cell array of equations compatible with tom2struct.
%
% sym2prob(type,f,c,x0,OPTIONS) where OPTIONS is a structure sets options.
% The options structure can have the following fields.
%
%   OPTIONS.name    - The name of the problem
%   OPTIONS.solver  - The solver that will be used (alows sym2prob to guess
%                     some settings.)
%   OPTIONS.use_g   - (boolean) true = generate symbolic gradient
%   OPTIONS.use_H   - (boolean) true = generate symbolic Hessian
%   OPTIONS.use_dc  - (boolean) true = symbolic constraint derivatives
%   OPTIONS.use_d2c - (boolean) true = symbolic second constr. derivs.
%   OPTIONS.scale   - 'auto' = autoscale the problem, 
%                     'man'  = generate code that works with scaleProb
%                     ''     = do not use scaling (default)
%   OPTIONS.checkguess - If 'true' print a warning if guess is out of
%                        bounds.
%   OPTIONS.Prob.*  - Any options that can be sent to the solver via the
%                     "Prob" struct can be used here.
%
% (Note that the use_X flags only determine whether thevsymbolic
% derivatives will be generated. Whether the solver will actually use them
% depends on the solver.)
%
% Any flags that are set in the substructure OPTIONS.Prob will be copied
% directly into the Prob struct that results.
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
