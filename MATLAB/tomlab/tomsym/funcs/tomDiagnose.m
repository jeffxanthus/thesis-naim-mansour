% tomDiagnose - Determine the type for of tomSym optimization problem.
%
% type = tomDiagnose(f,c) returns a text string, defining the problem type
% that is represented by the objective function f and the constraints c.
%
% Some possible return values:
%
%   'lp'    - Linear programming
%   'qp'    - Quadratic programming
%   'cls'   - Least squares with linear constraints
%   'nls'   - Nonlinear least squares
%   'con'   - Nonlinear programming
%   'qpcon' - Quadratic problem with nonlinear constraints
%   'minlp' - Mixed integer nonlinear programming
%   'sdp'   - Linear semidefinite programming
%   'bmi'   - Bilinear semidefinite programming
%
% The objective f must be a tomSym symbolic object, while the constraints
% list c should be a cell array of tomSym objects.
