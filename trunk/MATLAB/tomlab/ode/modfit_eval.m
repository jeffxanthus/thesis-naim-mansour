% function z = modfit_eval(flag, x, y, t, c, Prob)
%
% TOMLAB gateway routine for various function evaluations asked for by MODFIT
%
% INPUT:
%
% flag    Flags which function evaluation is requested by MODFIT
%         Cases supported:
%         1       z = dy/dt   Evaluate ODE. modfit_eval calls the routine
%                             ode_f as z=ode_f(t, y, Prob)
%         2       z = y0      Evaluate ODE initial values
%         3       z = y(Eeq)  Evaluate ymodel corresponding to data series
%
% x       Current iterate x
% y       Current ODE function values, y(t,x)
% t       Current t value
% c       Current concentration (not supported yet)
% Prob    Problem structure in TOMLAB format
%
% OUTPUT:
%
% z       Depends on the input flag, described above

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function z = modfit_eval(flag, x, y, t, c, Prob)

switch flag
    case 1,
        % Evaluate ODE
        Prob.ODE.X = x;
        z = ode_f(t, y, Prob);
    case 2,
        % Evaluate ODE initial values
        Y0Idx = Prob.ODE.Y0Idx;
        if isempty(Y0Idx)
           Prob.ODE.X         = x;
        else
           % Split ODE initial value parameters from other unknown parameters
           Prob.ODE.Y0(Y0Idx) = x(1:length(Y0Idx));
           Prob.ODE.X         = x(length(Y0Idx)+1:end);
        end
        z = Prob.ODE.Y0(:);
    case 3,
        % Evaluate fitting criteria (ymodel)
        z = y(Prob.ODE.Eeq);
    case 4,
        % Evaluate constraints
        error('Constraint evaluation not defined');
    case 5,
        % Evaluate ODE gradient (w.r.t. y and x separately)
        error('ODE gradient (w.r.t. y and x) evaluation not defined');
    case 6,
        % Evaluate ODE initial gradient
        error('ODE initial gradient evaluation not defined');
    case 7,
        % Evaluate model function gradient
        error('Model function evaluation not defined');
    case 8,
        % Evaluate constraint gradient
        error('Constraint gradient evaluation not defined');
    case 9,
        % Evaluate ODE gradient (w.r.t. only y)
         error('ODE gradient evaluation not defined');
    otherwise,
        error('Invalid flag number');
end

% MODIFICATION LOG:
%
% 050421 bkh Written
% 050428 bkh Works ok with call chain modfitTL->modfitmex->modfit_eval
%            and modfitodeTL->modfitmex->modfit_eval
% 050705 med Help updated