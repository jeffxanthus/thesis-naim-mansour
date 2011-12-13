function [solution, result]=ezsolve(f,c,x0,options)
% ezsolve - Solve a tomSym optimization problem.
%
% [solution, result] = ezsolve(f,c) returns the solution to the minimization
% problem that is defined by the objective function f and the constraints c.
% The result structure from tomRun is provided in a second output argument.
%
% Ezsolve can also be used to find least-square solutions to equations. If
% options.norm is set, then f can be a set of equations. (If there
% are equations both in f and c, then the ones in c are considered as
% strict, and will be solved to tolerances, while the ones in f are solved,
% for example, in a least-square sense. If f is a vector, then the norm of
% f is minimized. 
% The norm that is used is defined by options.norm. Allowed values are:
% - 'L2' (default) minimizes 0.5*sum(vec(f)) (A least-squares problem.) 
% - 'L1' minimizes sum(abs(f))
% - 'LInf' minimizes max(abs(f))
%
% Ezsolve is meant to be as simple to use as possible. It automatically
% determines the problem type and finds a suitable solver.
%
% The returned solution is a struct, where the fields represent the unknown
% variables. This struct can be used by subs to convert a tomSym to a
% numeric value.
%
% s = ezsolve(f,c,x0) uses the initial guess x0. The input argument x0 can
% be a struct containing fields names as the unknown symbols, for example a
% previously returned solution. Alternatively, x0 can be a cell array of
% simple tomSym equation that can be converted to a struct using tom2struct.
%
% s = ezsolve(f,c,x0,name) sets the problem name.
%
% s = ezsolve(f,c,x0,OPTIONS) where OPTIONS is a structure sets solver
% options. The options structure can have the following fields.
%
%   OPTIONS.name    - The name of the problem
%   OPTIONS.type    - The problem type, e.g. 'lp', 'qp', 'con', ...
%   OPTIONS.solver  - The solver to use, e.g. 'snopt', 'knitro', ...
%   OPTIONS.prilev  - The ezsolve and tomRun print level (default 1)
%   OPTIONS.use_d2c - (boolean) true = compute symbolic d2c
%   OPTIONS.use_H   - (boolean) true = compute symbolic H
%   OPTIONS.scale   - 'auto' = autoscale the problem, 
%                     'man'  = generate code that works with scaleProb
%                     ''     = do not use scaling
%
%   OPTIONS.noCleanup   - If true, then temporary files are not deleated,
%                         which is useful for debugging or profiling.
%                         (Default = true when profiler is active).
%   OPTIONS.docollocate - If true, call docollocate for delayed
%                         collocation.
%
%   OPTIONS.proptIterations - See "help proptIterate"
%
%   OPTIONS.feasibilitychecks - The number of times to pre-sovle for a
%   feasible solution, before starting optimization. Default = 0. Setting
%   this to 1 may display some helpful hints in diagnosing conflicting
%   constraints.
%
% Ezsolve works as a wrapper to tomRun. Any options that can be given to
% tomRun in the Prob stucture can also be used with ezsolve. Just set the
% corresponding field in the OPTIONS.Prob substructure.
%
% See also: tomDiagnose, sym2prob, tomRun

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-11-30 by rutquist for TOMLAB release 7.7

global lastTomSymSolution

if nargin<2
    c = {};
end
if nargin<3
    x0 = struct;
end
if nargin<4
    options = struct;
end

if ischar(options)
    n = options;
    options = struct;
    options.name = n;
end
if ~isfield(options,'name')
    options.name = '';
end
if ~isfield(options,'scale')
    options.scale = '';
end

if numel(f)==1
    [f,c] = rewriteV(f,c);
end

if ~isfield(options,'prilev') || isempty(options.prilev)
    options.prilev = 1;
end

PriLev = options.prilev;

if isfield(options,'proptIterations') && ~isempty(options.proptIterations)
    [solution, result]=proptIterate(f,c,x0,options);
    return
end

% Do delayed collocation, if applicable
global proptActivePhase
if isfield(options,'phase')
    phase = options.phase;
elseif ~isempty(proptActivePhase)
    phase = proptActivePhase;
else
    phase = [];
end

if ~isempty(phase) && ~(isfield(options,'docollocate') && ~options.docollocate)
    if isempty(phase.cpoints)
        proptIterations = 16;
        phase = tomPhase(phase.label, phase.t, phase.tstart, phase.tdelta, proptIterations);
    end
    [c,s] = docollocate(phase,c);
    [f,s] = docollocate(phase,f,s);
    x0 = docollocate(phase,x0,s);
    options.docollocate = false;
end

if iscell(x0) || isa(x0,'tomSym')
    x0 = tom2struct(x0);
end

% Extract constraints from any subjectTo calls
[f,c,x0] = extractConstraints(f,c,x0); 

% Find problem type
if ~isfield(options,'type');
    [options.type, options.complementary] = tomDiagnose(f,c,options);
    if PriLev > 0
        disp(['Problem type appears to be: ' options.type]);
    end
end

if ~isfield(options,'solver')
    options.solver = GetSolver(options.type,true);
    if isfield(options,'complementary') && options.complementary && ...
            checkLicense('knitro')
        options.solver = 'knitro';
    end
end

% There seems to be a problem with checkLicence and certain solvers.
% Skipping this check for the time being.
%if ~checkLicense(options.solver)
%    warning('tomSym:noLicense',...
%        ['Could not obtain a license for solver: ' options.solver '.']);
%    options.solver = GetSolver(options.type,true);
%    if PriLev > 0
%        disp(['Instead, using solver: ' options.solver]);
%    end
%end

if ~isfield(options,'use_d2c')
    if any(strcmpi(options.solver,{'snopt','snopt7','oqnlp','minos','npsol'}))
        options.use_d2c = false;
    end
end

if ~isfield(options,'use_H')
    if any(strcmpi(options.solver,{'snopt','snopt7','oqnlp','minos','npsol'}))
        options.use_H   = false;
    end
end

if isfield(options,'adObj') && options.adObj && ...
    isfield(options,'adCons') && options.adCons
    % The syntax for ezsolve is to set options.derivatives, but if someone
    % set options.adObj and options.adCons, we know what they meant.
    options.derivatives = 'automatic';
end

if isfield(options,'feasibilitychecks')
    for fcki = 1:options.feasibilitychecks
        x0 = feasibilitycheck(c,x0,options);
    end
end

% Note: proptDEBUG had been removed. Do "dbstop if caught error" instead.
try
    Prob = sym2prob(options.type,f,c,x0,options);
    Prob = transferOptions(Prob,options);
    if isfield(options,'derivatives') && ...
            (strcmpi(options.derivatives, 'automatic') ...
        || strcmpi(options.derivatives, 'auto') ...
        || strcmpi(options.derivatives, 'MAD'))
        madinitglobals;
        Prob.ADObj = 1;
        Prob.ADCons = 1;
    end
    if PriLev > 0
        disp('Starting numeric solver');
    end
    result = tomRun(options.solver,Prob,PriLev);
    solution = getSolution(result);
    if isfield(Prob.tomSym,'fScale') && Prob.tomSym.fScale~=1 && PriLev > 0
        disp(['Objective value, corrected for scaling: ' ...
            num2str(result.f_k./Prob.tomSym.fScale)]);
    end
catch le
    if exist('Prob','var')
        if ~isfield(options,'nocleanup') || ~options.nocleanup
            tomCleanup(Prob);
        end
    end
    rethrow(le);
end

if result.ExitFlag ~= 0
    warning('ezsolve:ExitFlag',['Solver returned ExitFlag = ' num2str(result.ExitFlag)])
    disp('The returned solution may be incorrect.')
elseif isfield(Prob,'MIP') && isfield(Prob.MIP,'IntVars') && any(Prob.MIP.IntVars) && ...
        any(abs(result.x_k(Prob.MIP.IntVars)-round(result.x_k(Prob.MIP.IntVars)))>1e-4)
    warning('ezsolve:integer','Solver returned non-integer solution to a mixed-integer problem.');
    disp('Make sure you have selected a solver that can handle mixed-integer problems.');
end

lastTomSymSolution = result; % Keep result for ezplot, etc.

if ~isfield(options,'noCleanup')
   pstatus = profile('status');
   options.noCleanup = strcmp(pstatus.ProfilerStatus,'on');
end
    
if ~options.noCleanup
    tomCleanup(Prob); % Delete temporary files.
end

function Prob = transferOptions(Prob,options)
% Set some options based on solver choice
if strcmp(options.solver,'snopt') && any(strcmp(options.type,{'qp','qpcon'}))
    Prob.QP.useQPobj = 1;
end
% Transfer options into the Prob struct
% Note - This is obsolete. Use options.Prob instead.
copylist = {'xInit','PriLevOpt'};
for i=1:length(copylist)
    if isfield(options,copylist{i})
        Prob.(copylist{i}) = options.(copylist{i});
    end
end
