function x=feasibilitycheck(c,x0,options)
% feasibilitycheck - Check for infeasible tomSym constraints.
%
% x=feasibilitycheck(c,x0) searches for a point satisfying the constraint
% list c by solving a relaxed problem starting with the guess x0.
%
% x=feasibilitycheck(c,x0,options) passes additional options to the solver.
%
% The constraint list can only contain equality and intequality
% constraints. Other types of constraints (integer, etc) are not
% allowed. (Constant constraints are allowed if they evaluate to "true".)
%
% If a point x is found that satisfies all constraints, then it is
% returned. Otherwise, a point which (locally) minimizes the sum of
% constraint violations is returned, and a list of constraints that
% were not satisfied is displayed.
%
% Often, the list of unsatisfied constraints will be realtively short, so
% running feasibilitycheck can significantly speed up the debugging
% process. (It is possible that the error is in one of the constraints not
% listed, which is infeasible in combination with one of the listed
% constraints. Still, the list of unsatisfied constraints might give a hint
% as to what is wrong.)
%
% options.prilev controls how much information is displayed. Default = 1.
%
% With nonlinear problems, it can happen that feasibilitycheck
% finds a feasible point even in cases where regular solvers fail, or vice
% versa. This is because feasibilitycheck minimizes the sum of errors while
% most solvers minimize the maximum error, so they take different paths
% through the search-space.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-11-20 by rutquist for TOMLAB release 7.7

if nargin<2
    x0 = {};
end

if nargin<3
    options = struct;
end

if ~isfield(options,'prilev') || isempty(options.prilev)
    options.prilev = 1;
end

PriLev = options.prilev;

[c,n] = flattencellarray(c,'c');

s = symbols(c,'struct');
sn = fieldnames(s);

if iscell(x0) || isa(x0,'tomSym')
    x0 = tom2struct(x0);
end

for i=1:length(sn)
    if ~isfield(x0,sn{i}) || ~isnumeric(x0.(sn{i}));
        if PriLev > 0
            disp(['Guess is missing for ' sn{i}]);
        end
        x0.(sn{i}) = zeros(size(s.(sn{i})));
    end
end

slack = cell(size(c));
txt = cell(size(c));
slacksum = 0;
isfeas = true;
for i=1:length(c)
    txt{i} = char(c{i});
    if tomCmp(c{i},'eq')
        slack{i} = tom([],size(c{i}));
        x0.(char(slack{i})) = abs(subs(operand(2,c{i}) - operand(1,c{i}), x0));
        c{i} = { -slack{i} <= operand(2,c{i}) - operand(1,c{i}) <= slack{i}, slack{i} >= 0 };
        slacksum = slacksum + sum(vec(slack{i}));
    elseif tomCmp(c{i},'le');
        slack{i} = tom([],size(c{i}));
        x0.(char(slack{i})) = max(0,subs(operand(1,c{i}) - operand(2,c{i}), x0));
        c{i} = { operand(1,c{i}) <= operand(2,c{i}) + slack{i}, slack{i} >= 0 };
        slacksum = slacksum + sum(vec(slack{i}));
    elseif tomCmp(c{i},'ge');
        slack{i} = tom([],size(c{i}));
        x0.(char(slack{i})) = max(0,subs(operand(2,c{i}) - operand(1,c{i}), x0));
        c{i} = { operand(1,c{i}) + slack{i} >= operand(2,c{i}), slack{i} >= 0 };
        slacksum = slacksum + sum(vec(slack{i}));
    elseif tomCmp(c{i},'positiveSemidefinite')
        M = operand(1,c{i});
        slack{i} = tom([],1,1);
        x0.(char(slack{i})) = max(0,-1.01*eigs(subs(0.5*(M+M'), x0), 1, 'sa'));
        c{i} = { positiveSemidefinite(M+slack{i}*eye(size(M))), slack{i} >= 0 };
        slacksum = slacksum + slack{i};
    elseif isa(c{i},'logical');
        if ~all(c{i}(:))
            error(['Constant constraint ' n{i} ' evaluated to false.']);
        end
    else
        if isa(c{i},tomSym)
            error(['Constraint ' n{i} '(' operator(c{i}) ') is not supported by feasibilitycheck']);
        else
            error(['Constraint ' n{i} 'is a ' class(c{i}) '. Expected a tomSym expression.']);
        end
    end
    if ~isempty(slack{i})
        if ~all(vec(x0.(char(slack{i}))<=1e-6))
            isfeas = false;
            if PriLev > 0
                disp(['Constraint ' n{i} ' not satisfied by initial guess.']);
            end
        end
    end
    if ~all(vec(isfinite(x0.(char(slack{i})))))
        error(['Initial guess resulted in inf/nan for constraint' n{i}]);
    end
end

if isfeas
    if PriLev > 0
        disp('Initial guess is a feasible point.')
    end
    solution = x0;
    
else
    if PriLev > 0
        disp('Solving the relaxed problem (where f_k=0 corresponds to a feasible point).')
    end
        
    if ~isfield(options,'scale')
        options.scale = 'auto';
    end
    
    o = options;
    o.feasibilitychecks = 0;
    [solution, result] = ezsolve(slacksum,c,x0,o);
    
    if PriLev > 0
        if result.ExitFlag ~= 0
            disp('The relaxed problem failed to give a reliable solution.');
        else
            if result.f_k < 5e-5
                disp('Found a feasible point.');
            else
                disp('The following constraints were not satisfied:');
                nes = 0;
                ncs = 0;
                for i=1:length(slack)
                    e = vec(solution.(char(slack{i})));
                    nes = nes+sum(e>=5e-6);
                    if ~all(e < 5e-6)
                        disp(['Constraint ' n{i} ' "' txt{i} '" ' ...
                            num2str(sum(e>=5e-6)) '/' ...
                            num2str(length(e)) ' elements, sum = ' ...
                            num2str(sum(max(e,0)))]);
                        ncs = ncs+1;
                    end
                end
                disp(['Total: ' num2str(ncs) ' unsatisfied constraints. (' num2str(nes) ' elements.)']);
            end
        end
    end        
end

if nargout > 0
    x = solution;
    for i=1:length(slack)
        x = rmfield(x,char(slack{i}));
    end
end

