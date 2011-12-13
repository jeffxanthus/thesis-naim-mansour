% minlpQGsym - TomSym version of minlpQG

Name='minlp1Demo - Kocis/Grossman.';

toms 2x1          x
toms 3x1 integer  y

objective = [2 3 1.5 2 -0.5]*[x;y];

constraints = {
    x(1) >= 0
    x(2) >= 1e-8
    x <= 1e8
    0 <= y <=1
    [1 0 1 0 0]*[x;y] <= 1.6
    1.333*x(2) + y(2) <= 3
    [-1 -1 1]*y <= 0
    x(1)^2+y(1) == 1.25
    sqrt(x(2)^3)+1.5*y(2) == 3};

guess = struct('x',ones(size(x)),'y',ones(size(y)));
options =  struct;
options.name = Name;
Prob = sym2prob('minlp',objective,constraints,guess,options);

Prob.DUNDEE.optPar(20) = 1;

% Get default TOMLAB solver for your current license, for "minlp" problems
% Solver = GetSolver('minlp');

% Call driver routine tomRun, 3rd argument > 0 implies call to PrintResult

Result  = tomRun('minlpBB',Prob,2);

tomCleanup(Prob);
