toms t
p = tomPhase('p', t, 0, 5, 60);
setPhase(p);

tomStates x1 x2 x3
tomControls u

% Initial guess
x0 = {icollocate({x1 == 0; x2 == 1; x3 == 0})
    collocate(u == -0.01)};

% Box constraints
cbox = {-10  <= icollocate(x1) <= 10
    -10  <= icollocate(x2) <= 10
    -10  <= icollocate(x3) <= 10
    -0.3 <= collocate(u)   <= 1};

% Boundary constraints
cbnd = initial({x1 == 0; x2 == 1; x3 == 0});

% ODEs and path constraints
ceq = collocate({dot(x1) == (1-x2.^2).*x1-x2+u
    dot(x2) == x1; dot(x3) == x1.^2+x2.^2+u.^2});

% Objective
objective = final(x3);

% Solve the problem
options = struct;
options.name = 'Van Der Pol';
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);