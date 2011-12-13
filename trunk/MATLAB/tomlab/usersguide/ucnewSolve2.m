x_0     = [-1.2;1];    % Starting values for the optimization.
x_L     = [-10;-10];   % Lower bounds for x.
x_U     = [2;2];       % Upper bounds for x.
Prob    = conAssign('ucnew_f', 'ucnew_g', [], ... % Setup Prob structure.
    [], x_L, x_U, 'ucNew',x_0);
Prob.P  = 18;                                      % Problem number.
Prob.Solver.Alg=1;                                 % Use quasi-Newton BFGS
Prob.user.uP = 100;                                % Set alpha parameter
Result  = tomRun('ucSolve',Prob,1);