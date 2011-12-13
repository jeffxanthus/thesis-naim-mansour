qpExample;

Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, 'qpExample');

Result = tomRun('qpSolve', Prob, 1);