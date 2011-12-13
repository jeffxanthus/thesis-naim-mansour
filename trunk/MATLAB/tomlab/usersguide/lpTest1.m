lpExample;

Prob   = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'lpExample');
Result = tomRun('lpSimplex', Prob, 1);