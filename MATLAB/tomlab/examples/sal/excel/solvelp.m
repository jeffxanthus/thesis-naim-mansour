function [x, Inform] = solvelp(c, A, x_L, x_U, b_L, b_U)

[x, slack, v, rc, f_k, ninf, sinf, Inform] = cplex(c, A, x_L, x_U, ...
                                                  b_L, b_U);