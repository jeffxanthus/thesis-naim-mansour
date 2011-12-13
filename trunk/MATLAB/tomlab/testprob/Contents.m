% TOMLAB - Test Problems
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% tomlab\testprob
%
% Predefined problems in TOMLAB
%
% contents     This file
%
% bmi_prob   Defines predefined BMI problems
%
% chs_prob   Hoch-Schittkowski constrained test problem set
% chs_c      Constraint residuals c(x) for HS problems
% chs_dc     Derivative of constraint residuals for HS problems
% chs_f      The objective function f(x) for HS problems
% chs_g      The gradient g(x) for HS problems
% chs_H      The Hessian H(x) of f(x) for HS problems
%
% cls_prob   Defines predefined Constrained NLLS problems
% cls_c      Nonlinear constraints c(x) for constrained LS problem
% cls_dc     Derivative of nonlinear constraints for CLS problem
% cls_J      Jacobian matrix dr_i/dx_j, i=1,...,m, j=1,...,n
% cls_r      Residuals r_i(x), i=1,...,m. x = (x_1,...,x_n)^T
%
% con_prob   Defines predefined constrained problems
% con_c      Constraint residuals c(x) for constrained problem
% con_d2c    2nd part of 2nd derivative of Lagrangian function
% con_dc     Derivative of constraint residuals for constrained problem
% con_f      The objective function f(x) for constrained problem
% con_g      The gradient g(x) for constrained problem
% con_H      The Hessian H(x) of f(x) for constrained problem
%
% exp_d2r    Compute d2r in Hessian H(x) = J^T * J + d2r (2nd der part)
% exp_J      Compute the Jacobian matrix dr_i/dx_j, i=1,...,m, j=1,...,n
% exp_prob   Defines problem data for EXPFIT. Name and series (t,y)
%            Data from several different empirical test series.
% exp_r      Compute the residuals r_i(x), i=1,...,m. x = (x_1,...,x_n)^T
%
% Helax AB two-term exponential test problems in five mat files
% xp020.mat, xp031.mat, xp200.mat, xp400.mat, yp050.mat, yp250.mat
% Used by exp_prob.m.
% funexpm2   Used as subfunction in initial value alg in exp.fitting.
% funm2exp   Used as subfunction in initial value alg in exp.fitting.
%
% geno_c     Constraint vector c(x) for GENO problems
% geno_dc    Derivative of constraint Jacobian for GENO problems
% geno_d2c   2nd part of 2nd derivative of Lagrangian function
% geno_f     The objective function f(x) for GENO problems
% geno_g     The gradient g(x) for GENO problems
% geno_H     The Hessian H(x) of f(x) for GENO problems
% geno_prob  Defines predefined problems from the GENO manual
%
% glb_f      Compute the objective function f(x) for global optimization
% glb_prob   Defines predefined box-bounded global optimization problems
%
% glc_c      Compute the constraint functions c(x).
% glc_f      Compute the objective function value f(x).
% glc_prob   Defines predefined constrained global optimization problems
%
% goals_c    Compute the constraint functions c(x).
% goals_dc   Compute the constraint Jacobian dc(x).
% goals_J    Compute the residuls Jacobian r(x).
% goals_prob Defines predefined goal attainment problems.
% goals_r    Compute the residuls r(x).
%
% lgo1_prob  1-dimensional global optimization problems.
% lgo1_f     Compute the objective function value f(x).
%
% lgo2_prob  higher-dimensional global optimization problems.
% lgo2_f     Compute the objective function value f(x).
%
% lls_prob   Test problems for linear least squares.
%
% lp_prob    Test problems for linear programming (LP)
% lp_probmat Mat-file for some of the LP problems
%
% ls_prob    Defines predefined NLLS problems
% ls_d2r     If possible, compute the second derivatives of the residual
% ls_J       Jacobian matrix J(x) = dr_i/dx_j, i=1,...,m, j=1,...,n
% ls_r       Residuals r_i(x), i=1,...,m. x = (x_1,...,x_n)^T
%
% mco_c      Compute the constraint functions c(x).
% mco_dc     Compute the constraint Jacobian dc(x).
% mco_J      Compute the residuls Jacobian r(x).
% mco_prob   Defines predefined goal attainment problems.
% mco_r      Compute the residuls r(x).
%
% mgh_prob   Predefined NLLS problems from More, Garbow, Hillstrom (MGH)
% mgh_J      Jacobian matrix J(x) = dr_i/dx_j, i=1,...,m, j=1,...,n (MGH)
% mgh_r      Residuals r_i(x), i=1,...,m. x = (x_1,...,x_n)^T (MGH)
%
% minlp_prob Defines predefined MINLP problems
% minlp_c    Constraint residuals c(x) for constrained problem
% minlp_d2c  2nd part of 2nd derivative of Lagrangian function
% minlp_dc   Derivative of constraint residuals for constrained problem
% minlp_f    The objective function f(x) for constrained problem
% minlp_g    The gradient g(x) for constrained problem
% minlp_H    The Hessian H(x) of f(x) for constrained problem
%
% mip_prob   Defines predefined MILP problems
% mip_probmat Mat-file for some of the MILP problems
%
% miqp_prob  Defines predefined MIQP problems
% miqp_probmat Mat-file for some of the MIQP problems
%
% miqq_prob  Defines predefined MIQQ problems
%
% mpex_prob  Defines predefined Mixed Complementarity example problems
%
% qp_prob    Defines predefined QP problems
%
% sdp_prob   Defines predefined SDP problems
%
% tsp_prob   Defines TSP (Travelling Salesman Problems) problems
% tsprun     Run the TSP problems with salesman method
%
% uc_prob    Defines predefined unconstrained problems
% uc_f       Compute the objective function f(x)
% uc_g       Compute the gradient g(x)
% uc_H       Compute the Hessian H(x)
%
% uhs_prob   Hoch-Schittkowski unconstrained test problem set
% uhs_f      The objective function f(x) for unconstrained HS problems
% uhs_g      The gradient g(x) for unconstrained HS problems
% uhs_H      The Hessian H(x) of f(x) for unconstrained HS problems
