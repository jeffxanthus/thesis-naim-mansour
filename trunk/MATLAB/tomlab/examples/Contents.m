% TOMLAB - Test Examples
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% tomlab\examples
%
% Test examples on how to use Tomlab, not using the GUI or Menu system
%
% Contents        This file
%
% CGOrun.m        Runs costly global optimization solvers for predefined probs
% -----------------------------------------------------------------------------
% conDemo.m       Solution of constrained optimization (nonlinear programming)
% con1_f.m        The function computation for the conDemo functions
% con1_g.m        The gradient for the conDemo functions
% con1_H.m        The Hessian for the conDemo functions
% con1_c.m        The constraints for the conDemo functions
% con1_dc.m       The constraint gradients for the conDemo functions
% con1_d2c.m      The 2nd order constraint information for the conDemo functions
% -----------------------------------------------------------------------------
% expDemo.m       Solution of exponential fitting problems
% -----------------------------------------------------------------------------
% glbDemo.m       Solution of global optimization (box-bounded)
% glb1_f.m
% glb4_f.m
% glb5_f.m
% -----------------------------------------------------------------------------
% glcDemo.m       Solution of global optimization (box-bounded,integer,constr.)
% glc4_c.m
% glc4_f.m
% glc5_c.m
% glc5_f.m
% -----------------------------------------------------------------------------
% L1Demo.m        Solution of L1 data fitting problems
% L1ex_r.m        The residual routine for L1Demo
% L1ex_J.m        The Jacobian routine for L1Demo
% -----------------------------------------------------------------------------
% llsDemo         Solution of linear least squares
% -----------------------------------------------------------------------------
% lpDemo.m        Solution of linear programming
% -----------------------------------------------------------------------------
% lsDemo.m        Solution of nonlinear least squares
% ls1_r.m         The residual routine for the lsDemo functions
% ls1_J.m         The Jacobian routine for the lsDemo functions
% -----------------------------------------------------------------------------
% mipDemo.m       Solution of mixed-integer programms
% -----------------------------------------------------------------------------
% minimaxDemo.m   Solution of minimax problems
% mima_r.m        Residual in the minimax test problem
% mima_J.m        Constraint Jacobian in the minimax test problem
% -----------------------------------------------------------------------------
% minlpDemo       Solution of mixed-integer nonlinear problems
% -----------------------------------------------------------------------------
% mpecDemo        Solution of mixed-complimentary problems
% mpd_f           Objective for MPEC
% mpd_g           Gradient for MPEC
% mpd_H           Hessian for MPEC
% mpd_c           Constraints for MPEC
% mpd_dc          Jacobian for MPEC
% mpd_d2c         Constraint Hessian for MPEC
% -----------------------------------------------------------------------------
%                 Test examples for interior point solver pdco
% -----------------------------------------------------------------------------
%                 Run as pdctest?(m,n), where ? is empty,2,3,4
%                 m = number of constraints, n = number of variables
%                 e.g pdcotest(50,100), bigger problem pdcotest(500,1000)
% pdcotest.m      Test on entropy problem
% pdcotest2.m     Test on entropy problem , implicit routine pdcoA.m used
% pdcoA.m         Implicit routine computing A*x,A'*x, using global A
% pdcotest3.m     Test LP problem
% pdcotest4.m     Test LP problem in two steps, faster than pdcotest3
% entropy         An example objective function
% linobj          An example linear objective function
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
%                 Test example for interior point solver pdsco
% -----------------------------------------------------------------------------
%                 Run as [x,y,z,inform]=entnet('entnet.dat')
% entnet.m        Test on network problem with entropy objective
% entnet.dat      List of i j pairs (easily changed)
% -----------------------------------------------------------------------------
% qpDemo.m        Solution of quadratic programming
% -----------------------------------------------------------------------------
% sdpDemo.m       Solution to semi-definite programming
% -----------------------------------------------------------------------------
% slsDemo.m       Solution of Sparse Least Squares (sls) problems
% -----------------------------------------------------------------------------
% ucDemo.m        Solution of unconstrained optimization
% uc1_f.m         The function computation for the ucDemo functions
% uc1_g.m         The gradient for the ucDemo functions
% uc1_H.m         The Hessian for the ucDemo functions
% uc3_f.m         A more advanced type of function computation for ucDemo
% uc4_f.m         The function computation for the ucDemo functions
% uc4_g.m         The gradient for the ucDemo functions
% uc4_H.m         The Hessian for the ucDemo functions
%
% -----------------------------------------------------------------------------
% psfDemo         Illustrates the use of partially separable functions
% -----------------------------------------------------------------------------
% diffDemo        Illustrates numerical differentiation
% -----------------------------------------------------------------------------
% democon         Demonstration file for constrained problems
% -----------------------------------------------------------------------------
% demouc          Demonstration file for unconstrained problems
%
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
% Test examples demonstration the use of the Optimization tbx 2.x interfaces
% -----------------------------------------------------------------------------
% testbintprog    Run test of the OPTIM TB 3.x bintprog equivalent
% -----------------------------------------------------------------------------
% testfgoalattain Run test of the OPTIM TB 3.x fgoalattain equivalent
% -----------------------------------------------------------------------------
% testfmincon     Run test of the OPTIM TB 2.x fmincon equivalent
% fmincon_f       Function value routine used by testfmincon
% fmincon_c       Constraint routine used by testfmincon
% -----------------------------------------------------------------------------
% testfminunc     Run test of the OPTIM TB 2.x fminunc equivalent
% fminunc_fg      Function and gradient routine used by testfminunc
% fminunc_fg2     Function and gradient routine used by testfminunc
% -----------------------------------------------------------------------------
% testlinprog     Run test of the OPTIM TB 2.x linprog equivalent
% -----------------------------------------------------------------------------
% testlsqcurvefit Run test of the OPTIM TB 2.x lsqcurvefit equivalent
% curve_rJ        Computes residual and Jacobian for testlsqcurvefit
% -----------------------------------------------------------------------------
% testlsqlin      Run test of the OPTIM TB 2.x lsqlin equivalent
% -----------------------------------------------------------------------------
% testlsqnonlin   Run test of the OPTIM TB 2.x lsqnonlin equivalent
% lsq_rJ          Computes residual and Jacobian for testlsqnonlin
% -----------------------------------------------------------------------------
% testlsqnonneg   Test Matlab lsqnonneg compared to Tlsqnonneg
% -----------------------------------------------------------------------------
% testquadprog    Run test of the OPTIM TB 2.x quadprog equivalent
% -----------------------------------------------------------------------------
% Trunfleq1       Optimization toolbox using the option HessMult
