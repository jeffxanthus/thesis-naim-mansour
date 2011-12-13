function OPTIONS=goptions(inopt)
% GOPTIONS Default parameters used by the optimization routines.
% Expanded version of the MATLAB parameter vector defined by foptions.
%	In TOMLAB: Used for some old parts
%	In MATLAB itself: FMIN and FMINS.
%	In the Optimization Toolbox:
%	    FMINU, CONSTR, ATTGOAL, MINIMAX, LEASTSQ, FSOLVE.
% The parameters are:
%     OPTIONS(1)  - Display parameter (Default 0). 1 displays some results
%                   1 = Result of optimization.  2 = Every iteration
%                   3 = Extra info, residual...  4 = Maximal info, Jacobians...
%     OPTIONS(2)  - Termination tolerance for X.(Default=1E-8).
%     OPTIONS(3)  - Termination tolerance on F.(Default=1E-8). Dir.derivative
%     OPTIONS(4)  - Termination criterion on constraint violation(Default=1E-6)
%     OPTIONS(5)  - Algorithm: Strategy:  Not always used.
%     OPTIONS(6)  - Algorithm: Optimizer: Not always used. 
%     OPTIONS(7)  - Algorithm: Line Search Algorithm. (Default =0 quad. 1=cubic)
%     OPTIONS(8)  - Function value. (Lambda in goal attainment. )
%     OPTIONS(9)  - Set to 1 if you want to check user-supplied gradients
%     OPTIONS(10) - Number of Function and Constraint Evaluations.
%     OPTIONS(11) - Number of Function Gradient Evaluations.
%     OPTIONS(12) - Number of Constraint Evaluations
%     OPTIONS(13) - Number of equality constraints. 
%     OPTIONS(14) - Maximum number of iterations. (Default 100*no. of variables)
%     OPTIONS(15) - Used in goal attainment for special objectives. 
%     OPTIONS(16) - Minimum change in variables for finite difference gradients.
%     OPTIONS(17) - Maximum change in variables for finite difference gradients.
%     OPTIONS(18) - Step length. (Default 1 or less). 
% Expanded parameters
%     OPTIONS(19) - Penalty parameter. Default 100.
%     OPTIONS(20) - eps_g. Termination tolerance on G.(Default=1E-6).
%     OPTIONS(21) - Line search accuracy sigma.  (Default=0.9)
%                   sigma = 0.9 inexact search. sigma = 0.1  exact search
%     OPTIONS(22) - f_Low. Lower bound on function value. Used in line search.
%     OPTIONS(23) - eps_rank.Rank test tolerance. Used in subspace minimization.
%     OPTIONS(24) - wait. Flag. If set use pause after iteration printout.
%     OPTIONS(25) - eps_fabs. Absolute convergence tolerance in function f.
%     OPTIONS(26) - iter k. Number of main (major) iterations.
%     OPTIONS(27) - Number of minor iterations.
%     OPTIONS(28) - EXIT flag, convergence to local min = 0. Otherwise errors.
%     OPTIONS(29) - INFORM, information parameter, type of convergence.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written May 17, 1995.  Last modified June 23, 1999.
%

OPTIONS=zeros(1,29);

if nargin < 1;
   inopt = []; 
else
   OPTIONS(1:length(inopt))=inopt;
end
def_opt=  [0,1E-8,1E-8,1E-6,0,0,0,0,0,0,0,0,0,0,0,1E-8,0.1,0,...
          100,1E-6,0.9,-1E-100,1E-10,0,-1E-300,0,0,0,0];
OPTIONS=OPTIONS+(OPTIONS==0).*def_opt;
