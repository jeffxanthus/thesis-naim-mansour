% nsdpAssign provides a direct way of setting up linear semidefinite
% programming (SDP) problems with linear or bilinear matrix objectives
% in the TOMLAB (TQ) format.
%
% The information is put into the Tomlab input problem structure Prob.
%
% Prob = nsdpAssign( .... )
%
% It is then possible to solve the SDP problem using PENNON
%
% Result = tomRun('pennon',Prob);
%
% -----------------------------------------------------
%
% The NLP-SDP minimization problem is defined as
%
%        min   f(x,Y)
%        x,Y
%
%        where x is a vector variable in R^n
%              Y are k matrix variables in S^pi, i = 1,...,k
%
%   subject to
%
%          x_L <=   x     <= x_U         x_L,x_U in R^n
%          Y_L <=   Y     <= Y_U         Y_L, Y_U in R^[p_i x p_i] i=1..k
%   lambda_L*I <=   Y     <= lambda_U*I  lambda_L_i and lambda_U_i are k
%                                        bounds on the eigenvalues of Y_i
%          b_L <= A*(x,Y) <= b_U         b_L,b_U in R^ml, A in R^[ml x n]
%          c_L <= c(x,Y)  <= c_U
%
% Equality equations: Set b_L==b_U
%                         c_L==c_U
% Fixed variables:    Set x_L==x_U
%                     Set Y_L==Y_U
%
% -----------------------------------------------------
%
% Syntax of nsdpAssign:
%
% function Prob = nsdpAssign(f, g, H, HessPattern, x_L, x_U, x_0, NSDP, ...
%        A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, Name, fLowBnd)
%
% INPUT
%
% O B J E C T I V E   F U N C T I O N
%
% f           Name of the function that computes the function value f(x,Y)
% g           Name of the function that computes the n x 1 gradient vector
% H           Name of the function that computes the n x n Hessian matrix
% HessPattern n x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the Hessian and ones indicate values that might
%             be non-zero. If empty indicates estimation of all elements
%             HessPattern is used when estimating the Hessian numerically.
%             Estimated before solve, if Prob.LargeScale==1, HessPattern==[]
%
% S T A N D A R D   V A R I A B L E S
%
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as a nx1  Inf vector.
% x_0         Starting values, default nx1 zero vector
%
% M A T R I X   V A R I A B L E S
%
% NSDP          Structure array describing the SDP matrix variables.
%               i = 1, ..., k matrix variables
%
% NSDP(i).Y     p_i x p_i matrix variables (or empty)
%               Variable elements are indicated as NaN.
%               Constants with their respective values.
%
% NSDP(i).Y_L   p_i x p_i matrices (or empty)
%               representing the lower bounds on the elements in Y
%               If [] set as a p_i x p_i -Inf matrix.
% NSDP(i).Y_U   p_i x p_i matrices (or empty)
%               representing the upper bounds on the elements in Y
%               If [] set as a p_i x p_i Inf matrix.
%
% NSDP(i).lambda_L  The lower bound on the eigenvalue of NSDP(i).Y
%                   If [] set as a p_i x 1 zero vector.
% NSDP(i).lambda_U  The upper bound on the eigenvalue of NSDP(i).Y
%                   If [] set as a p_i x 1 Inf vector.
%
% NSDP(i).Y_0 Starting values, default p x p zero Matrix
%
% NSDP(i).Type  Type of matrix constraint
%                0  Standard (Constraints may be infeasible)
%                1  Lower strict (The lower bounds will always be feasible)
%                2  Upper strict (The upper bounds will always be feasible)
%                3  Lower/upper strict (both bounds will always be
%                feasible)
%                4  Slack (The block is a slack variable).
%
% L I N E A R   C O N S T R A I N T S
% A           mA x n matrix A, linear constraints b_L <= A*x <= b_U. Dense or sparse
% b_L         Lower bound vector in linear constraints b_L <= A*x <= b_U.
% b_U         Upper bound vector in linear constraints b_L <= A*x <= b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints
% dc          Name of function that computes the constraint Jacobian mN x n
% d2c         Name of function that computes the second part of the
%             Lagrangian function (only needed for some solvers)
%             See the help gateway routine nlp_d2c for an explanation of d2c
%
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
%
% c_L         Lower bound vector in nonlinear constraints c_L <= c(x) <= c_U.
% c_U         Upper bound vector in nonlinear constraints c_L <= c(x) <= c_U.
%
%               b_L, b_U, c_L, c_U, x_L, x_U must either be
%               empty or of full length
%
% Name          The name of the problem (string)
%
% fLowBnd       A lower bound on the function value at optimum. Default -1E300
%               Use [] if not known at all. SDP solvers do not use fLowBnd

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written October 2, 2008. Last modified Oct 23, 2009.

function Prob = nsdpAssign(f, g, H, HessPattern, x_L, x_U, x_0, NSDP, ...
        A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
			  Name, fLowBnd)


if nargin < 20
  fLowBnd=[];
  if nargin < 19
    Name=[];
    if nargin < 18
      c_U=[];
      if nargin < 17
        c_L=[];
        if nargin < 16
          ConsPattern=[];
          if nargin < 15
            d2c=[];
            if nargin < 14
              dc=[];
              if nargin < 13
                c=[];
                if nargin < 12
                  b_U=[];
                  if nargin < 11
                    b_L=[];
                    if nargin < 10
                      A=[];
                      if nargin < 9
                        SDP = [];
                        if nargin < 8
                          x_0=[];
end, end, end, end, end, end, end, end, end, end, end, end, end

% Do we need a new probType for NSDP?
probType = checkType('sdp');

Prob          = ProbDef(1);
Prob.P        = 1;
Prob.probFile = 0;
Prob.probType = probType;
Prob.Name     = deblank(Name);

% Determine the number of standard (real) variables
n = max([length(x_L),length(x_U),length(x_0)]);

% Step through NSDP and extend x_0, x_L and x_U
Prob.PENOPT.NSDP = NSDP;
if isstruct(NSDP)
   nmv = length(NSDP);
   for i = 1:nmv
      
      % i as string for errors and warnings
      i_str = num2str(i);

      % Check for required field Y
      % Y defines the structure of matrix variable i
      if ~isfield(NSDP(i),'Y')
         error(strcat('Field Y missing in struct Prob.PENOPT.NSDP(',i_str,')'));
      else
         Y = NSDP(i).Y;
      end

      if isempty(Y)
         error(strcat('Field Y missing in struct Prob.PENOPT.NSDP(',i_str,')'));
      end
      % Y must be quadratic
      p = size(Y,1);
      if p ~= size(Y,2)
         error(['Field Y in struct Prob.PENOPT.NSDP(',num2str(i),') must be quadratic']);
      end

      % Y should be symmetric (but only upper triangle is extracted)
      Yt = Y';
      Y_num_ix = find(~isnan(Y));
      Y_NaN_ix = find(isnan(Y));  % Since NaN ~= NaN we need to compare indices
      Yt_NaN_ix = find(isnan(Yt));
      YT_NaN_ix = find(isnan(Y));
      if any(find(Y(Y_num_ix) ~= Yt(Y_num_ix)) | any(Y_NaN_ix ~= Yt_NaN_ix))
         warning(strcat('Field Y is not symmetric in struct Prob.PENOPT.NSDP(',i_str,')'));
      end

      % Find nonzero elements of matrix variables (both variables and nonzero constants)
      triu_Y = triu(Y);
      index_NaN = find(isnan(triu_Y));
      index_Nz = find(triu_Y ~= 0 & ~isnan(triu_Y));
      Y_var_index = sort([index_NaN; index_Nz]);
      
      n_Yi = length(Y_var_index);

      n = n + n_Yi;

      if n_Yi < 1
         error(strcat('Field Y in struct Prob.PENOPT.NSDP(',i_str,') is a zero matrix'));
      end    

      % Check Y_L and set -Inf if empty
      Y_L = DefPar(NSDP(i), 'Y_L');
      if isempty(Y_L)
         Y_L = Y;
         Y_L(index_NaN) = -Inf;
         Y_L(index_Nz) = Y(index_Nz);
      else
         if size(Y_L,1) ~= size(Y_L,2) | size(Y_L,1) ~= p
            error(['Field Y_L in struct Prob.PENOPT.NSDP(',i_str,') must be of the same size as Y']);
         end
      end
      
      % Check Y_U and set Inf if empty
      Y_U = DefPar(NSDP(i),'Y_U');
      if isempty(Y_U)
         Y_U = Y;
         Y_U(index_NaN) = Inf;
         Y_U(index_Nz) = Y(index_NaN);
      else
         if size(Y_U,1) ~= size(Y_U,2) | size(Y_U,1) ~= p
            error(['Field Y_U in struct Prob.PENOPT.NSDP(',i_str,') must be of the same size as Y']);
         end
      end
      
      % Check Y_0 and set 0 if empty
      Y_0 = DefPar(NSDP(i),'Y_0');
      if isempty(Y_0)
         Y_0 = Y;
         Y_0(Y_var_index) = min(max(0,Y_L(Y_var_index)),Y_U(Y_var_index));
      else
         if size(Y_0,1) ~= size(Y_0,2) | size(Y_0,1) ~= p
            error(['Field Y_0 in struct Prob.PENOPT.NSDP(',i_str,') must be of the same size as Y']);
         end
      end
      
      % Store variable bounds
      x_L = [x_L; Y_L(Y_var_index)];
      x_U = [x_U; Y_U(Y_var_index)];

      % Store initial variable values
      x_0 = [x_0; Y_0(Y_var_index)];
  end
end

% Redetermine the number of variables
n = max([length(x_L),length(x_U),length(x_0)]);

% Determine the number of linear constraints
mL = max(length(b_L),length(b_U));

% Determine the number of nonlinear constraints
mN = max(length(c_L),length(c_U));

% Store in Prob
Prob.P    = 1;
Prob.N    = n;

% Check size of A
if ~isempty(A)
  Am = size(A,1);
  An = size(A,2);
else
  Am = 0;
  An = 0;
end

Prob.HessPattern = HessPattern;
Prob.ConsPattern = ConsPattern;

if ~isempty(fLowBnd)
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
   
if mN > 0
   if isempty(c_L) 
      Prob.c_L=-Inf*ones(mN,1);
   else
      Prob.c_L=full(double(c_L(:)));
   end
   if isempty(c_U) 
      Prob.c_U=Inf*ones(mN,1);
   else
      Prob.c_U=full(double(c_U(:)));
   end
   if any(Prob.c_L>Prob.c_U)
      error('c_L and c_U have crossover values');
   end
end

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

Prob.USER.f   = f;
Prob.USER.g   = g;
Prob.USER.H   = H;
Prob.USER.c   = c;
Prob.USER.dc  = dc;
Prob.USER.d2c = d2c;

if ~isempty(NSDP)
  Prob = tomFiles(Prob, 'nsdp_f', 'nsdp_g', 'nsdp_H', 'nsdp_c', 'nsdp_dc', 'nsdp_d2c');
end

% MODIFICATION LOG
%
% 081002 bjo Written (based on conAssign bmiAssign and sdpAssign)
% 091023 med Corrected errors and removed code in checkAssign