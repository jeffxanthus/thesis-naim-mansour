% SOL solvers - Definition of default values
%
% -----------------------------------------------------------------------
% function [optPar] = SOLGet(Solver, optType, nObj, nJac, m)
% ------------------------------------------------------------------------
% INPUT:  
% Solver
% optType   String with type
% nObj      Number of nonlinear objective variables
% nJac      Number of nonlinear constraint variables
% m         Number of rows in the constraint matrix
%
% OUTPUT: 
% optPar

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., $Release: 4.9.0$
% Written Sep 3, 2000.    Last modified Aug 1, 2005.

function optPar = SOLGet(Solver, optType, nObj, nJac, m)

if nargin < 5
   m=[];
   if nargin < 4
      nJac=[];
      if nargin < 3
         nObj=[];
         if nargin < 2
            optType=[];
            if nargin < 1
              Solver=[];
end,end,end,end,end

if isempty(m),       m=70; end
if isempty(nJac),    nJac=70; end
if isempty(nObj),    nObj=70; end
if isempty(optType) 
   if strcmpi(Solver,'lssol')
      optType=checkType('lls');
   else
      optType=checkType('con');
   end
end
if isempty(Solver),  Solver='minos'; end

Solver=lower(Solver);

if ischar(optType)
   optType = checkType(optType);
elseif isempty(optType)
   optType = 3;
end

nL=max(nObj,nJac);

optPar=[];

switch lower(Solver)
  case 'snopt'
      e67 = eps^(0.67);     % e67 = 3.25E-11, eps^(2/3) approx
      e08 = 3.0E-13;        % eps^(0.8)
      m42 = eps^0.4;        %  =  m42 = e08^0.5; = 5.48E-7
      m43 = e08^0.333333;   % 6.69E-5
      % Now the stupid mistake is corrected, defalt SNOPT opt tol is 2E-6 now
      % m10 = 1.73E-6;      % max(1E-6,(10*e08)^0.5)
      m10 = 2E-6;           % max(2E-6,(10*e08)^0.5)
      ep10 = 10*eps;
      m25 = eps^0.25;       %  = 1.22E-4
      m71 = min(2000,1+nL);
      % SNOPT

      if nL <= 75
         m47=1;
      else
         m47=0;
      end
      m20=max(10000,20*m); 
      m48 = min(500,1+nL);
      m35=max(1000,3*max(m,nObj)); 
      big=99999999;
      %      [   1     2    3    4    5    6    7    8     9    10
      optPar=[   0     1   -1   -1    1    1    1    0  1E-6   m10 ...
              1E-6  1E-6    0    1 nObj    1 nJac    1   0.9     0 ...
               0.1   0.9    5    5  m25  e67 ep10    3   100 10000 ...
                 1     0    0  1E6  m35  500    2  big     3     1 ...
               e08   m42  m43    0 1E20 1E15  m47  m48    20   big ...
                60 10000   50  100    0    0    0    0     0     0 ...
                 0     0    0    0   99    0 1E-2  100     0   0.1 ...
               m71 ...
             ];
      % SNOPT for QP ? LP ? LLS ?
      if checkType('lp',optType)
         optPar(18)=2;
         optPar(23)=100;
         optPar(24)=10;
         optPar(31)=10;
      %elseif checkType('qp',optType)
         %optPar(35)=-900;
      end

  case 'minos'
      e67 = eps^(0.67);     % e67 = 3.25E-11, eps^(2/3) approx
      m35 = max(m,1000); 
      e08 = 3.0E-13;        % eps^(0.8)
      m42 = e08^0.5;        % 5.48E-8
      m43 = e08^0.333333;   % 6.69E-5
      m10 = 1.73E-6;        % max(1E-6,(10*e08)^0.5)
      % MINOS
      m30=3*max(1,m + any(optType==[2 7 8]))+10*nL;
      m25 = eps^0.25;       %  = 1.22E-4
      %      [-999  -999 -999 -999 -999 -999 -999 -999  -999  -999
      %      [   1     2    3    4    5    6    7    8     9    10
      optPar=[   0  -900   -1   -1  100  100    1    0  1E-6   m10 ...
              1E-6  -900    0    1 nObj    1 nJac    1   0.9     0 ...
               0.1   0.1    5    5  m25  e67  e67    3   0.0   m30 ...
                 1     0    1  1.0   50   40  2.0  2.0     3  0.01 ...
               e08   m42  m43    0 1E20 1E10   50   50  -900  -900 ...
                60 10000   50  100    0    0    0    0     0     0 ...
                 0     0    0    0    0    0 1E-10   0     0   0.5 ...
                 0 ...
             ];
      if checkType('lp',optType)
         optPar(6)=100;
         %optPar(9)=-900;
         optPar(10)=1E-6;
         %optPar(13:17)=-900;
         optPar(18)=2;

         %optPar(22)=-900; % Line search tol
         optPar(23)=100;
         optPar(24)=10;
         optPar(31)=10;
         %optPar(32:50)=-900;
         optPar(44)=1;
         optPar(53)=100;
      elseif checkType('qp',optType)
         %optPar([9,16:17,33:38,40])=-900;
      end

  case 'lp-minos'
      % MINOS
      e67 = eps^(0.67);     % e67 = 3.25E-11, eps^(2/3) approx
      m30 = 3*m;
      m25 = eps^0.25;       %  = 1.22E-4
      % must set element 33 for GUI
      %      [   1     2    3    4    5    6    7    8     9   10
      optPar=[   0  -900   -1   -1  100  100    1    0  -900 1E-6 ...
              1E-6  -900 -900 -900 -900 -900 -900    2   0.9    0 ...
               0.1  -900  100   10  m25  e67  e67    3   0.0  m30 ...
                10  -900    1 -900 -900 -900 -900 -900  -900 -900 ...
              -900  -900 -900    1 -900 -900 -900 -900  -900 -900 ...
                60 10000  100  100    0    0    0    0     0    0 ...
                 0     0    0    0    0    0 1E-10   0     0   0.5 ...
                 0 ...
             ];
  case 'qp-minos'
      e67 = eps^(0.67);     % e67 = 3.25E-11, eps^(2/3) approx
      % MINOS-QP
      m10 = 1.73E-6;        % max(1E-6,(10*e08)^0.5)
      m30 = 3*m+10*nL;
      e08 = 3.0E-13;        % eps^(0.8)
      m42 = e08^0.5;        % 5.48E-8
      m43 = e08^0.333333;   % 6.69E-5
      m25 = eps^0.25;       %  = 1.22E-4
      % must set element 33 for GUI
      %      [   1     2    3    4    5    6    7    8     9    10
      optPar=[   0  -900   -1   -1  100  100    1    0  -900   m10 ...
              1E-6  -900    0    1 nObj -900 -900    1   0.9     0 ...
               0.1   0.9    5    5  m25  e67  e67    3   0.0   m30 ...
                 1     0    1 -900 -900 -900 -900 -900     3  -900 ...
               e08   m42  m43    0 1E20 1E10   50   50  -900  -900 ...
                60 10000   50  100    0    0    0    0     0     0 ...
                 0     0    0    0    0    0 1E-10   0     0   0.5 ...
                 0 ...
             ];
  case {'qpopt','lpopt'}
      % QPOPT
      optPar(1:65)=-900;
      mxIter=max(50,5*(nObj+m));
      % macheps in QPOPT is different from SNOPT, MINOS and Matlab
      macheps = 2^(-53);  % eps in Matlab is = 2^(-52); 
      % macheps = 0.5*eps;  % eps in Matlab is = 2^(-52); 
      tolfea = sqrt(macheps);  % 1.05E-8, epspt5
      tolOpt = macheps^(-0.8); % 1.72E-13, epspt8
      tolrnk = 100*macheps;    % 1.11E-14
      optPar([1 10 11 21 27 30 33 36 45 48 51 52])= ...
         [ 0 tolfea tolOpt 0.01 tolrnk mxIter 0 mxIter 1E20 nObj 50 5];

      % Currently no difference between LP and QP problems

      % 1.  PRINT LEVEL  optPar(1)  = 10;
      % 10. OPTIMALITY TOLERANCE   optPar(10) = 1.72E-13; 
      % 11. FEASIBILITY TOLERANCE  optPar(11) = 1.05E-8;
      % 21. CRASH TOLERANCE        optPar(21) = 0.01; 
      % 27. RANK TOLERANCE         optPar(27) = 1.11E-14;
      % 30. ITERATION LIMIT        optPar(30) = max(50,5*(n+m));
      % 33. MIN SUM. IF SET TO 1, minimize infeasibilities before return (def 0)
      % 36. FEAS PHASE ITERATIONS LIMIT  optPar(36) = max(50,5*(n+m));
      % 45. INFINITE STEP SIZE optPar(45) = 1E20;
      % 48. MAX DEGREES OF FREEDOM, ONLY USED IF HESSIAN ROWS == N (Def = n)
      % 51. CHECK FREQUENCY  optPar(51) = 50;
      % 52. EXPAND FREQUENCY optPar(52) = 5;
      % 3,4,47 set by Matlab interface
      % 3.  PRINT FILE           0        0                 Fortran Unit #
      %           SET BY INTERFACE IF PrintFile is given
      % 4.  SUMMARY FILE         0        0                 Fortran Unit #
      %           SET BY INTERFACE IF SummFile  is given
      % 47. HESSIAN ROWS         0        n         n       0 if FP or LP
      %     IMPLICITLY GIVEN BY THE DIMENSIONS OF H IN THE CALL FROM MATLAB
  case 'npsol'
      % NPSOL
      % nnObj == n, nnJac == nonlin cons (ncnln), m = linear cons (nclin)

      seps=sqrt(eps);                    % 1.1E-8
      e08 = 3.0E-13;                     % eps^(0.8)
      m30=max(50,3*(nObj+m) + 10*nJac);  % 3*(n+nclin) + 10*ncnln 
      m36=max(50,3*(nObj + m + nJac));   % 3*(n+nclin+ncnln) 
      e09 = eps^(0.9);
      %      [   1     2    3    4    5    6    7    8     9    10
      optPar=[   0     0   -1   -1 -999 -999 -999 -999  seps   e08 ...
              seps  -999    0    1 nObj    1 nObj -999  -999  -999 ...
              0.01   0.9 -999 -999 -999 -999 -999 -999  -999   m30 ...
              -999  -999 -999 -999 -999  m36 -999 -999     3  -999 ...
               e08  -900 -900 -999 1E10 1E10 -999 -999  -999     0 ...
              -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
              -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
              -999  ...
             ];
      %optPar(51:65)=-999;
  case 'nlssol'
      % NLSSOL. Difference NPSOL is #47 and #48
      % nnObj == n, nnJac == nonlin cons (ncnln), m = linear cons (nclin)

      % Currently no difference between LLS, LP and QP problems

      seps=sqrt(eps);                    % 1.1E-8
      e08 = 3.0E-13;                     % eps^(0.8)
      m30=max(50,3*(nObj+m) + 10*nJac);  % 3*(n+nclin) + 10*ncnln 
      m36=max(50,3*(nObj + m + nJac));   % 3*(n+nclin+ncnln) 
      e09 = eps^(0.9);
      %      [   1     2    3    4    5    6    7    8     9    10
      optPar=[   0     0   -1   -1 -999 -999 -999 -999  seps   e08 ...
              seps  -999    0    1 nObj    1 nObj -999  -999  -999 ...
              0.01   0.9 -999 -999 -999 -999 -999 -999  -999   m30 ...
              -999  -999 -999 -999 -999  m36 -999 -999     3  -999 ...
               e08  -900 -900 -999 1E10 1E10    1    2  -999  -999 ...
              -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
              -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
              -999  ...
             ];
      %optPar(51:65)=-999;
  case 'lssol'
      % LSSOL  -   Linear least squares, LP and QP
      optPar(1:65)=-900;
      mxIter=max(50,5*(nObj+m));
      optPar([1 10 11 21 27 30 36 45 46])= ...
         [  0 1.72E-13 1.05E-8 0.01 1.11E-14 mxIter mxIter 1E20 1E20];


      if ~checkType('lls',optType)
         % QP or LP, lower rank tolerance
         optPar(27) = 1.05E-7;
      end
     

      % 1.  PRINT LEVEL  optPar(1)  = 10;
      % 10. OPTIMALITY TOLERANCE   optPar(10) = 1.72E-13;
      % 11. FEASIBILITY TOLERANCE  optPar(11)  = 1.05E-8;
      % 21. CRASH TOLERANCE optPar(21) = 0.01; 
      % 27. RANK TOLERANCE  optPar(27) = 1.11E-14;
      % 30. ITERATION LIMIT optPar(30) = max(50,5*(n+m));
      % 36. FEAS PHASE ITERATIONS LIMIT  optPar(36) = max(50,5*(n+m));
      % 45. INFINITE STEP SIZE optPar(45)  = 1E20;
      % 46. INFINITE BOUND SIZE optPar(46) = 1E20;

      % LP/QP 
      % 27. RANK TOLERANCE         optPar(27)  = 1.05E-7;

  case 'sqopt'
      % SQOPT

      e67 = eps^(0.67);     % e67 = 3.25E-11, eps^(2/3) approx
      m48 = min(500,1+nObj);
      m30 = 3*m;
      m25 = eps^0.25;       %  = 1.22E-4
      ep10 = 10*eps;
      %      [   1     2    3    4    5    6    7    8     9    10
      optPar=[   0     0   -1   -1  100  100    1    0  -900  -900 ...
              1E-6  1E-6 -900 -900 -900 -900 -900    1   0.9     0 ...
               0.1  -900    5    5  m25  e67 ep10    3     1   m30 ...
                 1     0 -900 -900 -900 -900 -900 -900     3  -999 ...
              -900  -900 -900 -900 1E20 -900 -900  m48     1     2 ...
                60 10000   50  100    0    0    0    0     0     0 ...
                 0     0    0    0    0    0    0    0     0     0 ...
                 0 ...
             ];
      if checkType('lp',optType)
         optPar(18)=2;
         optPar(23)=100;
         optPar(24)=10;
         optPar(31)=10;
      end
  otherwise
      %disp('ERROR in SOLGet. Wrong solver name')
      optPar=[];
end

% MODIFICATION LOG:
%
% 000902 hkh  Written
% 000915 hkh  Modified for npsol
% 001102 hkh  Added lssol
% 001110 hkh  Interchange #34 and #44 for snopt, gives 0-1 for #44, as minos
% 010403 hkh  Superbasics from #38 to #48 for snopt, sqopt
% 010903 hkh  Change default print level for snopt to 0, add optPar(63).
% 040929 hkh  Default line search tolerance in MINOS should be 0.1
% 041117 frhe Change default print level for all solvers to 0.
% 041221 hkh  Correct qpopt defaults for opt, feas and rank tol
% 050605 hkh  Correct snopt defaults
% 050606 hkh  Expand optPar to handle SNOPT 7
% 050606 hkh  Revision for update of MINOS to 71 parameters
% 050801 med  isstr renamed to ischar