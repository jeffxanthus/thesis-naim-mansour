% glbSolve implements a modified version of the algorithm DIRECT. Reference:
%
% D. R. Jones, C. D. Perttunen and B. E. Stuckman,
% Lipschitzian Optimization Without the Lipschitz Constant,
% JOTA  Vol. 79, No. 1, October 1993.
%
% glbSolve solves box-bounced global optimization problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U   , where x, x_L, x_U is in R^n
%
% Calling syntax:
%
% function Result = glbSolve(Prob,varargin)
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string.
%             A call to tomFiles.m or glcAssign.m is setting this field. 
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   PriLevOpt Print Level 
%             0 = silent. 1 = some printing. 2 = print each iteration
%   WarmStart If true, >0, glbSolve reads the output from the last run
%             from the mat-file glbSave.mat, and continues from the last run.
%   MaxCPU    Maximal CPU Time (in seconds) to be used
% optParam    Structure in Prob, Prob.optParam 
%             Defines optimization parameters. Fields used: 
%  IterPrint  Print one line each iteration
%  MaxIter    Maximal number of iterations, default max(5000,n*1000);
%  MaxFunc    Maximal number of function evaluations, default max(10000,n*2000) 
%  EpsGlob    Global/local weight parameter, default 1E-4.
%  fGoal      Goal for function value, if empty not used
%  eps_f      Relative accuracy fTol for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol, if fGoal~=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal==0
%
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with optimal points as columns.
%  f_k      The best function value found so far
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  maxTri   Maximum size of any triangle
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag 0 = Normal termination, max number of iterations /func.evals reached
%           1 = Some bound, lower or upper is missing
%           2 = Some bound is inf, must be finite
%           4 = Numerical trouble determining optimal rectangle, empty set
%               and cannot continue
%  Inform   0 = Normal Exit
%           1 = Function value f is less than fGoal
%           2 = Absolute function value f is less than fTol, only if fGoal = 0
%            or Relative error in function value f is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           9 = Max CPU Time reached
%
% To make a warm start possible, glbSolve saves the following information in
% the file glbSave.mat (for internal solver use only):
%   C       Matrix with all rectangle centerpoints, in [0,1]-space.
%   D       Vector with distances from centerpoint to the vertices.
%   DMin    Row vector of minimum function value for each distance 
%   DSort   Row vector of all different distances, sorted.
%   E       Computed tolerance in rectangle selection
%   F       Vector with function values.
%   L       Matrix with all rectangle side lengths in each dimension.
%   Name    Name of the problem. Used for security if doing warm start
%   glbfMin   Best function value found at a feasible point.
%   iMin    The index in D which has lowest function value, i.e. the
%           The rectangle which minimizes (F - glbfMin + E)./D where
%           E = max(EpsGlob*abs(glbfMin),1E-8)
%
% USAGE:
%
% Let the name of the problem be "GLBF Test"
% The function GLBF is best written as
%     function f = GLBF(x, Prob) 
% Then any information, say b and C is easily sent to GLBF in the
% Prob structure by the call (assume bounds x_L and x_U are initialized) 
%
%      Prob   = glcAssign('GLBF',x_L,x_U,'GLBF Test');
%      Prob.user.b = b; Prob.user.C=C;        % example of extra user data
%
%      % Default values are now set for PriLevOpt, and structure optParam
%      % To change a value, examples are shown on the two next lines
%      Prob.optParam.MaxFunc = 500; % Change max number of function evaluations 
%      Prob.optParam.MaxIter = 100; % Change the number of iterations to 100
%
% Driver call, including printing with level 2:
%      Result = tomRun('glbSolve',Prob,2);
%
% Direct solver call:
%      Result = glbSolve(Prob);
%      PrintResult(Result,2);
%           
% The user function GLBF is written
%
%      function f = GLBF(x, Prob) 
%      b = Prob.user.b; C = Prob.user.C;
%      f = "some function of x, b and C"
%
% It is also possible to use the function format
%      function f = GLBF(x) 
% or sending the extra parameters as additional input 
%
%      Result = glbSolve(Prob,b,C);
%
%      function f = GLBF(x,Prob,b,C) 
%
% NOTE! If additional parameters are sent, Prob must be the second input 
% parameter to GLBF
%
% To make a restart, just set the restart flag, and call glbSolve once again:
%
%      Prob.WarmStart = 1;
%      Result = glbSolve(Prob);  % Assuming no extra parameters in the call
%      PrintResult(Result);
%
% IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM:
%
% The DIRECT algorithm only reaches the variable bounds in the limit.
% Therefore convergence for global optimum where components are on the bounds
% is slow. 
% One remedy is to reduce lower bounds with a tolerance, say 1E-4, and add
% a similar tolerance 1E-4 to the upper bounds that might be reached.
% Another possibility is to fix a variable on its bound by setting the lower
% and upper bounds equal.
%
% Always try to reduce the dimension as much as possible when using the
% DIRECT algorithm, and try to shrink the box, defined by the lower and
% upper bounds, as much as possible.
%
% DIRECT has good global properties, but is slow in achieving local
% convergence with high accuracy. It is advised to use the best point found
% by DIRECT as initial value for a local search to increase accuracy.
% e.g.:
%      Result   = tomRun('glbSolve',Prob,2);
%      Prob.x_0 = Result.x_k(:,1);          
%      Result2  = tomRun('npsol',Prob,2);

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.6.0$
% Written Sep 1, 1998.    Last modified Oct 3, 2006.

function Result = glbSolve(Prob, varargin)

if nargin < 1
   error('glbSolve needs input structure Prob');
end

solvType=checkType('glb');

Prob=ProbCheck(Prob,'glbSolve',solvType);

Prob = iniSolve(Prob,solvType,0,0);

MaxCPU = Prob.MaxCPU;

if isempty(Prob.x_L) | isempty(Prob.x_U)
   disp('glbSolve requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.Inform   = 0;
   Result.ExitText = 'glbSolve requires both lower and upper variable bounds';
   Result=endSolve(Prob,Result);
   return;
end

PriLev    = Prob.PriLevOpt;          % Print level
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds
MaxIter   = Prob.optParam.MaxIter;   % Number of iterations
MaxFunc   = Prob.optParam.MaxFunc;   % Number of function evaluations
EpsGlob   = Prob.optParam.EpsGlob;   % global/local weight parameter. 
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration
fGoal     = Prob.optParam.fGoal;     % Goal for f(x).
fTol      = Prob.optParam.eps_f;     % Relative tolerance for fGoal
n         = Prob.N;

% Safeguard
if isempty(MaxIter) | MaxIter < 0
   MaxIter = max(5000,n*1000);
end
if isempty(MaxFunc) | MaxFunc < 0
   MaxFunc = max(10000,n*2000);
end

nFunc     = 0;
convFlag  = 0;

Result                 = ResultDef(Prob);
Result.Solver          = 'glbSolve';
Result.SolverAlgorithm = 'DIRECT - Lipschitzian Optimization';

if any(isinf(x_L)) | any(isinf(x_U))
   disp('glbSolve solves box-bounded problems.');
   disp('Found some bound to be Inf');
   Result.ExitFlag = 2;
   Result.Inform   = 0;
   Result.ExitText =str2mat('glbSolve solves box-bounded problems' ...
                           ,'Found some bound to be Inf');
   Result=endSolve(Prob,Result);
   return
end

x_D = x_U - x_L;

n   = length(x_L);  % Problem dimension

%tol1  = 1E-16;
tol2 = 1E-8;

%
%  STEP 1, Initialization
%

if Prob.WarmStart
   % Restart with values from previous run.

   load('glbSave.mat','Name','C','F','D','L','DSort','DMin','glbfMin','E','iMin')
   Name1 = Prob.Name;               % Name for the problem
   if strcmp(Name1,Name)

      m = length(F);
   
      if PriLev > 0
         fprintf('\n Restarting with %d sampled points from previous run\n',m);
      end
   else
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files glbSave.mat?\n');
      end
   end
end

if ~Prob.WarmStart
   Name = deblank(Prob.Name);  % Problem name
   % No restart, set first point to center of the unit hypercube.
   m = 1;             % Current number of rectangles
   %C = ones(n,1)./2;  % Matrix with all rectangle centerpoints
   % All C_coordinates refers to the n-dimensional hypercube. 
   mxDim = min(MaxFunc+n*20,9999);
   C = [0.5*ones(n,1), zeros(n,mxDim)]; 
   
   %x_m = x_L + C.*x_D;   % Transform C to original search space
   x_m = x_L + 0.5*x_D;   % Start with mid point
   glbfMin = nlp_f( x_m, Prob, varargin{:});  % Function value at x_m
   if isnan(glbfMin), glbfMin = realmax; end
   f_0 = glbfMin;
   nFunc=nFunc+1;
   iMin = 1; % The rectangle which minimizes (F - glbfMin + E)./D where
             % E = max(EpsGlob*abs(glbfMin),1E-8)

% Matrix with all rectangle side lengths in each dimension
   L = [0.5*ones(n,1), zeros(n,mxDim)];  
   % Vector with distances from centerpoint to the vertices
   %D = sqrt(sum(L.^2));  
   D = [sqrt(n/4),zeros(1,mxDim)];  
   % Vector with function values
   F = [glbfMin,zeros(1,mxDim)];            

   DSort = D(1);     % Row vector of all different distances, sorted
   DMin  = glbfMin(1);  % Row vector of minimum function value for each distance   
end

c23 = 2/3;

% ITERATION LOOP
Iter   = 1;           % Iter is the iteration counter
cpumax = 0;
TIME0  = Prob.TIME0;

while Iter <= MaxIter & nFunc < MaxFunc & convFlag == 0 

   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end
   
   %
   %  STEP 2  Identify the set S of all potentially optimal rectangles
   %
   S = [];  % Set of all potentially optimal rectangles
   
   idx = find(DSort==D(iMin));
   if isempty(idx)
      if PriLev >= 0
         fprintf('\n WARNING: Numerical trouble when determining S_1\n');
      end
      Result.ExitFlag = 4;
      Result.Inform   = 0;
      Result.ExitText =str2mat('Numerical trouble determining S1' ...
                           ,'Cannot continue');
      Result=endSolve(Prob,Result);
      return;
   end

   S_1 = zeros(1,length(DSort)-idx+1);
   mz = 0;
   for i = idx : length(DSort)
      idx2 = find(  F(1:m)==DMin(i)  &  D(1:m)==DSort(i)  );
      S_1(mz+1:mz+length(idx2)) = idx2;
      mz = mz + length(idx2);
   end
   %S_1 = [];
   %for i = idx : length(DSort)
   %   idx2 = find(  F(1:m)==DMin(i)  &  D(1:m)==DSort(i)  );
   %    S_1 = [S_1 idx2];
   %end
   % S_1 now includes all rectangles i, with D(i) >= D(iMin)
   % and F(i) is the minimum function value for the current distance.
     
   % Pick out all rectangles in S_1 which lies below the line passing through
   % the points: ( D(iMin), F(iMin) ) and the lower rightmost point.
   S_2 = [];
   if length(DSort)-idx > 1
      a1 = D(iMin);
      b1 = F(iMin);
      a2 = DSort(length(DSort));
      b2 = DMin(length(DSort));
      % The line is defined by: y = slope*x + const
      slope = (b2-b1)/(a2-a1);
      const = b1 - slope*a1;
      S_2 = [iMin];
      for i = 1 : length(S_1)
         j = S_1(i);
         if j ~= iMin
            if F(j) <= slope*D(j) + const + tol2 
               S_2 = [S_2 j];
            end   
         end
      end
      % S_2 now contains all points in S_1 which lies on or below the line
        
      % Find the points on the convex hull defined by the points in S_2
      xx = D(S_2);
      yy = F(S_2);
      %h  = conhull(xx,yy); % conhull is an internal subfunction
      h  = tomsol(11,xx,yy); % conhull is an internal subfunction
      S_3 = S_2(h);
   else
      S_3 = S_1;
   end
   S = S_3;
   
   
   %  STEP 3, 5  Select any rectangle j in S
   for jj = 1:length(S)  % For each potentially optimal rectangle     
      j = S(jj);
      
      %
      %  STEP 4  Determine where to sample within rectangle j and how to
      %          divide the rectangle into subrectangles. Update glbfMin
      %          and set m=m+delta_m, where delta_m is the number of new
      %          points sampled.
      
      %  4:1 Identify the set I of dimensions with the maximum side length.
      %      Let delta equal one-third of this maximum side length.

      [max_L, I] = max(L(:,j));
      I = find( L(:,j)==max_L );
      %NOT NEEDED I = find( abs( L(:,j) - max_L ) < tol1);    
      delta = c23*max_L;
      
      %  4:2 Sample the function at the points c +- delta*e_i for all
      %      i in I.
      %w=[];
      w=zeros(length(I),1);
      % HKH
      mz=m;
      if m+2*length(I) > size(C,2)
         C  = [C, zeros(n,MaxFunc)];
         F  = [F, zeros(1,MaxFunc)];
         L  = [L, zeros(n,MaxFunc)];
      end
%fprintf('size C %d %d, size I %d\n',size(C),length(I));
%fprintf('size F %d %d, size I %d\n',size(F),length(I));
      %NEWc_m = C(:,j);    % Centerpoint for new rectangle

      for ii = 1:length(I) % for each dimension with maximum side length
         i = I(ii);
         %e_i = [zeros(i-1,1);1;zeros(n-i,1)];
         %c_m1 = C(:,j) + delta*e_i;       % Centerpoint for new rectangle
         c_m1 = C(:,j);    % Centerpoint for new rectangle
         c_m1(i)=c_m1(i)+delta;

         %NEWc_m(i)=c_m(i)+delta;

         % Transform c_m1 to original search space
         %x_m1 = x_L + c_m1.*x_D;  

         x_m1 = tomsol(9,x_L, c_m1,x_D); 
         f_m1 = nlp_f(x_m1, Prob, varargin{:});   % Function value at x_m1
         if isnan(f_m1), f_m1 = realmax; end

         %c_m2 = C(:,j) - delta*e_i;    % Centerpoint for new rectangle
         c_m2 = C(:,j);                 % Centerpoint for new rectangle
         c_m2(i)=c_m2(i)-delta;
         x_m2 = tomsol(9,x_L, c_m2,x_D); % Transform c_m2 to original search space
         %x_m2 = x_L + c_m2.*x_D; % Transform c_m2 to original search space
         f_m2 = nlp_f(x_m2, Prob, varargin{:});   % Function value at x_m2
         if isnan(f_m2), f_m2 = realmax; end

         nFunc=nFunc+2;
         
         w(ii) = min(f_m1,f_m2);

         if convFlag == 0
            convFlag = isClose(fGoal,w(ii),fTol,nFunc,Iter,EpsGlob,PriLev);
         end
         
         %C = [C c_m1 c_m2];  % Matrix with all rectangle centerpoints
         %F = [F f_m1 f_m2];  % Vector with function values
         mz = mz+2;
         C(:,mz-1:mz) = [c_m1 c_m2];   % Matrix with all rectangle centerpoints
         F(:,mz-1:mz) = [f_m1 f_m2];   % Vector with function values
         
      end
      
      %  4:3 Divide the rectangle containing C(:,j) into thirds along the
      %      dimension in I, starting with the dimension with the lowest
      %      value of w(ii)
      [a b] = sort(w);
      glbfMin = min(glbfMin,a(1)); % Best new value in a(1)
%L = [L, zeros(size(L,1),2*length(I))];
%fprintf('size L %d %d, size I %d\n',size(L),length(I));
      %ix1 = m;
      for ii = 1:length(I)
         i = I(b(ii));
         
         ix1 = m + 2*b(ii)-1; % Index for new rectangle
         %ix2 = m + 2*b(ii);   % Index for new rectangle
         
         ix2 = ix1 + 1;
         
         L(i,j) = delta/2;
         
         %L(:,ix1) = L(:,j);
         %L(:,ix2) = L(:,j);
         L(:,[ix1,ix2]) = L(:,[j,j]);
         
         %D(j)   = sqrt(sum(L(:,j).^2));
         %D(ix1) = D(j);
         %D(ix2) = D(j);
         D([j ix1 ix2])   = norm(L(:,j));
         %zz   = norm(L(:,j));
         %D(j)   = zz;
         %D(ix1) = zz;
         %D(ix2) = zz;
         
      end
      m = m + 2*length(I);
   end
   
   % UPDATE:
   %glbfMin = min(F(1:m));  

   E = max(EpsGlob*abs(glbfMin),1E-8);
   DSort = D(1:m); % Use DSort instead of D before sorting in min computation
   [dummy iMin] = min( (F(1:m) - glbfMin + E)./DSort );

% [dummy iMinold] = min( (F - glbfMin + EpsGlob)./D );
% fprintf('\n iMin = %d     iMinold = %d',iMin,iMinold);
   
   i = 1;
   while 1
      d_tmp = DSort(i);
      idx = find(DSort~=d_tmp);
      DSort = [d_tmp DSort(idx)];
      if i==length(DSort)
         break;
      else
         i = i + 1;
      end
   end
   DSort = sort(DSort);
   
   DMin = zeros(1,length(DSort));
   for i = 1:length(DSort);
      idx1 = find(D(1:m)==DSort(i));
      %idx1 = find( abs( D-DSort(i) ) <= tol1 );
      DMin(i) = min(F(idx1));
   end
   
   if PriLev > 1 | IterPrint
      fprintf('\n Iteration: %d   glbfMin: %15.10f   Sampled points: %d',...
         Iter,glbfMin,nFunc);
   end
   
   Iter = Iter + 1;   
end % ITERATION LOOP

if PriLev > 1 | IterPrint
   fprintf('\n');
end


% SAVE RESULTS

%best=find(F==glbfMin)

C = C(:,1:m);    % All lengths
F = F(1:m);      % All function values computed
D = D(1:m);      % All distances
L = L(:,1:m);    % All lengths

WarmStartInfo.Name    = Name;
WarmStartInfo.C       = C;
WarmStartInfo.F       = F;
WarmStartInfo.D       = D;
WarmStartInfo.L       = L;
WarmStartInfo.DSort   = DSort;
WarmStartInfo.DMin    = DMin;
WarmStartInfo.glbfMin = glbfMin;      
WarmStartInfo.E       = E;
WarmStartInfo.iMin    = iMin;

SaveWarmStartFile(WarmStartInfo);
Result.glbSolve.WarmStartInfo = WarmStartInfo;

% All points i with F(i)=glbfMin
Result.x_k      = tomsol(9,x_L,C(:,find(F==glbfMin)),x_D);    

% One column for each optimal solution x

Result.f_k      = glbfMin;  % Best function value
Result.Iter     = Iter;     % Number of iterations
Result.FuncEv   = nFunc;
Result.ExitFlag = 0;
Result.ExitText = ['Tried ' num2str(m) ' function values in total '];
if cpumax
   Result.Inform   = 9;
   Result.ExitText = [Result.ExitText '. Max CPU reached. '];
else
   Result.Inform   = convFlag;
end
Result.maxTri   = max(L(:,m));
Result          = endSolve(Prob,Result);

function convFlag = isClose(fGoal,f,fTol,nFunc,Iter,EpsGlob,PriLev)

convFlag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convFlag = 1;
elseif fGoal == 0
   if abs(f-fGoal) < fTol
      convFlag = 2;
   end
elseif abs(f-fGoal) <= abs(fGoal) * fTol
   convFlag = 3;
end

if convFlag > 0 & PriLev >= 0 
   if convFlag == 1
      fprintf('\n\nFunction value %f is less than fGoal %f \n',f,fGoal);
   elseif convFlag == 2
      fprintf('\n\nError in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal),fTol);
   elseif convFlag == 3
      fprintf('\n\nRelative error in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal)/abs(fGoal),fTol);
   end
   fprintf('Number of function evaluations:  %d\n',nFunc);
   fprintf('Number of iterations:            %d\n',Iter);
   fprintf('Epsilon:                         %f\n',EpsGlob);
end
if convFlag == 3, convFlag = 2; end



% Save warmstart information to glbSave.mat - using try-catch to not loose
% information if the save fails. 
function SaveWarmStartFile(WarmStartInfo)

Name    = WarmStartInfo.Name;
C       = WarmStartInfo.C;
F       = WarmStartInfo.F;
D       = WarmStartInfo.D;
L       = WarmStartInfo.L;
DSort   = WarmStartInfo.DSort;
DMin    = WarmStartInfo.DMin;
glbfMin = WarmStartInfo.glbfMin;      
E       = WarmStartInfo.E;
iMin    = WarmStartInfo.iMin;

try
   save('glbSave.mat','Name','C','F','D','L','DSort','DMin','glbfMin','E','iMin');
catch
   warning('Failed to save warmstart information to glbSave.mat');
   disp(lasterr);
   disp('Warmstart information is available in Result.glbSolve.WarmStartInfo');
end


% MODIFICATION LOG
%
% 980915  mbk  Convergence test if known global optima are available.
% 980922  hkh  Change to use structure optParam. Use nFunc to count f calls.
% 980924  mbk  New definition of iMin:
%              [dummy iMin] = min( (F - glbfMin + E)./D ) instead of
%              [dummy iMin] = min( (F - glbfMin + EpsGlob)./D ).
% 980927  mbk  Wrong number of iterations performed.
% 981005  mbk  Changes in comment rows.
% 981013  hkh  Added call to iniSolve and endSolve
% 981026  hkh  Use SolverAlgorithm
% 981128  hkh  Change some printing levels for output
% 991220  hkh  Improve finding of convex hull
% 000830  hkh  Speedups
% 000923  hkh  Revised and simplified input format, use mat-file to save results
% 000927  hkh  Add IterPrint
% 001007  hkh  Do conhull in mex
% 001011  hkh  Add safety with Name in glbSave.mat
% 010414  hkh  Define isClose, and target value convergence
% 010715  hkh  Name changes. Influences glbSave.mat load / save.
% 010726  hkh  Minor speedups.
% 011031  hkh  Added Result.maxTri as output. Clean up comments.
% 020110  hkh  Change name fMin to glbfMin. 
% 020506  hkh  Improving comments about DIRECT algorithm
% 040111  hkh  Change call to inisolve
% 040327  hkh  Revised warm start comments
% 040327  hkh  Add Inform value, revised ExitFlag values
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 041123  hkh  Change call to tomRun in help
% 051005  hkh  Safe handling of MaxFunc, new default
% 060814  med  FUNCS used for callbacks instead
% 061003  ango Safe try-catch saving to glcSave.mat
