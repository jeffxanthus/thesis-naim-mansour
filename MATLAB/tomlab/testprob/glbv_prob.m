% Defines variable dimension global optimization problems with simple bounds.
%
% glbv_prob: 
%
% function [probList, Prob] = glbv_prob(P, n, varargin)
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%    n      Dimension of problem
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2008-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 24, 2008.   Last modified Jun 24, 2008.

function [probList, Prob] = glbv_prob(P, n, varargin)

if nargin < 2
   n=[];
   if nargin < 1
       P=[];
end, end

probList=str2mat(...
     'Sphere model'...
    ,'Griewanks function'...
    ,'Shekels foxholes (n=2-10)'...
    ,'Michalewiczs function'...
    ,'Langermans function (n=2-10)'...
    ,'Easom (ES)'...
    ,'De Joung (DJ) '...
    ,'Rosenbrock'...
    ,'Zakharov'...
    ,'Griewank (GR) '...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return;
end

IntVars = []; uP = []; fGoal = []; nLocal = [];
nGlobal = []; f_Low = []; x_opt = []; f_opt = [];
x_0 = [];
b_L=[]; b_U=[]; A=[];

if n < 1, n = []; end

if  P == 1
    if isempty(n), n = 2; end
    Name = ['Sphere model ',num2str(n)];
    fGoal = 1e-6;
    f_Low = 0;
    f_opt = 0;
    x_opt = ones(1,n);
    x_L = -5*ones(n,1);
    x_U =  5*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 2
    if isempty(n), n = 2; end
    Name = ['Griewanks function ',num2str(n)];
    fGoal = 1e-4;
    f_Low = 0;
    f_opt = 0;
    x_opt = 100*ones(1,n);
    x_L   = -600*ones(n,1);
    x_U   =  600*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 3
    if isempty(n), n = 2; end
    n     = min(n,10);
    Name  = ['Shekels foxholes ',num2str(n)];
    x_opt = [];
    f_opt = [];
    if n == 2
       fGoal = -12.1190083798;
       f_opt = fGoal; % Guess
       x_opt = [8.024065,9.146534]; % from glbSolve / npsol
    elseif n == 5
       % from glbSolve / npsol
       fGoal = -10.4039206;
       f_opt = fGoal; % Guess
       x_opt = [8.024917,9.151728,5.113927,7.620861,4.564085];
    elseif n == 10
       %fGoal = -9; % Really OK?
       fGoal = -1.47736831291; % Obtained from glbSolve + npsol
       f_opt = fGoal; % Guess
       x_opt = [8.628544,4.409955,4.831266,5.767775,7.047738,6.712125, ...
                1.715764,4.323849, 4.405513,4.591235];
    end
    f_Low = -20;
    x_L   = zeros(n,1);
    x_U   = 10*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 4 
    if isempty(n), n = 2; end
    Name  = ['Michalewiczs function ',num2str(n)];
    if n == 2
       fGoal = -1.80130341008983;
       f_opt = fGoal; % Guess
       x_opt = [2.202905,1.570796]; % from glbSolve / npsol
    elseif n == 5
       fGoal = -4.687658179;  % from glbSolve + npsol
       f_opt = fGoal; % Guess
       x_opt = [2.202906,1.570796,1.284992,1.923058,1.720470];
    elseif n == 10
       fGoal = -7.84595327227;  % from glbSolve + npsol
       f_opt = fGoal; % Guess
       x_opt = [2.202906,1.570796,1.284992,2.934435,2.627180,1.570976, ...
                1.877235,2.604394 1.958965,1.570796];
    end
    x_L   = zeros(n,1);
    x_U   = pi*ones(n,1);
    f_Low = -20;
    x_min = x_L;
    x_max = x_U;
elseif P == 5
    if isempty(n), n = 2; end
    n     = min(n,10);
    Name  = ['Langermans function ',num2str(n)];
    f_opt = [];
    x_opt = [];
    if n == 2
       fGoal = -3.0677475617235372; % Accuracy from snopt
       f_opt = fGoal; % Guess
       x_opt = [7.735257,8.860257];
    elseif n == 5
       Name  = 'Langermans function 5';
       fGoal = -1.499943823715334600; % From snopt
       f_opt = fGoal; % Guess
       x_opt = [ 8.02329762255323 9.15653111755802 5.11632952767342 ...
                 7.61983157672815 4.56378371087901];
    elseif n == 10
       Name  = 'Langermans function 10';
       fGoal = -0.96399998152060;
       f_opt = fGoal; % Guess
       x_opt = [ 6.30598735426646 8.58299141102869 6.08399254062668 ...
                 1.13800014809038 4.35000836618424 3.13399328043645 ...
                 7.85300210477163 6.06100677670264 7.45700606167984 ...
                 2.25799751564348 ];
    end
    f_Low = -5;
    x_L   = zeros(n,1);
    x_U   = 10*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 6 % Easom (ES) 
    if isempty(n), n = 2; end
    Name  = ['Easom ',num2str(n)];
    f_opt = -1;
    x_opt = [];
    Name  = 'Easom (ES) ';
    x_L   = -100*ones(n,1);
    x_U   = 100*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 7	% De Joung (DJ) 3
    if isempty(n), n = 2; end
    Name  = ['De Joung (DJ) ',num2str(n)];
    f_opt = 0;
    x_opt = zeros(1,n);
    x_L   = -5.12*ones(n,1);
    x_U   = 5.12*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 8
    if isempty(n), n = 2; end
    Name  = ['Rosenbrock ',num2str(n)];
    x_L   = -5*ones(n,1);
    x_U   = 10*ones(n,1);
    x_opt = ones(1,n);
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 9
    if isempty(n), n = 2; end
    Name  = ['Zakharov ',num2str(n)];
    x_L   = -5*ones(n,1);
    x_U   = 10*ones(n,1);
    x_opt = zeros(1,n);
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 10	% Griewank (GR) 6
    if isempty(n), n = 6; end
    Name  = ['Griewank ',num2str(n)];
    x_L = -1*ones(n,1);
    x_U = -x_L;
    x_opt = zeros(1,n);
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
else
    error('glb_prob: Illegal problem number');
end
% Set x_0 to zeros (dummy for GUI)
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

% Define the Prob
Prob = glcAssign('glbv_f', x_L, x_U, Name, A, b_L, b_U, ...
    [], [], [], x_0, IntVars, [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.uP = uP;
Prob.MIP.nLocal  = nLocal;
Prob.MIP.nGlobal = nGlobal;
Prob.MIP.fGoal=fGoal(:);

% MODIFICATION LOG:
%
% 080624  hkh  Written, based on some problems from glb_prob
