function varargout = expDemo(varargin)
% EXPDEMO M-file for expDemo.fig
%      EXPDEMO, by itself, creates a new EXPDEMO or raises the existing
%      singleton*.
%
%      H = EXPDEMO returns the handle to a new EXPDEMO or the handle to
%      the existing singleton*.
%
%      EXPDEMO('Property','Value',...) creates a new EXPDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to expDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      EXPDEMO('CALLBACK') and EXPDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in EXPDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help expDemo

% Last Modified by GUIDE v2.5 27-Oct-2003 09:45:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @expDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @expDemo_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before expDemo is made visible.
function expDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for expDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes expDemo wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = expDemo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
val = get(hObject,'Value');
switch val
    case 2
        exp1Demo
    case 3
        exp2Demo
    case 4
        exp3Demo
    case 5
        exp4Demo
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --------------------------------------------------------------------
function expHelp_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
expHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function expHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file expDemo.m in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'the  TOMLAB Quick (TQ) Format to formulate and solve ls problems.  ' ...
,' ' ...
,'The TOMLAB Quick Format (TQ) is illustrated in the beginning of each of'...
,'the examples.' ...
,' ' ...
,'The code for the examples are all in this file (ls1Demo, ls2Demo etc.)' ...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'A real life parameter estimation problem with time values t and' ...
,'observation values Y(t) are used.' ...
,'The residual is computed in ls1_r and the Jacobian in ls1_J.' ...
,'The problem illustrates how to send extra information to the r and J' ...
,'routines, as well as communication between the r and J routine to save' ...
,'computational time. It is often the case that savings can be made using' ...
,'parts of the computation in r, especially exponential computations.' ...
,'Using the predefined global TOMLAB variable for the communication ' ...
,'between the r and J routine makes recursive calls possible because TOMLAB'...
,'can save the global variables before the recursive call.' ...
,' ' ...
,'The second example is the same problem solved with numerical Jacobian.' ...
,' ' ...
,'The third example shows how to send parameters down to the user routine'...
,'in a hidden way, in this case the weight parameter K. ' ...
,'The example also shows how to solve a sequence of similar problems,' ...
,'avoiding to redefine the problem structure in the TQ format each time.' ...
,' ' ...
,'Note that the call to a TOMLAB solver is always very simple, in this case'...
,'the call to clsSolve is just:  Result = clsSolve(Prob);' ...
,'We may call the TOMLAB driver routine instead. Then the call is'...
,'   Result = tomRun(''clsSolve'',Prob);  still very simple.' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function exp1Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Steyn and Wyk time series, 2+2 unknowns\n');
fprintf('=====================================================\n');

echo on

Name='Steyn and Wyk ';

% Time
t=[30:20:1170]';

% Observations
y=[18299 15428 13347 11466 10077 8729 7382 6708 5932 5352 4734 4271 ...
   3744 3485 3111 2950 2686 2476 2278 2107 1867 1768 1593 1504 1353 ...
   1389 1197 1080 1027 949 877 801 758 695 650 592 555 566 493 394 392 ...
   381 349 352 301 270 300 266 200 229 180 189 181 158 126 117 132 95]';

% Scale the problem
t=t/1000;     % Scale to seconds. Gives lambda*1000, of order 1
y=y/10000;    % Scale function values. Avoid large alpha

% Optimal values found:
x_opt=[3.831509;13.350236];  % Found by tests with wType=1 (1/y)
x_opt=[3.619756;11.548920;0.848085;1.507830];  % Found by tests, wType=0

% The routines exp_r and exp_J computes the exponential fitting residual
% and Jacobian for the given type of exponential model (eType)

% No special pattern is the Jacobian, it is a full matrix
JacPattern = [];

% Give number of terms, p
p = 2;

% Number of unknowns
n = 2*p;

% Just assign correct lengths for initial x_0, and bounds x_L and x_U
x_L = zeros(n,1);
x_U = ones(n,1);
x_0 = zeros(n,1);

% For exponential problem, use routine ExpFitW for weighting
weightType=3;
weightY='ExpFitW';  % Define function to compute weights

% Lower bound on optimal function value
fLowBnd = 0;

% SepAlg = 1 implies the use of a separable least squares algorithm
SepAlg  = 0; % Use ordinary least squares, no separation

% If linear constraints are present, set these in A, b_L and b_U
A = []; b_L = []; b_U = [];

% If nonlinear constraints are present, set these the bounds in c_L and c_U
c_L = []; c_U = [];
% The nonlinear constraint routine is c, contraint Jacobian is dc and
% the derivative pattern is set in ConsPattern
c = []; dc = []; ConsPattern = [];

f_opt =[];
x_min =x_L;
x_max =x_U;

Prob = clsAssign('exp_r', 'exp_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);

% Now reset the problem type, so Tomlab knows this is a exponential problem
% global probType
probType=checkType('exp');
Prob.probType=probType;     

% Set exponential parameters
lambda=[]; alpha=[]; beta=[];

% Select type of exponential model (see User's Guide)
eType=1;
% Select weighting type (see User's Guide)
wType=1; % Weight with data

% The following parameters are used by the initial value algorithm
% Normally do not change these values.

x0Type=0;
sumType=0;
infCR=0;
dType=0;
geoType=0;
qType=0;
sigType=0;

Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
                qType, sigType, lambda, alpha, beta, x0Type, sumType);

ask = 0; 

% If ask = 1, interactive setting of:	
% SepAlg (ordinary versus separable nonlinear least squares strategy)
% wType  (no weighting, or weighting with data)
% p      (number of exponential terms)
% All initial values for lambda, one by one. Suggestion given by initial
% value algorithm in expInit.

% Find initial values
Prob=expInit(Prob,ask);

x_0 = Prob.x_0;

% Add linear constraints to robustify the optimization, avoiding
% singularities because of exponential components being too close. 

Prob = ExpLinCon(Prob);


Result  = tomRun('clsSolve',Prob,2);

[TomV,os,TV] = tomlabVersion;

if TV(3)
   % License for NLSSOL is available
   disp(' ')
   disp('Run nlssol as well')
   disp(' ')
   disp('Press return to continue')
   pause
   % Here we may set parameters for the SOL solvers into the structure
   % See help nlssolTL for more help on the parameters
   Prob.SOL.optPar(1) = 11;  % Increase print level
   Prob.SOL.PrintFile = 'nlssol.txt';   % Name of text log print file
   Prob.SOL.SummFile  = 'nlssols.txt';  % Name of text log summary file

   % Special parameters for least squares problems are:
   % JTJ Initial Hessian (often best to have as true)
   Prob.SOL.optPar(47) = 1;  % Default 1, other unit Hessian

   % RESET Frequency . When Gauss-Newton works fine, often for small
   % residual problems, then one may raise this value
   Prob.SOL.optPar(48) = 2;  % Default 2, Reset each 2nd step

   Result  = tomRun('nlssol',Prob,2);
end

echo off

% ---------------------------------------------------------------------
function exp2Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Steyn and Wyk time series, Separable NLLS\n');
fprintf('=====================================================\n');

echo on

Name='Steyn and Wyk ';

% Time
t=[30:20:1170]';

% Observations
y=[18299 15428 13347 11466 10077 8729 7382 6708 5932 5352 4734 4271 ...
   3744 3485 3111 2950 2686 2476 2278 2107 1867 1768 1593 1504 1353 ...
   1389 1197 1080 1027 949 877 801 758 695 650 592 555 566 493 394 392 ...
   381 349 352 301 270 300 266 200 229 180 189 181 158 126 117 132 95]';

% Scale the problem
t=t/1000;     % Scale to seconds. Gives lambda*1000, of order 1
y=y/10000;    % Scale function values. Avoid large alpha

% Optimal values found:
x_opt=[3.831509;13.350236];  % Found by tests with wType=1 (1/y)
x_opt=[3.619756;11.548920;0.848085;1.507830];  % Found by tests, wType=0

% The routines exp_r and exp_J computes the exponential fitting residual
% and Jacobian for the given type of exponential model (eType)

% No special pattern is the Jacobian, it is a full matrix
JacPattern = [];

% Give number of terms, p
p = 2;

% Number of unknowns
n = 2*p;

% Just assign correct lengths for initial x_0, and bounds x_L and x_U
x_L = zeros(n,1);
x_U = ones(n,1);
x_0 = zeros(n,1);

% For exponential problem, use routine ExpFitW for weighting
weightType=3;
weightY='ExpFitW';  % Define function to compute weights

% Lower bound on optimal function value
fLowBnd = 0;

% SepAlg = 1 implies the use of a separable least squares algorithm
SepAlg  = 1; % Use ordinary least squares, no separation

% If linear constraints are present, set these in A, b_L and b_U
A = []; b_L = []; b_U = [];

% If nonlinear constraints are present, set these the bounds in c_L and c_U
c_L = []; c_U = [];
% The nonlinear constraint routine is c, contraint Jacobian is dc and
% the derivative pattern is set in ConsPattern
c = []; dc = []; ConsPattern = [];

f_opt =[];
x_min =x_L;
x_max =x_U;

Prob = clsAssign('exp_r', 'exp_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);

% Now reset the problem type, so Tomlab knows this is a exponential problem
% global probType
probType=checkType('exp');
Prob.probType=probType;     

% Set exponential parameters
lambda=[]; alpha=[]; beta=[];

% Select type of exponential model (see User's Guide)
eType=1;
% Select weighting type (see User's Guide)
wType=1; % Weight with data

% The following parameters are used by the initial value algorithm
% Normally do not change these values.

x0Type=0;
sumType=0;
infCR=0;
dType=0;
geoType=0;
qType=0;
sigType=0;

Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
                qType, sigType, lambda, alpha, beta, x0Type, sumType);

ask = 0; 

% If ask = 1, interactive setting of:	
% SepAlg (ordinary versus separable nonlinear least squares strategy)
% wType  (no weighting, or weighting with data)
% p      (number of exponential terms)
% All initial values for lambda, one by one. Suggestion given by initial
% value algorithm in expInit.

% Find initial values, expInit will set correct lengths of x_0, x_L, x_U
Prob=expInit(Prob,ask);

x_0 = Prob.x_0;

% Add linear constraints to robustify the optimization, avoiding
% singularities because of exponential components being too close. 

Prob = ExpLinCon(Prob);

Result  = tomRun('clsSolve',Prob,2);

[TomV,os,TV] = tomlabVersion;

if TV(3)
   % License for NLSSOL is available
   disp(' ')
   disp('Run nlssol as well')
   disp(' ')
   disp('Press return to continue')
   pause

   Result  = tomRun('nlssol',Prob,2);
end

echo off

% ---------------------------------------------------------------------
function exp3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Steyn and Wyk time series, Interactive user settings\n');
fprintf('=====================================================\n');

echo on

Name='Steyn and Wyk ';

% Time
t=[30:20:1170]';

% Observations
y=[18299 15428 13347 11466 10077 8729 7382 6708 5932 5352 4734 4271 ...
   3744 3485 3111 2950 2686 2476 2278 2107 1867 1768 1593 1504 1353 ...
   1389 1197 1080 1027 949 877 801 758 695 650 592 555 566 493 394 392 ...
   381 349 352 301 270 300 266 200 229 180 189 181 158 126 117 132 95]';

% Scale the problem
t=t/1000;     % Scale to seconds. Gives lambda*1000, of order 1
y=y/10000;    % Scale function values. Avoid large alpha

% Optimal values found:
x_opt=[3.831509;13.350236];  % Found by tests with wType=1 (1/y)
x_opt=[3.619756;11.548920;0.848085;1.507830];  % Found by tests, wType=0

% The routines exp_r and exp_J computes the exponential fitting residual
% and Jacobian for the given type of exponential model (eType)

% No special pattern is the Jacobian, it is a full matrix
JacPattern = [];

% Give number of terms, p
p = 2;

% Number of unknowns
n = 2*p;

% Just assign correct lengths for initial x_0, and bounds x_L and x_U
x_L = zeros(n,1);
x_U = ones(n,1);
x_0 = zeros(n,1);

% For exponential problem, use routine ExpFitW for weighting
weightType=3;
weightY='ExpFitW';  % Define function to compute weights

% Lower bound on optimal function value
fLowBnd = 0;

% SepAlg = 1 implies the use of a separable least squares algorithm
SepAlg  = 0; % Use ordinary least squares, no separation

% If linear constraints are present, set these in A, b_L and b_U
A = []; b_L = []; b_U = [];

% If nonlinear constraints are present, set these the bounds in c_L and c_U
c_L = []; c_U = [];
% The nonlinear constraint routine is c, contraint Jacobian is dc and
% the derivative pattern is set in ConsPattern
c = []; dc = []; ConsPattern = [];

f_opt =[];
x_min =x_L;
x_max =x_U;

Prob = clsAssign('exp_r', 'exp_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);

% Now reset the problem type, so Tomlab knows this is a exponential problem
% global probType
probType=checkType('exp');
Prob.probType=probType;     

% Set exponential parameters
lambda=[]; alpha=[]; beta=[];

% Select type of exponential model (see User's Guide)
eType=1;
% Select weighting type (see User's Guide)
wType=1; % Weight with data

% The following parameters are used by the initial value algorithm
% Normally do not change these values.

x0Type=0;
sumType=0;
infCR=0;
dType=0;
geoType=0;
qType=0;
sigType=0;

Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
                qType, sigType, lambda, alpha, beta, x0Type, sumType);

ask = 1; % NOTE - NOW THE USER WILL SET WHAT VALUES TO USE 

% If ask = 1, interactive setting of:	
% SepAlg (ordinary versus separable nonlinear least squares strategy)
% wType  (no weighting, or weighting with data)
% p      (number of exponential terms)
% All initial values for lambda, one by one. Suggestion given by initial
% value algorithm in expInit.

% Find initial values, expInit will set correct lengths of x_0, x_L, x_U
Prob=expInit(Prob,ask);

x_0 = Prob.x_0;

% Add linear constraints to robustify the optimization, avoiding
% singularities because of exponential components being too close. 

Prob = ExpLinCon(Prob);

Result  = tomRun('clsSolve',Prob,2);

Result  = tomRun('clsSolve',Prob,2);

[TomV,os,TV] = tomlabVersion;

if TV(3)
   disp(' ')
   disp('Run nlssol as well')
   disp(' ')
   disp('Press return to continue')
   pause

   Result  = tomRun('nlssol',Prob,2);
end

echo off

% ---------------------------------------------------------------------
function exp4Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Steyn and Wyk time series, 2+2 unknowns, IF format\n');
fprintf('=====================================================\n');

echo on

Name='Steyn and Wyk ';

% Steyn and Wyk problem is predefined as problem 17.
%
% The user may edit the file exp_prob.m in tomlab\testprob
% create his own version of exp_prob (as described in the User's Guide)
%
% When this editing is done for a new problem, it is easy to use the
% init file (IF) format:

Prob = probInit('exp_prob',17);

% If linear constraints are present, they should normally be set in
% the Init File. Otherwise, set them here into the structure. 
% After that, expInit MUST BE CALLED.
% Prob.A = []; Prob.b_L = []; Prob.b_U = [];

% If nonlinear constraints are present, their bounds should normally be set in
% the Init File. Otherwise, set them here into the structure. 
% Prob.c_L = []; Prob.c_U = [];
% The nonlinear constraint routine name c and contraint Jacobian routine dc
% must be set into the Prob structure, if there are nonlinear constraints.
% The value of Prob.FUNCS.c should be a string with the Matlab name without
% any extension
% Prob.FUNCS.c = []; Prob.FUNCS.dc = []; 


% Type of exponential model is normally set in the Init File (see User's Guide)
eType=1;
% Weighting type should be set in the Init File (see User's Guide)
wType=1; % Weight with data

% Changing to 1 term exponential
p = 1;

% SepAlg = 1 implies the use of a separable least squares algorithm
SepAlg  = 1; % Use separable nonlinear least squares

% Calling expSet1 sets the above parameters into the Prob structure
Prob = expSet1(Prob, SepAlg, p, wType, eType);

% For the initial value algorith, there is a similar routine expSet2, but
% these parameters are seldom changed.

% If the user now has changed the problem, e.g the number of terms, or
% separable or non-separable strategy, then expInit should be called again, 
% to estimate new initial values x_0, and bounds x_L and x_U.

ask = 0; 

% If ask = 1, interactive setting of:	
% SepAlg (ordinary versus separable nonlinear least squares strategy)
% wType  (no weighting, or weighting with data)
% p      (number of exponential terms)
% All initial values for lambda, one by one. Suggestion given by initial
% value algorithm in expInit.

% Find initial values
Prob=expInit(Prob,ask);

x_0 = Prob.x_0;

% Add linear constraints to robustify the optimization, avoiding
% singularities because of exponential components being too close. 
%
% Must be called again if p has been changed.

Prob = ExpLinCon(Prob); % In this case no linear constraints are created
                        % because the number of terms are 1 !!!

Result  = tomRun('clsSolve',Prob,2);

Result  = tomRun('clsSolve',Prob,2);

[TomV,os,TV] = tomlabVersion;

if TV(3)
   disp(' ')
   disp('Run nlssol as well')
   disp(' ')
   disp('Press return to continue')
   pause

   Result  = tomRun('nlssol',Prob,2);
end

echo off
