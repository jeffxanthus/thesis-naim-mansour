function varargout = mipDemo(varargin)
% MIPDEMO M-file for mipDemo.fig
%      MIPDEMO, by itself, creates a new MIPDEMO or raises the existing
%      singleton*.
%
%      H = MIPDEMO returns the handle to a new MIPDEMO or the handle to
%      the existing singleton*.
%
%      MIPDEMO('Property','Value',...) creates a new MIPDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to mipDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MIPDEMO('CALLBACK') and MIPDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MIPDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mipDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 13:51:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mipDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @mipDemo_OutputFcn, ...
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


% --- Executes just before mipDemo is made visible.
function mipDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for mipDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mipDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mipDemo_OutputFcn(hObject, eventdata, handles)
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
        mip0
    case 3
        mip1
    case 4
        mip2
    case 5
        mip3
    case 6
        mip4
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function mipHelp_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mipHelp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function mipHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file mipDemo in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'the  TOMLAB Quick (TQ) Format to formulate and solve MIP problems.  ' ...
,' ' ...
,'Four different ways to solve knapsack problems using TOMLAB are shown:' ...
,'With a standard branch and bound, adding a knapsack heuristic, adding' ...
,'priorities on the variables, or both approaches. ' ...
,' ' ...
,'Functions mip1, mip2, mip3 and mip4 in this file implements the four ways' ...
,' ' ...
,'The TOMLAB Quick Format (TQ) is illustrated in the routine Weingartner,'...
,'at the end of this file.' ...
,' ' ...
,'After defining the TOMLAB problem structure called Prob, is easy to set'...
,'different variables in Prob before actually solving the problem.' ...
,' ' ...
,'The call to a TOMLAB solver is always very simple, in this case'...
,'the call to mipSolve is just:  Result = mipSolve(Prob);' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function mip0
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Run very simple IP defined in the TOMLAB Quick format\n');
fprintf('=====================================================\n');

Name = 'Simple mip example'

A=[8 5 1 0
  -1 2 0 1 ];

disp('The matrix A')
A

disp('The right hand side b')
b=[31 6]';       % right-hand side

b

disp('The objective function coefficients c')
c=[-30  -18  0 0]'; % cost vector 

c

disp('Set lower bounds for x as zero and upper bounds infinite')

[m,n] = size(A);

x_L = zeros(n,1);       % NOTE! Must be set as 0, otherwise assumed -Inf
x_U = [];               % Default set as infinite

x_0 = [];

% All variables should be integer
IntVars=[1:n];    

b_U   = b;         % Upper bounds set as b
b_L   = b;         % Lower bounds equal to upper bounds, now when slacks
                   % have been added

x_min = [];        % Bounds for plotting, normally not used for MIP
x_max = [];        % Bounds for plotting, normally not used for MIP

f_opt = [];        % The function value at the optimal point
                   % Used in the printout. Assumed not known
x_opt = [];        % The optimal integer solution is assumed not to be known

f_Low = [];        % f_Low <= f_optimal must hold
                   % Set as -realmax, lowest possible number

fIP   = [];        % Do not use any prior knowledge about function values
                   % at some point where the solution is integer
xIP   = [];        % Do not use any prior knowledge about integer points

setupFile = []; % Just define the Prob structure, not any permanent Init File,
                % i.e. do not create a file according to the
                % TOMLAB Init File (IF) format

nProblem  = []; % Number of problems in the IF format is 0

VarWeight = []; % No variable priorities, largest fractional part will be used
Knapsack  = []; % No knapsack problem.

% Use the TOMLAB Quick (TQ) format
%
% Call the TOMLAB mipAssign routine that creates a structure with all
% problem information.

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 IntVars, VarWeight, Knapsack, fIP, xIP, ...
                 f_Low, x_min, x_max, f_opt, x_opt);

format compact
% Call the TOMLAB MIP solver directly
Result = mipSolve(Prob);

PrintResult(Result,2);

% ---------------------------------------------------------------------
function mip1
% ---------------------------------------------------------------------

% Use no knapsack heuristic, and no priorities

Knapsack  = [];

% Define the Weingartner 1 - 2/28 0-1 knapsack problem

Prob = Weingartner(Knapsack);

% Depth First, then Breadth search (Default is Depth First)
Prob.Solver.Alg=2;

% Must increase number of iterations from default 500
Prob.optParam.MaxIter=5000;

fprintf('No priority weighting used\n');

% Avoid the one line print each iteration
Prob.optParam.IterPrint = 0;

% Call the TOMLAB MIP solver directly
Result = mipSolve(Prob);

PrintResult(Result,2);

% ---------------------------------------------------------------------
function mip2
% ---------------------------------------------------------------------

% Add priorities on the variables

Knapsack  = 0;
'mip2'

% Define the Weingartner 1 - 2/28 0-1 knapsack problem
Prob = Weingartner(Knapsack);

% Set the priority weights equal to the coefficients in the objective
% Lower values give higher priority
Prob.MIP.VarWeight = Prob.QP.c;

fprintf('Using priority weighting on variables\n');

% Avoid the one line print each iteration
Prob.optParam.IterPrint = 0;

% Call the TOMLAB MIP solver directly
Result = mipSolve(Prob);

PrintResult(Result,2);

% ---------------------------------------------------------------------
function mip3
% ---------------------------------------------------------------------

% Use the knapsack heuristic, but not priorities

Knapsack  = 1;

% Define the Weingartner 1 - 2/28 0-1 knapsack problem
Prob = Weingartner(Knapsack);

% Depth First, then Breadth search (Default is Depth First)
Prob.Solver.Alg=2;

% Must increase number of iterations from default 500
Prob.optParam.MaxIter=5000;

fprintf('No priority weighting used\n');

% Avoid the one line print each iteration
Prob.optParam.IterPrint = 0;

% Call the TOMLAB MIP solver directly
Result = mipSolve(Prob);

PrintResult(Result,2);

% ---------------------------------------------------------------------
function mip4
% ---------------------------------------------------------------------

% Now use both the knapsack heuristic and priorities

Knapsack  = 1;

% Define the Weingartner 1 - 2/28 0-1 knapsack problem
Prob = Weingartner(Knapsack);

% Depth First, then Breadth search (Default is Depth First)
Prob.Solver.Alg=2;

% Also print one line / iteration
Prob.optParam.IterPrint=1;

% Set the priority weights equal to the coefficients in the objective
% Lower values give higher priority
Prob.MIP.VarWeight = Prob.QP.c;

fprintf('Using priority weighting on variables\n');

% Avoid the one line print each iteration
Prob.optParam.IterPrint = 0;

% Call the TOMLAB MIP solver directly
Result = mipSolve(Prob);

PrintResult(Result,2);



function Prob = Weingartner(Knapsack)

% Define the problem Weingartner, using the TOMLAB Quick (TQ) format

if nargin < 1, Knapsack = []; end

Name='Weingartner 1 - 2/28 0-1 knapsack';

% Problem formulated as a minimum problem

A = [ 45      0     85     150     65     95     30      0    170  0 ...
      40     25     20       0      0     25      0      0     25  0 ...
      165     0     85       0      0      0      0    100 ; ...
      30     20    125       5     80     25     35     73     12  15 ...
      15     40      5      10     10     12     10      9      0  20 ...
      60     40     50      36     49     40     19    150]; 

b_U = [600;600];  % // 2 knapsack capacities

c   = [1898  440  22507  270  14148   3100   4650  30800   615   4975 ...
  1160   4225    510   11880    479    440    490    330    110    560 ...
  24355   2885  11748    4550    750   3720   1950  10500]';  % 28 weights

% Make problem on standard form for mipSolve

[m,n]=size(A);

A=[A eye(m,m)]; % Add an identity matrix for the slack variables

% Change sign to make minimum problem
c=[-c;zeros(m,1)];

[mm nn]=size(A);
x_L=zeros(nn,1);
x_U=[ones(n,1);b_U];
x_0=[zeros(n,1);b_U];

fprintf('Knapsack problem. Variables %d. Knapsacks %d\n',n,m);

fprintf('Making standard form with   %d variables\n',nn);

if ~isempty(Knapsack) & Knapsack
   fprintf('Using knapsack heuristic\n');
else
   fprintf('No knapsack heuristic used\n');
end

% All original variables should be integer, but also slacks in this case
IntVars=[1:nn]; 

n   = length(c);
b_L = b_U;         % Lower bounds equal to upper bounds, now when slacks
                   % have been added

x_min=x_L;         % Bounds for plotting, normally not used for MIP
x_max=x_U;         % Bounds for plotting, normally not used for MIP

f_opt=-141278;     % The function value at the optimal point
                   % Used in the printout

f_Low=-1E7;        % f_Low <= f_optimal must hold
                   % Could be set as -realmax, lowest possible number
                   % This is default if set as empty, f_Low = [];

% This problem is number 7 in the tomlab\testprob\mip_prob.m file

fIP       = []; % Do not use any prior knowledge about function values
                % at some point where the solution is integer
xIP       = []; % Do not use any prior knowledge about integer points

setupFile = []; % Just define the Prob structure, not any permanent Init File,
                % i.e. do not create a file according to the
                % TOMLAB Init File (IF) format
nProblem  = []; % Number of problems in the IF format is 0

x_opt     = []; % The optimal integer solution is assumed not to be known
VarWeight = []; % No variable priorities, largest fractional part will be used

% Use the TOMLAB Quick (TQ) format
%
% Call the TOMLAB mipAssign routine that creates a structure with all
% problem information.

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 IntVars, VarWeight, Knapsack, fIP, xIP, ...
                 f_Low, x_min, x_max, f_opt, x_opt);
