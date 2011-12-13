function varargout = llsDemo(varargin)
% LLSDEMO M-file for llsDemo.fig
%      LLSDEMO, by itself, creates a new LLSDEMO or raises the existing
%      singleton*.
%
%      H = LLSDEMO returns the handle to a new LLSDEMO or the handle to
%      the existing singleton*.
%
%      LLSDEMO('Property','Value',...) creates a new LLSDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to llsDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LLSDEMO('CALLBACK') and LLSDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LLSDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help llsDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 15:07:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @llsDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @llsDemo_OutputFcn, ...
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


% --- Executes just before llsDemo is made visible.
function llsDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for llsDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes llsDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = llsDemo_OutputFcn(hObject, eventdata, handles)
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
        lls1Demo
    case 3
        lls2Demo
    case 4
        lls3Demo
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function llsHelp_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
llsHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function llsHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '--------------------------------------------- ' ...
,'This is the file llsDemo.m in tomlab\examples ' ...
,'--------------------------------------------- ' ...
,' ' ...
,'Here are shown different examples that illustrates the use of'...
,'Tomlab to formulate and solve linear least squares problems.  ' ...
,' ' ...
,'The code for the examples are all in this file (lls1Demo, lls2Demo etc.)' ...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The lls problem if formulated very general, with linear equalities and' ...
,'inequalities, and simple bounds. One sparse linear least squares solver,' ...
,'Tlsqr, only solves large and sparse unconstrained linear least squares'...
,'The LSEI solver solves the general problem. ' ...
,'The Tlsqr and LSEI solver are part of the Tomlab Base Module. ' ...
,' ' ...
,'If formulating the problem in the TQ or IF format, a general solver' ...
,'like clsSolve or conSolve (minos) could solve the problem as well' ...
,' ' ...
,'The first example shows a linear least squares problem with linear ' ...
,'inequalities and simple bounds' ...
,'The example is a test example in the LSSOL Users Guide'...
,' ' ...
,'In the first demo the Tomlab Quick (TQ) Format is used to get a' ...
,'problem structure and then lsei is called to solve the problem.' ...
,' ' ...
,'The second demo shows how to call LSSOL directly to solve the same'...
,'problem. This example only works for Tomlab/NPSOL or Tomlab/SOL,' ...
,'where the LSSOL solver is one of the solvers.' ...
,' ' ...
,'The third example is similar to the first example, but calls the'...
,'constrained nonlinear least squares solver clsSolve instead.'...
,'It is simple to call any other solver, it is just to change the name'...
,'in the call to the driver routine tomRun.'...
,' ' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function lls1Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Least squares example from LSSOL Users Guide\n');
fprintf('=====================================================\n');

echo on

% Test of a Least Squares example from LSSOL User's Guide

Name='LSSOL test example';

% In Tomlab it is best to use Inf and -Inf, not big numbers.
n = 9;  % Number of unknown parameters are 9
x_L = [-2 -2 -Inf, -2*ones(1,6)]';
x_U = 2*ones(n,1);

A   = [ ones(1,8) 4; 1:4,-2,1 1 1 1; 1 -1 1 -1, ones(1,5)];
b_L = [2    -Inf -4]';
b_U = [Inf    -2 -2]';

bl = [x_L;b_L];
bu = [x_U;b_U];

y = ones(10,1); % Number of observations are 10.
C = [ ones(1,n); 1 2 1 1 1 1 2 0 0; 1 1 3 1 1 1 -1 -1 -3; ...
      1 1 1 4 1 1 1 1 1;1 1 1 3 1 1 1 1 1;1 1 2 1 1 0 0 0 -1; ...
      1 1 1 1 0 1 1 1 1;1 1 1 0 1 1 1 1 1;1 1 0 1 1 1 2 2 3; ...
      1 0 1 1 1 1 0 2 2];

x_0 = 1./[1:n]';

% We may optionally set the optimal values, if known
% x_opt estimated from LSSOL, rather similar to User's Guide
x_opt = [2 1.57195927 -1.44540327 -0.03700275 0.54668583 0.17512363 ...
          -1.65670447 -0.39474418  0.31002899]; 
f_opt = 0.1390587318; % Estimated from LSSOL, wrong in User's Guide

% Use the llsAssign routine to make a Prob structure

t          = [];   % No time set for y(t) (used for plotting)
weightY    = [];   % No weighting
weightType = [];   % No weighting type set
x_min      = [];   % No lower bound for plotting
x_max      = [];   % No upper bound for plotting

Prob = llsAssign(C, y, x_L, x_U, Name, x_0, t, weightType, weightY, ...
                 A, b_L, b_U,  x_min, x_max, f_opt, x_opt);

Result  = tomRun('lsei',Prob,2);

%PrintResult(Result,2);

% ---------------------------------------------------------------------
function lls2Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Least squares example from LSSOL Users Guide\n');
[TomV,os,TV] = tomlabVersion;

if TV(3)
   fprintf('Solve by a direct call to LSSOL\n');
else
   fprintf('No license for the LSSOL solver!\n');
   return
end
fprintf('=====================================================\n');

% Test of a Least Squares example from LSSOL User's Guide

% Direct call to LSSOL. 

% Note that when calling the LSSOL MEX interface directly, avoid using
% Inf and -Inf. Instead use big numbers that indicate Inf.
% The standard for the MEX interfaces is 1E20 and -1E20, respectively.

n = 9; % There are nine unknown parameters, and 10 equations
x_L = [-2 -2 -1E20, -2*ones(1,6)]';
x_U = 2*ones(9,1);

A   = [ ones(1,8) 4; 1:4,-2,1 1 1 1; 1 -1 1 -1, ones(1,5)];
b_L = [2    -1E20 -4]';
b_U = [1E20    -2 -2]';

% Must put lower and upper bounds on variables and constraints together
bl = [x_L;b_L];  
bu = [x_U;b_U];

y = ones(10,1);
H = [ ones(1,n); 1 2 1 1 1 1 2 0 0; 1 1 3 1 1 1 -1 -1 -3; ...
      1 1 1 4 1 1 1 1 1;1 1 1 3 1 1 1 1 1;1 1 2 1 1 0 0 0 -1; ...
      1 1 1 1 0 1 1 1 1;1 1 1 0 1 1 1 1 1;1 1 0 1 1 1 2 2 3; ...
      1 0 1 1 1 1 0 2 2];

x_0 = 1./[1:n]';

% Set empty indicating default values for most variables
c          = [];          % No linear coefficients, they are for LP/QP
Warm       = [];          % No warm start
iState     = [];          % No warm start
Upper      = [];          % C is not factorized
kx         = [];          % No warm start
SpecsFile  = [];          % No parameter settings in a SPECS file
PriLev     = [];          % PriLev is not really used in LSSOL
ProbName   = [];          % ProbName is not really used in LSSOL

optPar(1)  = 50;          % Set print level at maximum 
PrintFile  = 'lssol.txt'; % Print result on the file with name lssol.txt

z0 = (y-H*x_0);

f0 = 0.5*z0'*z0;
fprintf('Initial function value %f\n',f0);

[x, Inform, iState, cLamda, Iter, fObj, r, kx] = ...
    lssol( A, bl, bu, c, x_0, optPar, H, y, Warm, ...
          iState, Upper, kx, SpecsFile, PrintFile, PriLev, ProbName );

% We could equally well call with the following shorter call:
% [x, Inform, iState, cLamda, Iter, fObj, r, kx] = ...
%     lssol( A, bl, bu, c, x, optPar, H, y); 

z = (y-H*x);

f = 0.5*z'*z;
fprintf('Optimal function value %f\n',f);

fprintf('Optimal parameters:\n');
xprint(x,'x:');

echo off

% ---------------------------------------------------------------------
function lls3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Least squares example from LSSOL Users Guide\n');
fprintf('Call genereal cls solver clsSolve\n');
fprintf('=====================================================\n');

echo on

% Test of a Least Squares example from LSSOL User's Guide

Name='LSSOL test example';

% In Tomlab it is best to use Inf and -Inf, not big numbers.
n = 9;
x_L = [-2 -2 -Inf, -2*ones(1,6)]';
x_U = 2*ones(9,1);

A   = [ ones(1,8) 4; 1:4,-2,1 1 1 1; 1 -1 1 -1, ones(1,5)];
b_L = [2    -Inf -4]';
b_U = [Inf    -2 -2]';

bl = [x_L;b_L];
bu = [x_U;b_U];

y = ones(10,1);
C = [ ones(1,n); 1 2 1 1 1 1 2 0 0; 1 1 3 1 1 1 -1 -1 -3; ...
      1 1 1 4 1 1 1 1 1;1 1 1 3 1 1 1 1 1;1 1 2 1 1 0 0 0 -1; ...
      1 1 1 1 0 1 1 1 1;1 1 1 0 1 1 1 1 1;1 1 0 1 1 1 2 2 3; ...
      1 0 1 1 1 1 0 2 2];

x_0 = 1./[1:n]';

% We may optionally set the optimal values, if known
% x_opt estimated from LSSOL, rather similar to User's Guide
x_opt = [2 1.57195927 -1.44540327 -0.03700275 0.54668583 0.17512363 ...
          -1.65670447 -0.39474418  0.31002899]; 
f_opt = 0.1390587318; % Estimated from LSSOL, wrong in User's Guide

% Use the llsAssign routine to make a Prob structure

t          = [];   % No time set for y(t) (used for plotting)
weightY    = [];   % No weighting
weightType = [];   % No weighting type set
x_min      = [];   % No lower bound for plotting
x_max      = [];   % No upper bound for plotting

Prob = llsAssign(C, y, x_L, x_U, Name, x_0, t, weightType, weightY, ...
                 A, b_L, b_U,  x_min, x_max, f_opt, x_opt);

% Now use the general nonlinear solver clsSolve to solve the problem

Result  = tomRun('clsSolve',Prob,2);
