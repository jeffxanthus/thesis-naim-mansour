function varargout = lsDemo(varargin)
% LSDEMO M-file for lsDemo.fig
%      LSDEMO, by itself, creates a new LSDEMO or raises the existing
%      singleton*.
%
%      H = LSDEMO returns the handle to a new LSDEMO or the handle to
%      the existing singleton*.
%
%      LSDEMO('Property','Value',...) creates a new LSDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to lsDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LSDEMO('CALLBACK') and LSDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LSDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lsDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 14:33:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lsDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @lsDemo_OutputFcn, ...
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


% --- Executes just before lsDemo is made visible.
function lsDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for lsDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lsDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lsDemo_OutputFcn(hObject, eventdata, handles)
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
        ls1Demo
    case 3
        ls2Demo
    case 4
        ls3Demo
    case 5
        ls4Demo
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --------------------------------------------------------------------
function lsHelp_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lsHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function lsHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file lsDemo.m in tomlab\examples ' ...
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
function ls1Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Nonlinear parameter estimation, 3 unknowns\n');
fprintf('=====================================================\n');

echo on

Name='Gisela';

% Time values
t  = [0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
      96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];

% Observations
y = [30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
      19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];

x_0 = [6.8729,0.0108,0.1248]'; % Initial values for unknown x


% Generate the problem structure using the TOMLAB Quick format (short call)
Prob = clsAssign('ls1_r', 'ls1_J', [], [], [], Name, x_0, y, t);

% See next example for a more general call to clsAssign

%Prob = clsAssign(r, J, JacPattern, x_L, x_U, Name, x_0, ...
%                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
%                 A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
%                 x_min, x_max, f_opt, x_opt);

% Generate the problem structure using the TOMLAB Quick format (short call)
%Prob    = probAssign('ls', x_L, x_U, Name, x_0, fLowBnd);

% Update the Prob structure with the names of files
%Prob    = tomFiles(Prob,[],[],[],[],[],[],'ls1_r', 'ls1_J');

%Prob.LS.y = y;
%Prob.LS.t  = t;

% Weighting parameter K in model is sent to r and J computation using Prob
Prob.userParam.K = 5;  

Result  = clsSolve(Prob);

PrintResult(Result,2);

echo off

% ---------------------------------------------------------------------
function ls2Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Nonlinear parameter estimation, 3 unknowns, numerical J\n');
fprintf('=====================================================\n');

echo on

Name='Gisela';

% Time values
t  = [0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
      96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];

% Observations
y = [30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
      19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];

x_0  = [6.8729,0.0108,0.1248]';
x_L        = [-Inf,-Inf,-Inf]';       % Could be set as empty, default is -Inf
x_U        = [Inf,Inf,Inf]';          % Could be set as empty, default is Inf
JacPattern = [];                      % No sparsity pattern set for Jacobian

%x_min=[6 0.005 0.1]';
%x_max=[7 0.05  0.2]';

fLowBnd = 0;         % Lower bound on function.

% Generate the problem structure using the TOMLAB Quick format (short call)
Prob = clsAssign('ls1_r', 'ls1_J', [], [], [], Name, x_0, y, t);

%Prob = clsAssign(r, J, JacPattern, x_L, x_U, Name, x_0, ...
%                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
%                 A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
%                 x_min, x_max, f_opt, x_opt);


% Weigthing parameter K in model is sent to r and J computation using Prob
Prob.userParam.K        = 5;  
Prob.NumDiff            = 1; % Use standard numerical differences
Prob.optParam.IterPrint = 1; % Print one line each iteration

Result  = tomRun('clsSolve',Prob,2);

echo off

% ---------------------------------------------------------------------
function ls3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Nonlinear parameter estimation, 3 unknowns, sequence of calls\n');
fprintf('=====================================================\n');

echo on

Name='Gisela';

% Time values
t  = [0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
      96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];

% Observations
y = [30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
      19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];

x_0  = [6.8729,0.0108,0.1248]';

% Generate the problem structure using the TOMLAB Quick format (short call)
Prob = clsAssign('ls1_r', 'ls1_J', [], [], [], Name, x_0, y, t);

% Weigthing parameter K in model is sent to r and J computation using Prob

for i=1:6

    Prob.userParam.K = 3.8 + 0.2*i;  

    Result(i)  = tomRun('clsSolve',Prob,2);

    fprintf('\nWEIGHT PARAMETER K is %9.3f\n\n\n',Prob.userParam.K);
    pause(3)
end

echo off

% ---------------------------------------------------------------------
function ls4Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Nonlinear parameter estimation, 3 unknowns\n');
fprintf('Calling general NLP solver\n');
fprintf('=====================================================\n');

echo on

Name='Gisela';

% Time values
t  = [0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
      96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];

% Observations
y = [30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
      19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];

x_0  = [6.8729,0.0108,0.1248]';

% Generate the problem structure using the TOMLAB Quick format (short call)
Prob = clsAssign('ls1_r', 'ls1_J', [], [], [], Name, x_0, y, t);

% Weigthing parameter K in model is sent to r and J computation using Prob
Prob.userParam.K = 5;  

% The use of a general solver is automatic, because TOMLAB is defining
% general routines to compute f(x) (ls_f.m ) g(x) (ls_g.m) and an
% approximation of the Hessian H(x) = J'(x)*J(x) as routine (ls_H.m)

Result  = tomRun('conSolve',Prob,2);

echo off
