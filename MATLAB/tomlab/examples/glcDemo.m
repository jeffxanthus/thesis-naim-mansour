function varargout = glcDemo(varargin)
% GLCDEMO M-file for glcDemo.fig
%      GLCDEMO, by itself, creates a new GLCDEMO or raises the existing
%      singleton*.
%
%      H = GLCDEMO returns the handle to a new GLCDEMO or the handle to
%      the existing singleton*.
%
%      GLCDEMO('Property','Value',...) creates a new GLCDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to glcDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GLCDEMO('CALLBACK') and GLCDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GLCDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help glcDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 15:01:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @glcDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @glcDemo_OutputFcn, ...
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

% --- Executes just before glcDemo is made visible.
function glcDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for glcDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes glcDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = glcDemo_OutputFcn(hObject, eventdata, handles)
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
        glc1Demo
    case 3
        glc2Demo
    case 4
        glc3Demo
    case 5
        glc4Demo
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function glcHelp_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glcHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function glcHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file glcDemo.m in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'the TOMLAB Quick (TQ) and Init File (IF) format to solve glc problems  ' ...
,'(global box-bounded integer, linear and nonlinearly constrained problems).'...
,' ' ...
,'The code for the examples are all in this file (glc1Demo, glc2Demo etc.)'...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The Floudas-Pardalos 3.3 test problem is used as the main example' ...
,' ' ...
,'Note that the call to a TOMLAB solver is always very simple:'...
,'   Result = glcSolve(Prob);' ...
,' ' ...
,'We may call the TOMLAB driver routine instead. Then the call is'...
,'   Result = tomRun(''glcSolve'',Prob);' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function glc1Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Floudas-Pardalos 3.3 with TQ format, 10 restarts\n');
fprintf('=====================================================\n');

% Example of constrained global optimization problems.
%
% This example is number 16 in glc_prob.m

Name   = 'Floudas-Pardalos 3.3';

x_L   = [ 0  0 1 0 1 0]';  % Lower bounds on x
A =     [ 1 -3 0 0 0 0     
         -1  1 0 0 0 0
          1  1 0 0 0 0];   % Linear equations
b_L   = [-inf -inf  2 ]';  % Upper bounds for linear equations
b_U   = [  2    2   6 ]';  % Lower bounds for linear equations
%
% Original problem has no bounds on x(1) and x(2)
%
% x(1) >= 0, x(2) >= 0 & linear equation 3:  x(1) + x(2) <= 6  ==>
%
% x(1) <= 6 and x(2) <=6. This is inserted to get a box-bounded problem

x_U   = [6 6 5 6 5 10]';  % Upper bounds after x(1),x(2) values inserted
c_L   = [4 4]';           % Lower bounds on two nonlinear constraints
c_U   = [];               % Upper bounds are infinity for nonlinear constraints


x_opt = [5 1 5 0 5 10]'; 
f_opt = -310;
x_min = x_L;
x_max = x_U;

% Set the rest of the arguments as empty
IntVars   = []; VarWeight = [];
fIP       = []; xIP       = []; fLowBnd = [];
x_0       = [];

%IntVars = [1:5];

Prob = glcAssign('glc4_f', x_L, x_U, Name, ...
                 A, b_L, b_U, 'glc4_c', c_L, c_U, x_0, ... 
                 IntVars, VarWeight, fIP, xIP, ...
                 fLowBnd, x_min, x_max, f_opt, x_opt);

% Increase the default max number of function evaluations in glcSolve

Prob.optParam.MaxFunc = 500;

Result = glcSolve(Prob);

PrintResult(Result,3);

Prob.WarmStart = 1;

% Do 10 restarts, call tomRun, PriLev = 2 gives call to PrintResult
for i=1:10
    Result = tomRun('glcSolve',Prob,2);
end

% ---------------------------------------------------------------------
function glc2Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Floudas-Pardalos 3.3 with Init File format I\n');
fprintf('=====================================================\n');

% Generate the problem structure using the TOMLAB Init File format
%
% Floudas-Pardalos 3.3 is number 16 in glc_prob.m

Prob   = probInit('glc_prob', 16);

% Increase the default max number of function evaluations in glcSolve

Prob.optParam.MaxFunc = 2000;

Result = glcSolve(Prob);

PrintResult(Result);

% ---------------------------------------------------------------------
function glc3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Floudas-Pardalos 3.3 with Init File format II\n');
fprintf('=====================================================\n');

% Generate the problem structure using the TOMLAB Init File format


% Increase the default max number of function evaluations in glcSolve

Prob.optParam.MaxFunc = 2000;

% Use driver routine to pick up problem
% Floudas-Pardalos 3.3 is number 16 in glc_prob.m
% Increase print level to get a printing from PrintResult

Result = tomRun('glcSolve','glc_prob',16,Prob,2);


function glc4Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Sending information to f(x) and c(x) in problem HS332\n');
fprintf('=====================================================\n');

% This is problem 12 in glc_prob.m

Name = 'Hock - Schittkowski 332';
x_L = [0 0]'; 
x_U = [1.5 1.5]'; 
b_L = []; b_U = []; A = [];  % No linear constraints
c_L = -30;                   % Lower bound, one nonlinear constraint
c_U = 30;                    % Upper bound, one nonlinear constraint


% x_0 = [0.75 0.75]'; % If running local solver

% Generate the problem structure using the TOMLAB Quick format (short call)

Prob = glcAssign('glc5_f', x_L, x_U, Name, A, b_L, b_U, 'glc5_c', c_L, c_U); 

% Set optimal points directly in structure, shorter call to glcAssign then
Prob.x_opt = [0.9114 0.02928];
Prob.f_opt = 29.92437939227878;

% f_opt = 114.95; % This is f_opt given in the book "More testexamples .."
% Obviously wrong !!!


% Add information to be sent to glc5_f and glc5_c
tt             = [1:100]';
Prob.user.t    = pi*(1/3+(tt-1)/180);      % Used in glc5_f and glc5_c

% For illustration, use multi-solver driver routine tomRun, instead of
% direct call to glcSolve. Set print level to 2, gives call to PrintResult

Result = tomRun('glcSolve',Prob,2);
