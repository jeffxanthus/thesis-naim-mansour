function varargout = ucDemo(varargin)
% UCDEMO M-file for ucDemo.fig
%      UCDEMO, by itself, creates a new UCDEMO or raises the existing
%      singleton*.
%
%      H = UCDEMO returns the handle to a new UCDEMO or the handle to
%      the existing singleton*.
%
%      UCDEMO('Property','Value',...) creates a new UCDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ucDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      UCDEMO('CALLBACK') and UCDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in UCDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ucDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 14:12:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ucDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @ucDemo_OutputFcn, ...
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


% --- Executes just before ucDemo is made visible.
function ucDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ucDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ucDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ucDemo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function ucHelp_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ucHelp;

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
        uc1Demo
    case 3
        uc2Demo
    case 4
        uc3Demo
    case 5
        uc4Demo
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function ucHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '-------------------------------------------- ' ...
,'This is the file ucDemo.m in tomlab\examples ' ...
,'-------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'the TOMLAB Format to formulate and solve uc problems.  ' ...
,' ' ...
,'The TOMLAB Format is illustrated in the beginning of each of'...
,'the examples.' ...
,' ' ...
,'The code for the examples are all in this file (uc1Demo, uc2Demo etc.)' ...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The solution of the Rosenbrocks banana problem is used as main example' ...
,'First in the case where the function, the gradient and the Hessian is ' ...
,'written in the routines uc1_f.m, uc1_g.m and uc1_H.m.' ...
,'Then Newtons method can be applied, and the convergence is rapid.' ...
,' ' ...
,'The second example shows how to solve the problem when the user is only' ...
,'able to give the function, not any gradients or the Hessian' ...
,' ' ...
,'The third example shows how to send parameters to the user routine'...
,'in a hidden way, in this case a steepness parameter in the more general' ...
,'formulation of the Rosenbrocks banana function.' ...
,'The example also shows how to solve a sequence of similar problems,' ...
,'avoiding to redefine the problem structure in the TQ format each time.' ...
,' ' ...
,'The preferred way to call any TOMLAB solver is through the tomRun' ...
,'driver routine: ' ...
,'   Result = tomRun(''ucSolve'',Prob,1); ' ...
,' ' ...
,'To use a general NLP solver, just change the solver name:' ...
,'   Result = tomRun(''npsol'',Prob,1);' ...
,'   Result = tomRun(''snopt'',Prob,1);' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function uc1Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Rosenbrocks banana with explicit f(x), g(x) and H(x)\n');
fprintf('=====================================================\n');

echo on

Name    = 'RB Banana';
x_0     = [-1.2 1]';   % Starting values for the optimization.
x_L     = [-10;-10];   % Lower bounds for x.
x_U     = [2;2];       % Upper bounds for x.
fLowBnd = 0;         % Lower bound on function.

% Generate the problem structure using the TOMLAB format
Prob = conAssign('uc1_f', 'uc1_g', 'uc1_H', [], x_L, x_U, Name, x_0, ...
           [], fLowBnd);

Result  = tomRun('ucSolve',Prob);

PrintResult(Result);

echo off

% ---------------------------------------------------------------------
function uc2Demo
% ---------------------------------------------------------------------

format compact
fprintf('========================================================\n');
fprintf('Rosenbrocks banana with explicit f(x) and numerical g(x)\n');
fprintf('========================================================\n');

echo on

Name    = 'RB Banana';
x_0     = [-1.2 1]';   % Starting values for the optimization.
x_L     = [-10;-10];   % Lower bounds for x.
x_U     = [2;2];       % Upper bounds for x.
fLowBnd = 0;           % Lower bound on function

% For illustration, initialize the constraint part and
% call glcAssign with these parameters

% No linear constraints
b_L = []; b_U = []; A = []; 
% No linear constraints
c_L = []; c_U = [];

% For illustration, initialize parameters just used for plotting and prints
x_opt = [1 1];       % Known optimal point (optional).
f_opt = 0;           % Known optimal function value (optional).
x_min = [-1.1 -0.2]; % Plot region parameters.
x_max = [ 1.3  1.3]; % Plot region parameters.

% Give full call to glcAssign, showing the full TOMLAB format
% Only give the function. TOMLAB then estimates any derivatives
% automatically

Prob = glcAssign('uc1_f', x_L, x_U, Name, A, b_L, b_U, ...
                  [], c_L, c_U, x_0, [], [], [], [], ...
                  fLowBnd, x_min, x_max, f_opt, x_opt);

% Illustrate a call to the driver routine instead of direct call to ucSolve
Result = tomRun('ucSolve',Prob);

PrintResult(Result);

% ---------------------------------------------------------------------
function uc3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Sequence of Rosenbrocks banana functions\n');
fprintf('=====================================================\n');

echo on

Name    = 'Generalized RB Banana';
x_0     = [-1.2 1]';   % Starting values for the optimization.
x_L     = [-10;-10];   % Lower bounds for x.
x_U     = [2;2];       % Upper bounds for x.
fLowBnd = 0;           % Lower bound on function (optional).


% Generate the problem structure using the TOMLAB format
% Only use function values.
% ucSolve will use safeguarded BFGS as default
Prob = glcAssign('uc3_f', x_L, x_U, Name, [], [], [], ...
                  [], [], [], x_0, [], [], [], [], ...
                  fLowBnd);                        

% The different steepness parameters to be tried
Steep = [100 500 1000 10000];

for i = 1:4
    % Sending the parameter down to the function routine in Prob.UserParam
    % Any field not used by Tomlab may be used, but it is recommended to
    % use subfields of Prob.user or Prob.userParam

    Prob.userParam.alpha = Steep(i);

    Result = tomRun('ucSolve',Prob);
    PrintResult(Result,1);
    fprintf('Steepness parameter in Rosenbrock Banana function %12.0f\n',...
             Steep(i));
    fprintf('\n\n\n');
end

% ---------------------------------------------------------------------
function uc4Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Running simple 4th order polynomial with Newtons method\n');
fprintf('=====================================================\n');

echo on

Name     = 'Simple 4th order polynomial';
x_0      = [0 0]';        % Starting values for the optimization.
x_0      = [5 2]';        % Starting values for the optimization.
x_L      = [-10;-10];     % Lower bounds for x.
x_U      = [10;10];       % Upper bounds for x.
fLowBnd  = 0;             % Lower bound on function (optional).
pSepFunc = [];            % No separable function assumed
HessPattern = ones(2,2);  % All elements in Hessian are possibly nonzero

% Generate the problem structure using the Tomlab Quick format (short call)
% Use the more general conAssign this time
Prob = conAssign('uc4_f','uc4_g','uc4_H', HessPattern, ...
                  x_L, x_U, Name, x_0, pSepFunc, fLowBnd);

Prob.optParam.IterPrint = 1;
Prob.Solver.Alg         = 1; % Use Newtons Method

Result = tomRun('ucSolve',Prob);
PrintResult(Result,3);
fprintf('\n\n\n');
