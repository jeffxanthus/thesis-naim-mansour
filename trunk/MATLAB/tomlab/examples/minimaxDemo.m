function varargout = minimaxDemo(varargin)
% MINIMAXDEMO M-file for minimaxDemo.fig
%      MINIMAXDEMO, by itself, creates a new MINIMAXDEMO or raises the existing
%      singleton*.
%
%      H = MINIMAXDEMO returns the handle to a new MINIMAXDEMO or the handle to
%      the existing singleton*.
%
%      MINIMAXDEMO('Property','Value',...) creates a new MINIMAXDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to minimaxDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MINIMAXDEMO('CALLBACK') and MINIMAXDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MINIMAXDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help minimaxDemo

% Last Modified by GUIDE v2.5 27-Oct-2003 09:23:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @minimaxDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @minimaxDemo_OutputFcn, ...
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


% --- Executes just before minimaxDemo is made visible.
function minimaxDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for minimaxDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes minimaxDemo wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = minimaxDemo_OutputFcn(hObject, eventdata, handles)
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
val = get(hObject,'Value');
switch val
    case 2
        mima1Demo
    case 3
        mima2Demo
    case 4
        mima3Demo
end
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --------------------------------------------------------------------
function mimaHelp_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mimaHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function mimaHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file minimaxDemo.m in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown examples that illustrates the use of'...
,'Tomlab to formulate and solve minimax problems.' ...
,' ' ...
,'The code for the examples are all in this file (mima1Demo, etc.)' ...
,'The names are listed in the beginning of this file (minimaxDemo.m)' ...
,' ' ...
,'The first example shows how to solve a minimax problem '...
,'where the maximum is taken over a set of absolute values of residuals.'...
,'If Tomlab /SOL, then the SOL solver SNOPT is used to solve the problem.'...
,'If Tomlab /MINOS, then the SOL solver MINOS is used to solve the problem.'...
,'Otherwise conSolve is used to solve the problem.'...
,' ' ...
,'The second example is similar to the first.'...
,'A linear constraint is added to the problem formulation.'...
,' ' ...
,'The third example is the same as the second.'...
,'But numerical differences is used instead of an analytic Jacobian.'...
,' ' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function mima1Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Minimax example, 2 variables, 3 residual functions\n');
fprintf('================================================================\n');

echo on

Name = 'Madsen-Tinglett 1';
x_0  = [1;1];        % Initial value
x_L  = [-10;-10];    % Lower bounds on x
x_U  = [10;10];      % Upper bounds on x

% Solve the problem min max |r_i(x)|, where i = 1,2,3
% Solve the problem by eliminating abs, doubling the residuals, reverse sign
% i.e. min max [r_1; r_2; r_3; -r_1; -r_2; -r_3];
y = [-1.5; -2.25; -2.625]; % The data values
y = [y ; -y];           % Eliminate abs, double the residuals, reverse sign
t = [];                 % No time vector used


% Generate the problem structure using the Tomlab Quick format
%
% The part in the residuals dependent on x are defined in mima_r.m 
% The Jacobian is defined in mima_J.m

Prob = clsAssign('mima_r', 'mima_J', [], x_L, x_U, Name, x_0, y, t);

% Set JacPattern, infSolve will compute correct ConsPattern
Prob.JacPattern  = ones(6,2);

% Set the optimal values into the structure (for nice result presenting)

Prob.x_opt=[3, 0.5];
Prob.f_opt=0;

% Use the standard method in infSolve
Prob.InfType = 1;

% Get the default solver for the current license
Solver = GetSolver('con',1,0);

Prob.SolverInf = Solver;

% One may set the solver directly:
%Prob.SolverInf = 'snopt';
%Prob.SolverInf = 'npsol';
%Prob.SolverInf = 'minos';
%Prob.SolverInf = 'conSolve';

% These statements generate ASCII log files when running SNOPT
if strcmpi(Solver,'snopt') 
   Prob.SOL.PrintFile = 'snopt.txt';
   Prob.SOL.SummFile  = 'snopts.txt';
   Prob.SOL.optPar(1) = 10;   % Print level in SNOPT
   %Prob.SOL.optPar(13) = 3;   % Verify level in SNOPT
end

% Set print level 2 to get output from PrintResult at the end
PriLev = 2;
Prob.PriLevOpt = 0;

Result  = infSolve(Prob, PriLev);

echo off

% ---------------------------------------------------------------------
function mima2Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Minimax example, 2 variables, 3 residuals, 1 linear inequality\n');
fprintf('================================================================\n');

echo on

Name = 'Madsen-Tinglett 2';
x_0  = [1;1];        % Initial value
x_L  = [-10;-10];    % Lower bounds on x
x_U  = [10;10];      % Upper bounds on x

% Solve the problem min max |r_i(x)|, where i = 1,2,3
% Solve the problem by eliminating abs, doubling the residuals, reverse sign
% i.e. min max [r_1; r_2; r_3; -r_1; -r_2; -r_3];
y = [-1.5; -2.25; -2.625]; % The data values
y = [y ; -y];           % Eliminate abs, double the residuals, reverse sign
t = [];                 % No time vector used

%
% Add the linear constraint -x(1) + x(2) + 2 >= 0
% Write the constraint as x(1) - x(2) <= 2

% The A matrix could be specified dense or sparse
% A   = sparse([1 -1]);

A   = [1 -1];
b_L = -inf;
b_U = 2;

% Generate the problem structure using the Tomlab Quick format
%
% The part in the residuals dependent on x are defined in mima_r.m 
% The Jacobian is defined in mima_J.m

Prob = clsAssign('mima_r', 'mima_J', [], x_L, x_U, Name, x_0, y, t, ...
                 [],[],[],[],A,b_L,b_U);

% Set the optimal values into the structure (for nice result presenting)

Prob.x_opt=[2.3660254038, 0.3660254038];
% Optimal residuals:
r_opt = [0; 0.2009618943;0.375];

% Compute optimal function value:
Prob.f_opt=max(abs(r_opt));

% Use the standard method in infSolve
Prob.InfType = 1;

% Get the default solver
Solver = GetSolver('con',1,0);

Prob.SolverInf = Solver;

% One may set the solver directly:
%Prob.SolverInf = 'snopt';
%Prob.SolverInf = 'npsol';
%Prob.SolverInf = 'minos';
%Prob.SolverInf = 'conSolve';

% These statements generate ASCII log files when running SNOPT
if strcmpi(Solver,'snopt') 
   Prob.SOL.PrintFile = 'snopt.txt';
   Prob.SOL.SummFile  = 'snopts.txt';
   Prob.SOL.optPar(1) = 10;   % Print level in SNOPT
   %Prob.SOL.optPar(13) = 3;   % Verify level in SNOPT
end

% Set print level 2 to get output from PrintResult at the end
PriLev = 2;
Prob.PriLevOpt = 0;

Result  = infSolve(Prob, PriLev);

echo off

% ---------------------------------------------------------------------
function mima3Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Minimax example, 2 variables, 3 residuals, 1 linear inequality\n');
fprintf('Use numerical differences\n');
fprintf('================================================================\n');

echo on

Name = 'Madsen-Tinglett 2 - NumDiff';
x_0  = [1;1];        % Initial value
x_L  = [-10;-10];    % Lower bounds on x
x_U  = [10;10];      % Upper bounds on x

% Solve the problem min max |r_i(x)|, where i = 1,2,3
% Solve the problem by eliminating abs, doubling the residuals, reverse sign
% i.e. min max [r_1; r_2; r_3; -r_1; -r_2; -r_3];
y = [-1.5; -2.25; -2.625]; % The data values
y = [y ; -y];           % Eliminate abs, double the residuals, reverse sign
t = [];                 % No time vector used

%
% Add the linear constraint -x(1) + x(2) + 2 >= 0
% Write the constraint as x(1) - x(2) <= 2

% The A matrix could be specified dense or sparse
% A   = sparse([1 -1]);

A   = [1 -1];
b_L = -inf;
b_U = 2;

% Generate the problem structure using the Tomlab Quick format
%
% The part in the residuals dependent on x are defined in mima_r.m 
% The Jacobian is defined in mima_J.m

Prob = clsAssign('mima_r', [], [], x_L, x_U, Name, x_0, y, t, ...
                 [],[],[],[],A,b_L,b_U);

% Set the optimal values into the structure (for nice result presenting)

Prob.x_opt=[2.3660254038, 0.3660254038];
% Optimal residuals:
r_opt = [0; 0.2009618943;0.375];

% Compute optimal function value:
Prob.f_opt=max(abs(r_opt));

% Use the standard method in infSolve
Prob.InfType = 1;

% Get the default solver
Solver = GetSolver('con',1,0);

Prob.SolverInf = Solver;

% One may set the solver directly:
%Prob.SolverInf = 'snopt';
%Prob.SolverInf = 'npsol';
%Prob.SolverInf = 'minos';
%Prob.SolverInf = 'conSolve';

% Here NumDiff = 1,2,3,4,5,6 is possible
% 1   = Finite differences 
% 2-4 = spline methods, if Spline TB is installed, and Spline TB enabled in
%       startup.m io tomlab directory
% 5   = Complex variable nethod
% 6   = If using a SOL solver, the SOL solver estimates derivatives internally

Prob.NumDiff = 1;

% These statements generate ASCII log files when running SNOPT
if strcmpi(Solver,'snopt') 
   Prob.NumDiff         = 6;  % use internal snopt derivative estimation
   Prob.SOL.PrintFile   = 'snopt.txt';
   Prob.SOL.SummFile    = 'snopts.txt';
   Prob.SOL.optPar(1)   = 10;   % Print level in SNOPT
   %Prob.SOL.optPar(13) = 3;   % Verify level in SNOPT
end

% Set print level 2 to get output from PrintResult at the end
PriLev         = 2;
Prob.PriLevOpt = 0;

% infSolve will reformulate the problem, putting the residuals and
% the nonlinear constraints all in the constraints.
% Therefore the gradient is known, and only the constraints will be
% estimated with numerical differences. infSolve will set
% Prob.NumDiff = 0;  % Because infSolve knows exactly the simple gradient
% Prob.ConsDiff = max(Prob.ConsDiff, Prob.NumDiff);  
%
% If a SOL solver estimates either Jacobian (NumDiff==6 or ConsDiff==6), 
% it must estimate both Jacobians
%
% If both ConsDiff and NumDiff <= 5, then Tomlab will use the original
% given NumDiff for the Jacobian of the residuals, and the original ConsDiff
% for the Jacobian of the constraints.


Result  = infSolve(Prob, PriLev);

echo off
