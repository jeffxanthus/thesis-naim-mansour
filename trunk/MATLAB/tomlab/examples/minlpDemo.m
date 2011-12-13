function varargout = minlpDemo(varargin)
% MINLPDEMO M-file for minlpDemo.fig
%      MINLPDEMO, by itself, creates a new MINLPDEMO or raises the existing
%      singleton*.
%
%      H = MINLPDEMO returns the handle to a new MINLPDEMO or the handle to
%      the existing singleton*.
%
%      MINLPDEMO('Property','Value',...) creates a new MINLPDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to minlpDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MINLPDEMO('CALLBACK') and MINLPDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MINLPDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help minlpDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 14:40:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @minlpDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @minlpDemo_OutputFcn, ...
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


% --- Executes just before minlpDemo is made visible.
function minlpDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for minlpDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes minlpDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = minlpDemo_OutputFcn(hObject, eventdata, handles)
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
        minlp1Demo
    case 3
        minlp2Demo
    case 4
        minlp3Demo
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
function minlpHelp_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minlpHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function minlpHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '----------------------------------------------- ' ...
,'This is the file minlpDemo.m in tomlab\examples ' ...
,'----------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'TOMLAB to formulate and solve constrained mixed-integer nonlinear problems (minlp).' ...
,' ' ...
,'The TOMLAB Quick Format (TQ) is illustrated in the beginning of the three'...
,'examples.' ...
,' ' ...
,'The code for the examples are all in this file (minlp1Demo, minlp2Demo etc.)' ...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The two first examples show setting up and solving a mixed-integer'...
,'problem with both linear and nonlinear constraints. '...
,' '...
,'The second example solves the same problem as the first, but with numerical'...
,'differentiation for both function and constraints gradients and Hessians.'...
,' ' ...
,'The third example shows a pure integer problem.'...
,' ' ...
,'Note that the call to a TOMLAB solver is always very simple. In this case'...
,'we first use GetSolver(''minlp'') to obtain the name of a suitable solver.'...
,' '...
,'We then call the TOMLAB driver routine with: '...
,'   Result = tomRun(Solver,Prob,2);' ...
,' ' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function minlp1Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Linear problem, 2 continuous and 3 integer variables,\n');
fprintf('3 linear and 2 nonlinear constraints.\n');
fprintf('================================================================\n');

echo on

% This is problem 1 in minlp_prob.m, the standard minlp test problems

Name='minlp1Demo - Kocis/Grossman.';

IntVars   = logical([ 0 0 1 1 1 ]); % Integer variables: x(3)-x(5)
VarWeight = [ ];                    % No priorities given 

% There are divisions and square roots involving x(2), so we must
% have a small but positive value for the lower bound on x(2). 
BIG = 1E8;

x_L = [ 0   1/BIG   0 0 0 ]';    % Lower bounds on x
x_U = [ BIG  BIG    1 1 1 ]'; 	 % Upper bounds on x

% Three linear constraints
A = [1 0     1  0 0 ; ...
     0 1.333 0  1 0 ; ...
     0 0    -1 -1 1 ];

b_L = [];               % No lower bounds
b_U = [1.6 ; 3 ; 0];    % Upper bounds

c_L = [1.25;3];         % Two nonlinear constraints
c_U = c_L;              % c_L==c_U implies equality

x_0 = ones(5,1);        % Initial value
  
x_opt = [1.12,1.31,0,1,1]';  % One optimum known
f_opt = 7.6672;              % Value f(x_opt)

x_min = [-1 -1 0 0 0];  % Used for plotting, lower bounds
x_max = [ 1  1 1 1 1];	% Used for plotting, upper bounds

HessPattern = spalloc(5,5,0);  % All elements in Hessian are zero. 
ConsPattern = [ 1 0 1 0 0; ...  % Sparsity pattern of nonlinear
	        0 1 0 1 0 ];    % constraint gradient

fIP = [];        % An upper bound on the IP value wanted. Makes it possible
xIP = [];        % to cut branches. xIP: the x value giving fIP

% Generate the problem structure using the TOMLAB Quick format
Prob = minlpAssign('minlp_f', 'minlp_g', 'minlp_H', HessPattern, ...
                  x_L, x_U, Name, x_0, ...
                  IntVars, VarWeight, fIP, xIP, ...
                  A, b_L, b_U, 'minlp_c', 'minlp_dc', 'minlp_d2c', ...
		  ConsPattern, c_L, c_U, ... 
                  x_min, x_max, f_opt, x_opt);

Prob.P = 1;   % Needed in minlp_xxx files

% Get default TOMLAB solver for your current license, for "minlp" problems
Solver = GetSolver('minlp');

% Call driver routine tomRun, 3rd argument > 0 implies call to PrintResult
Result  = tomRun(Solver,Prob,2);

fprintf('Note that minlpBB does not converge to the global optimum from the ');
fprintf('initial point x_0 = (1,1,1,1,1)\n');
fprintf('\n');
fprintf('The reason is that (1,1,1,1,1) is not feasible with respect ');
fprintf('to the linear constraints, and then\n');
fprintf('minlpBB does a Phase I and finds a feasible point (0.6,1,1,1,1) \n');
fprintf('From this starting point minlpBB does not converge for this');
fprintf('non-convex problem');
fprintf('\n');
fprintf('It claims that it is optimal, but it has only found a local ');
fprintf('optimum');
fprintf('\n');
fprintf('Let us first start from (0.6,1,1,1,1) to see that minlpBB ');
fprintf('only finds the local \noptimum from this feasible point');
fprintf('\n');
fprintf('Press return to continue\n');
pause
% Set new initial point
Prob.x_0 = [0.6, 1, 1, 1, 1]';

xprint(Prob.x_0,'x_0');

Result  = tomRun(Solver,Prob,2);

fprintf('minlpBB did not converge to the global optimum ');
fprintf('from the point (0.6,1,1,1,1)\n');
fprintf('\n');
fprintf('Setting Prob.DUNDEE.optPar(20) = 1 makes minlpBB treat all ');
fprintf('constraints as nonlinear');
fprintf('\n');
fprintf('Then minlpBB will NOT try to make linear constraints feasible\n');
fprintf('\n');
fprintf('This will make minlpBB converge to the correct global optimum !!!\n');
fprintf('\n');
fprintf('Press return to continue\n');
pause
% Set original initial point
Prob.x_0 = [1, 1, 1, 1, 1]';
xprint(Prob.x_0,'x_0');

% Set option Nonlin (optPar(20)) that makes minlpBB skip feasibility tests
% for linear constraints
Prob.DUNDEE.optPar(20) = 1;

Result  = tomRun(Solver,Prob,2);
echo off


% ---------------------------------------------------------------------
function minlp2Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Linear problem, 2 continuous and 3 integer variables,\n');
fprintf('3 linear and 2 nonlinear constraints.\n');
fprintf('Use numerical differentiation everywhere\n');
fprintf('================================================================\n');

echo on

% This is problem 1 in minlp_prob.m, the standard minlp test problems

Name='minlp2Demo - K/G, num. diff.';

% IntVars   = logical([ 0 0 1 1 1 ]); % Integer variables: x(3)-x(5)
IntVars   = [ 3 4 5];               % Integer variables: x(3)-x(5)
VarWeight = [ ];           % No priorities given 


% There are divisions and square roots involving x(2), so we must
% have a small but positive value for the lower bound on x(2). 
BIG = 1E8;

x_L = [ 0   1/BIG   0 0 0 ]';    % Lower bounds on x
x_U = [ BIG  BIG    1 1 1 ]'; 	 % Upper bounds on x

% Three linear constraints
A = [1 0     1  0 0 ; ...
     0 1.333 0  1 0 ; ...
     0 0    -1 -1 1 ];

b_L = [];               % No lower bounds
b_U = [1.6 ; 3 ; 0];    % Upper bounds

c_L = [1.25;3];         % Two nonlinear constraints
c_U = c_L;              % c_L==c_U implies equality

x_0 = ones(5,1);        % Initial value
  
x_opt = [1.12,1.31,0,1,1]';  % One optimum known
f_opt = 7.6672;              % Value f(x_opt)

x_min = [-1 -1 0 0 0];  % Used for plotting, lower bounds
x_max = [ 1  1 1 1 1];	% Used for plotting, upper bounds

HessPattern = spalloc(5,5,0);  % All elements in Hessian are zero. 
ConsPattern = [ 1 0 1 0 0; ...  % Sparsity pattern of nonlinear
	        0 1 0 1 0 ];    % constraint gradient

fIP = [];        % An upper bound on the IP value wanted. Makes it possible
xIP = [];        % to cut branches. xIP: the x value giving fIP

% Generate the problem structure using the TOMLAB Quick format
% Giving ONLY minlp_f and minlp_c, NO derivate routines, c.f. minlp1Demo
Prob = minlpAssign('minlp_f', 'minlp_g', '', HessPattern, ...
                  x_L, x_U, Name, x_0, ...
                  IntVars, VarWeight, fIP, xIP, ...
                  A, b_L, b_U, 'minlp_c', 'minlp_dc', '', ...
		  ConsPattern, c_L, c_U, ... 
                  x_min, x_max, f_opt, x_opt);

Prob.P = 1;   % Needed in minlp_xxx files

% Get default TOMLAB solver for your current license, for "minlp" problems
Solver = GetSolver('minlp');

% Call driver routine tomRun, 4th argument > 0 implies call to PrintResult
Result  = tomRun(Solver,Prob,2);

echo off



% ---------------------------------------------------------------------
function minlp3Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Pure integer problem, 2 variables.\n');
fprintf('3 linear and 2 nonlinear constraints.\n');
fprintf('================================================================\n');

echo on


% This is problem 4 in minlp_prob.m, the standard minlp test
% problems.

Name = 'minlp3Demo - Pörn 1997';

IntVars   = [1:2];  % Make x(1)-x(2) integers 
VarWeight = [];     % Give no priorities

x_L = [1 1]';     % Variable bounds, lower
x_U = [5 5]';     % ...and upper

x_0 = x_L;        % Initial value
 

% Three linear constraints
A = [ -1 -2  ; ...
      -3  1  ; ...
       4 -3  ];


b_L = -inf*ones(3,1);  % Linear constraints bounds, no lower bounds
b_U = [ 5 1 11 ]';     % Upper bounds

c_U = -24;             % Nonlinear constraint is -inf <= c(x) <= -24
c_L = -inf;
  
ConsPattern = [];             % All elements of grad c(x) are nonzero
HessPattern = spalloc(2,2,0); % Linear obj.function => zero Hessian

f_opt = 31;        % One optimal point known
x_opt = [3 1]'; 

x_min = [];       % Limits for plotting
x_max = [];

fIP   = [];       % An upper bound on the IP value
xIP   = [];       % wanted. Makes it possible to cut branches

% Generate the problem structure using the TOMLAB Quick format
Prob = minlpAssign('minlp_f', 'minlp_g', 'minlp_H', HessPattern, ...
		   x_L, x_U, Name, x_0, ...
		   IntVars, VarWeight, fIP, xIP, ...
		   A, b_L, b_U, 'minlp_c', 'minlp_dc', 'minlp_d2c', ...
		   ConsPattern, c_L, c_U, ... 
		   x_min, x_max, f_opt, x_opt);

Prob.P = 4;     % Needed in minlp_xxx files

% Get default TOMLAB solver for your current license, for "minlp" problems
Solver = GetSolver('minlp');

% Call driver routine tomRun, 3rd argument > 0 implies call to PrintResult
Result = tomRun(Solver,Prob,2);

echo off
