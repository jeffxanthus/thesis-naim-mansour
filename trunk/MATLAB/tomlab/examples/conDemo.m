function varargout = conDemo(varargin)
% CONDEMO M-file for conDemo.fig
%      CONDEMO, by itself, creates a new CONDEMO or raises the existing
%      singleton*.
%
%      H = CONDEMO returns the handle to a new CONDEMO or the handle to
%      the existing singleton*.
%
%      CONDEMO('Property','Value',...) creates a new CONDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to conDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CONDEMO('CALLBACK') and CONDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CONDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help conDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 14:21:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @conDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @conDemo_OutputFcn, ...
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


% --- Executes just before conDemo is made visible.
function conDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for conDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes conDemo wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = conDemo_OutputFcn(hObject, eventdata, handles)
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
        con1Demo
    case 3
        con2Demo
    case 4
        con3Demo
    case 5
        con4Demo
    case 6
        con5Demo
    case 7
        con6Demo
    case 8
        con7Demo
    case 9
        con8Demo
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --------------------------------------------------------------------
function conHelp_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
conHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function conHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '--------------------------------------------- ' ...
,'This is the file conDemo.m in tomlab\examples ' ...
,'--------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'TOMLAB to formulate and solve constrained nonlinear problems (con).' ...
,' ' ...
,'The TOMLAB Quick Format (TQ) is illustrated in the beginning of four of'...
,'the examples.' ...
,' ' ...
,'The code for the examples are all in this file (con1Demo, con2Demo etc.)' ...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The solution of an exponential function is the main example. One linear'...
,'equality constraint and two nonlinear inequality constraints are added.' ...
,'The first example shows a call to conSolve, using the function con1_f,' ...
,'the gradient con1_g, the Hessian con1_H, the constraint routine con1_c ' ...
,'and the constraint gradient matrix con1_dc. The linear constraint is' ...
,'given explicitely. The problem is called EXP' ...
,' ' ...
,'The second example shows how to solve the EXP problem when the user is' ...
,'able to give the function and the constraints, and no derivatives.' ...
,'There are five methods to estimate the derivatives numerically.' ...
,'For three of them, the Spline toolbox is needed. In this example the' ...
,'other two methods are used. ' ...
,'If the Spline toolbox is present, changes should be made to startup.m.' ...
,' ' ...
,'The third example shows how to solve the EXP problem for a number of'...
,'different initial values on x. The initial values are stored in the matrix'...
,'X0, and in each loop step Prob.x_0 is set to one of the columns in X0.' ...
,'In a similar way any of the values in the Prob structure may be changed' ...
,'in a loop step, if e.g. the loop is part of a control loop.' ...
,'The Prob structure only needs to be defined once.' ...
,' ' ...
,'The fourth example shows the same example as the third, but instead'...
,'calling nlpSolve. This routine can utilize second order information'...
,'in the constraints, and therefore con1_d2c is also used.'...
,' ' ...
,'In the fifth example the problem is solved 100 times for different'...
,'starting points, randomly generated. This is one way to secure that'...
,'the correct solution has been found. Another way is to use a global'...
,'optimization solver that handles constrained problems to crunch the'...
,'problem and find the global solution. In TOMLAB the glcSolve solver'...
,'could be used for such a procedure if the number of unknowns are not'...
,'too many. An example calling glcSolve is shown in the sixth example.'...
,'Also see the demonstration for constrained global optimization, glcDemo.m.'...
,' ' ...
,'Note that the call to a TOMLAB solver is always very simple, in this case'...
,'the call to conSolve is just:  Result = conSolve(Prob);' ...
,'We may call the TOMLAB driver routine instead. Then the call is'...
,'   Result = tomRun(''conSolve'',Prob);' ...
,' ' ...
,'To call minos (in TOMLAB /MINOS) use' ...
,'   Result = tomRun(''minos'',Prob);' ...
,' ' ...
,'If trying the global optimization approach, to use the glcSolve routine' ...
,'it is best to use glcAssign, not conAssign,' ...
,'to generate the structure, and then call glcSolve with' ...
,'   Result = tomRun(''glcSolve'',Prob);' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function con1Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('================================================================\n');

echo on

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

% This problem probably have several minima.

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10];    % Two nonlinear inequality constraints
c_U     = [];           % Empty means Inf (default)
x_0     = [-5;5];       % Initial value
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds

x_opt=[-3.162278, 3.162278; -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];

HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable

% Generate the problem structure using the TOMLAB Quick format
Prob = conAssign('con1_f', 'con1_g', 'con1_H', HessPattern, ...
                  x_L, x_U, Name, x_0, ...
                  pSepFunc, fLowBnd, ...
                  A, b_L, b_U, 'con1_c', 'con1_dc', [], ConsPattern, ...
                  c_L, c_U, ... 
                  x_min, x_max, f_opt, x_opt);

% Get default TOMLAB solver for your current license, for "con" problems
Solver = GetSolver('con');

% Call driver routine tomRun, 4th argument > 0 implies call to PrintResult
Result  = tomRun(Solver,Prob,2);

echo off

% ---------------------------------------------------------------------
function con2Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('Numerical estimate of gradient and constraint Jacobian\n');
fprintf('================================================================\n');

echo on

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

% This problem probably have several minima.

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10];    % Two nonlinear inequality constraints
c_U     = [];           % Empty means Inf (default)
x_0     = [-5;5];       % Initial value
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds

x_opt=[-3.162278, 3.162278; -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];

HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable


% Generate the problem structure using the TOMLAB Quick format
Prob = conAssign('con1_f', [], [], HessPattern, ...
                 x_L, x_U, Name, x_0, ...
                 pSepFunc, fLowBnd, ...
                 A, b_L, b_U, 'con1_c', 'con1_dc', [], ConsPattern, ...
                 c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);


% Use standard finite difference method to estimate function derivatives
Prob.NumDiff  = 1;

% Use complex variable method to estimate constraint derivatives
Prob.ConsDiff = 5;

% Default algorithm choice - see help conSolve for more information
Prob.Solver.Alg = 0;

% Call solver through tomRun driver routine
Result = tomRun('conSolve',Prob,2);

% Same as doing
%   Result  = conSolve(Prob);
%   PrintResult(Result);

echo off

% ---------------------------------------------------------------------
function con3Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('Testing start with different initial values\n');
fprintf('================================================================\n');

echo on

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10]     % Two nonlinear inequality constraints
c_U     = [];           % Empty means Inf (default)
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds


X0      = [ -1 -5  1 0 -5 ;
             1  7 -1 0  5];

x_0     = X0(:,1);

% Starting with different initial points reveals that this problem is nasty.
% There are different points that all fulfill the convergence criteria.
%
% From [-5;5], the routines converge to [-3.16 3.16]
% From [-1;1], the routines converge to [-1.22 1.22]

x_opt=[-3.162278, 3.162278; ...
       -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];


HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable


% Generate the problem structure using the TOMLAB Quick format
Prob = conAssign('con1_f', 'con1_g', 'con1_H', HessPattern, ...
                 x_L, x_U, Name, x_0, ...
                 pSepFunc, fLowBnd, ...
                 A, b_L, b_U, 'con1_c', 'con1_dc', [], ConsPattern, ...
                 c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);

Prob.Solver.Alg = 0;

[TomV,os,TV] = tomlabVersion;

for i = 1:size(X0,2)

    Prob.x_0 = X0(:,i);
    
    if TV(4)
       % If TOMLAB /SOL or TOMLAB /SNOPT
       Result   = tomRun('snopt',Prob);
    elseif TV(3)
       % If TOMLAB /NPSOL
       Result   = tomRun('npsol',Prob);
    elseif TV(2)
       % If TOMLAB /MINOS
       Result   = tomRun('minos',Prob);
    else
       % conSolve is part of TOMLAB base module, always present 
       Result   = tomRun('conSolve',Prob);
    end

    PrintResult(Result);
    if i < size(X0,2)
       disp('Press ENTER to continue');
       pause
    end

end

echo off

% ---------------------------------------------------------------------
function con4Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('Call to nlpSolve, using also  2nd order constraint information\n');
fprintf('================================================================\n');

echo on

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

% This problem probably have several minima.

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10]     % Two nonlinear inequality constraints
c_U     = [Inf;Inf];    % Upper bounds are Inf, i.e. no restrictions
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds

X0      = [ -1 -5  1 0 -5 ;
             1  7 -1 0  5];

x_0     = X0(:,1);

% Starting with different initial points reveals that this problem is nasty.
% There are different points that all fulfill the convergence criteria.
%
% From [-5;5], the routines converge to [-3.16 3.16]
% From [-1;1], the routines converge to [-1.22 1.22]

x_opt=[-3.162278, 3.162278; ...
       -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];

HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable

% Generate the problem structure using the TOMLAB Quick format
% Also using con1_d2c, because nlpSolve can utilize this information
Prob = conAssign('con1_f', 'con1_g', 'con1_H', HessPattern, ...
                  x_L, x_U, Name, x_0, ...
                  pSepFunc, fLowBnd, ...
                  A, b_L, b_U, 'con1_c', 'con1_dc', 'con1_d2c', ConsPattern,...
                  c_L, c_U, ... 
                  x_min, x_max, f_opt, x_opt);

for i = 1:size(X0,2)

    Prob.x_0 = X0(:,i);

    Result  = tomRun('nlpSolve',Prob,2);

    if i < size(X0,2)
       disp('Press ENTER to continue');
       pause
    end

end

echo off

% ---------------------------------------------------------------------
function con5Demo
% ---------------------------------------------------------------------

[TomV,os,TV] = tomlabVersion;

format compact
fprintf('================================================================\n');
fprintf('The Schittkowski 66 test problem\n');
if TV(3)
   fprintf('Call to npsol\n');
elseif TV(4)
   fprintf('Call to snopt\n');
else
   fprintf('Call to conSolve\n');
end
fprintf('The problem is solved 50 times, with random starting points\n');
fprintf('in the box [0,10]x[0,10]x[0,10]\n');
fprintf('================================================================\n');


Prob = probInit('con_prob', 10);

for i=1:50
    Prob.x_0 = [10 10 10]'.*rand(3,1);
    fprintf('\nTest %d  ---  x_0 = %f %f %f\n',i,Prob.x_0);
    if TV(3)
       Result = tomRun('npsol', Prob, 1); 
    elseif TV(4)
       Result = tomRun('snopt', Prob, 1); 
    else
       Result = tomRun('conSolve', Prob, 1); 
    end
end

% ---------------------------------------------------------------------
function con6Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('Use global optimization solver to find the global optimum\n');
fprintf('================================================================\n');

echo on

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

%
% This example is nearly identical with con2Demo, it is just the call
% to glcSolve, not conSolve, that differs 
%

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

% This problem probably have several minima.

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10]     % Two nonlinear inequality constraints
c_U     = [];           % Empty means Inf (default)
x_0     = [-5;5];       % Initial value
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds

x_opt=[-3.162278, 3.162278; ...
       -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];

HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable


% Generate the problem structure using the TOMLAB Quick format
Prob = conAssign('con1_f', [], [], HessPattern, ...
                 x_L, x_U, Name, x_0, ...
                 pSepFunc, fLowBnd, ...
                 A, b_L, b_U, 'con1_c', 'con1_dc', [], ConsPattern, ...
                 c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);

% Use very many function evaluations. 
% See glcDemo on how to do warm starts instead, enabling checks on the
% intermediate results.

fprintf('\nA first call to glcSolve with the default number of function ');
fprintf('evaluations (200)\n');

Result  = tomRun('glcSolve',Prob,2);

fprintf('Note that it is better to use glcAssign, not conAssign, ');
fprintf('to make the structure\n');
fprintf('But it still works with conAssign\n\n');

Prob.optParam.MaxFunc = 2000;
Prob.optParam.MaxIter = 5000; % Will never be reached, MaxFunc == 2000 will stop
Prob.optParam.IterPrint = 1;  % One line / iteration

fprintf('Increase the number of function evaluations to %d\n',...
         Prob.optParam.MaxFunc);

fprintf('Still the same accuracy as conSolve/nlpSolve will not be achieved\n');

fprintf('\nThis will take some time :-). Please wait...\n');
fprintf('\nPrint one line per iteration.\n');
disp('Press ENTER to continue');
pause

Result  = glcSolve(Prob);

PrintResult(Result);

echo off

% ---------------------------------------------------------------------
function con7Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('Solve with SOL solvers NPSOL and SNOPT\n');
fprintf('================================================================\n');

echo on

[TomV,os,TV] = tomlabVersion;

if all( TV(3:4)==0 )
   fprintf('No license for TOMLAB /NPSOL, TOMLAB /SNOPT or TOMLAB /SOL\n');
   return
end

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10]     % Two nonlinear inequality constraints
c_U     = [];           % Empty means Inf (default)
x_0     = [-5;5];       % Initial value
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds

x_opt=[-3.162278, 3.162278; -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];

HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable

% Generate the problem structure using the TOMLAB Quick format
Prob = conAssign('con1_f', 'con1_g', 'con1_H', HessPattern, ...
                  x_L, x_U, Name, x_0, ...
                  pSepFunc, fLowBnd, ...
                  A, b_L, b_U, 'con1_c', 'con1_dc', [], ConsPattern, ...
                  c_L, c_U, ... 
                  x_min, x_max, f_opt, x_opt);

% Set Prob.LargeScale = 1 to avoid any bookkeeping in TOMLAB of search steps
Prob.LargeScale = 1;

if TV(3)
   % Increase Print Level for NPSOL to 1 (maximum is 50)
   Prob.SOL.optPar(1) = 1;

   % To save the output from NPSOL on a text file, define

   Prob.SOL.PrintFile = 'npsol.txt';   % The print file
   Prob.SOL.SummFile  = 'npsols.txt';  % The summary file
 
   % Call NPSOL using TOMLAB driver routine

   % Add 2 as input to tomRun, to get a call: PrintResult(Result,2);

   Result  = tomRun('npsol',Prob,2);

   disp('Press ENTER to continue');
   pause
end

if TV(4)
   % To save the output from SNOPT on a text file, define

   Prob.SOL.PrintFile = 'snopt.txt';   % The print file
   Prob.SOL.SummFile  = 'snopts.txt';  % The summary file
   
   % Increase Print Level for SNOPT to maximum, default is 0, no output.
   Prob.SOL.optPar(1) = 111111;
   
   % Call SNOPT using TOMLAB driver routine
   
   Result  = tomRun('snopt',Prob,2);
end

echo off

% ---------------------------------------------------------------------
function con8Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Exponential problem, 1 linear equality, 2 nonlinear inequalities\n');
fprintf('Numerical estimate of gradient and constraint Jacobian\n');
fprintf('Solve with SOL solvers NPSOL and SNOPT\n');
fprintf('================================================================\n');

[TomV,os,TV] = tomlabVersion;

if all( TV(3:4)==0 )
  fprintf('No license for TOMLAB /NPSOL, TOMLAB /SNOPT or TOMLAB /SOL\n');
  return
end

echo on

% This is problem 3 in con_prob.m, the standard con test problems
% 2 nonlinear inequalities + 1 linear equality. Not x_L >=0, leads to problems

Name='Exponential problem 3. 1 linear eq.+ 2 nonlinear ineq.';

A       = [1 1];        % One linear constraint
b_L     = 0;            % Lower bound on linear constraint
b_U     = 0;            % b_L == b_U implies equality
c_L     = [1.5;-10]     % Two nonlinear inequality constraints
c_U     = [];           % Empty means Inf (default)
x_0     = [-5;5];       % Initial value
x_L     = [-10;-10];    % Lower bounds on x
x_U     = [10;10];      % Upper bounds on x
fLowBnd = 0;            % A lower bound on the optimal function value
x_min   = [-2;-2];      % Used for plotting, lower bounds
x_max   = [4;4];        % Used for plotting, upper bounds

x_opt=[-3.162278, 3.162278; -1.224745, 1.224745];
f_opt=[1.156627; 1.8951];

HessPattern = [];  % All elements in Hessian are nonzero. 
ConsPattern = [];  % All elements in the constraint Jacobian are nonzero. 
pSepFunc    = [];  % The function f is not defined as separable


% Generate the problem structure using the TOMLAB Quick format
Prob = conAssign('con1_f', [], [], HessPattern, ...
                 x_L, x_U, Name, x_0, ...
                 pSepFunc, fLowBnd, ...
                 A, b_L, b_U, 'con1_c', 'con1_dc', [], ConsPattern, ...
                 c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);


% If using the internal derivatives in SNOPT and NPSOL, set
% NumDiff and ConsDiff = 6, otherwise 1-5 is possible to use.
Prob.NumDiff  = 6;
Prob.ConsDiff = 6;

% Tell NPSOL and SNOPT that neither function gradients, nor constraint
% gradients are available
%
% Prob.SOL.optPar(39) = 0;  No derivative information available
% Prob.SOL.optPar(39) = 1;  Constraint derivatives are missing
% Prob.SOL.optPar(39) = 2;  Function derivatives are missing
% Prob.SOL.optPar(39) = 3;  Both types of derivatives are available

Prob.SOL.optPar(39) = 0;

% Set Prob.LargeScale = 1 to avoid any bookkeeping in TOMLAB of search steps
Prob.LargeScale = 1;

% To change the feasibility tolerance for NPSOL and SNOPT, one way is to
% directly change the input in the vector Prob.SOL.optPar:
% Higher than 1E-5 is not recommended.

Prob.SOL.optPar(9) = 1E-5; % Major feasibility tolerance nonlinear constraints

% To change the convergence tolerance in x for NPSOL and SNOPT, one way is to
% directly change the input in the vector Prob.SOL.optPar:

Prob.SOL.optPar(10) = 1E-5; % Major feasibility tolerance in x

if TV(3)
   % To save the output from NPSOL on a text file, define

   Prob.SOL.PrintFile = 'npsol.txt';   % The print file
   Prob.SOL.SummFile  = 'npsols.txt';  % The summary file

   % Increase Print Level for NPSOL to 10 (default 0, maximum is 50)
   Prob.SOL.optPar(1) = 10;


   % Call NPSOL using TOMLAB driver routine

   % Add 2 as input to tomRun, to get a call: PrintResult(Result,2);

   Result  = tomRun('npsol',Prob,2);

   disp('Press ENTER to continue');
   pause
end

if TV(3)
   % To save the output from SNOPT on a text file, define

   Prob.SOL.PrintFile = 'snopt.txt';   % The print file
   Prob.SOL.SummFile  = 'snopts.txt';  % The summary file

   % Increase Print Level for SNOPT to maximum, default is 0, no output.
   Prob.SOL.optPar(1) = 111111;

   % Call SNOPT using TOMLAB driver routine

   Result  = tomRun('snopt',Prob,2);
end

echo off
