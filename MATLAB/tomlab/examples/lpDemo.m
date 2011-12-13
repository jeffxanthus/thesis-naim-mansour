function varargout = lpDemo(varargin)
% LPDEMO M-file for lpDemo.fig
%      LPDEMO, by itself, creates a new LPDEMO or raises the existing
%      singleton*.
%
%      H = LPDEMO returns the handle to a new LPDEMO or the handle to
%      the existing singleton*.
%
%      LPDEMO('Property','Value',...) creates a new LPDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to lpDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LPDEMO('CALLBACK') and LPDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LPDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lpDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 13:27:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lpDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @lpDemo_OutputFcn, ...
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


% --- Executes just before lpDemo is made visible.
function lpDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for lpDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lpDemo wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lpDemo_OutputFcn(hObject, eventdata, handles)
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
        lp1
    case 3
        lp2
    case 4
        lp3
    case 5
        lp4
    case 6
        lp5
    case 7
        lp6
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function lpHelp_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lpHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function lpHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file lpDemo in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'TOMLAB to formulate and solve LP problems.  ' ...
,' ' ...
,'It is simple to formulate the problem with the TOMLAB Quick (TQ) Format.'...
,' ' ...
,'After defining the TOMLAB problem structure called Prob, is easy to set'...
,'different variables in Prob before actually solving the problem.' ...
,' ' ...
,'The call to a TOMLAB solver is always very simple, in this case'...
,'the call to lpSimplex is just:  Result = lpSimplex(Prob);' ...
,' ' ...
,'In the Tomlab Base Module the QLD solver is another alternative' ...
,' ' ...
,'The MINOS solver in Tomlab /MINOS is recommended for large problems and ' ...
,'LPOPT or MINOS for small or medium size problems.' ...
,'There is a special LP-interface to MINOS called LP-MINOS.' ...
,' ' ...
,'A more general call is to use the driver routine tomRun.'...
,'Selecting different solvers is easy: just change the name of the solver:'...
,'    Result = tomRun(''lpSimplex'',Prob,1);' ...
,'    Result = tomRun(''qld'',Prob,1);' ...
,'    Result = tomRun(''lp-minos'',Prob,1);' ...
,'    Result = tomRun(''lpopt'',Prob,1);' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function lp1
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Run very simple LP defined in the TOMLAB Quick format\n');
fprintf('=====================================================\n');

Name = 'Simple lp example'

A=[8 5 1 0
  -1 2 0 1 ];    % Matrix A including slack variables

b=[31 6]';       % right-hand side of linear constraints

c=[-30;-18;0;0]; % cost vector, coefficients in linear objective function 

[m,n] = size(A);

x_L = zeros(n,1);  % Lower bounds must be set as 0, otherwise assumed -Inf
x_U = [];          % Upper bounds, default set as infinite

x_0 = [];          % No initial guess of x given

b_U = b;           % Upper bounds set as b
b_L = b;           % Lower bounds equal to upper bounds, now when slacks
                   % have been added

x_min = [];        % Bounds for plotting, normally not used for LP
x_max = [];        % Bounds for plotting, normally not used for LP

f_opt = [];        % The function value at the optimal point
                   % Used in the printout. Assumed not known
x_opt = [];        % The optimal solution is assumed not to be known

f_Low = [];        % f_Low <= f_optimal must hold
                   % Set as -realmax, lowest possible number

setupFile = [];    % Define the Prob structure, but not any permanent Init File,
                   % i.e. do not create a file according to the
                   % TOMLAB Init File (IF) format

nProblem  = [];    % Number of problems in the IF format is 0

% Use the TOMLAB Quick (TQ) format
%
% Call the TOMLAB lpAssign routine that creates a structure with all
% problem information.

Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 f_Low, x_min, x_max, f_opt, x_opt);

% Print out the problem
disp('The matrix A')
A
disp('The right hand side b')
b
disp('The objective function coefficients c')
c
disp('Set lower bounds for x as zero and upper bounds infinite')

% Call the TOMLAB LP solver directly
Result = lpSimplex(Prob);

PrintResult(Result,2);

% ---------------------------------------------------------------------
function lp2
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Run simple LP, no slacks, defined in the TOMLAB Quick format\n');
fprintf('Call both lpSimplex and qpSolve\n');
fprintf('=====================================================\n');

% Define the LP problem, without slack variables

Name  = 'lptest';
c     = [-7 -5]';      % Coefficients in linear objective function 
A     = [ 1  2
             4  1 ];   % Matrix defining linear constraints
b_U   = [ 6 12 ]';     % Upper bounds on the linear inequalities
x_L   = [ 0  0 ]';     % Lower bounds on x

% b_L, x_U and x_0 have default values and need not be defined.
% It is possible to call lpAssign with empty [] arguments instead 
% b_L   = [-inf -inf]';
% x_U   = [];
% x_0   = [];


Prob = lpAssign(c, A, [], b_U, x_L, [], [], 'lpExample');

Result = lpSimplex(Prob);

PrintResult(Result,2); 

% It is easy to call another solver to solve the problem, like a QP solver
Result = qpSolve(Prob);

PrintResult(Result,2); 


% ---------------------------------------------------------------------
function lp3
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Solve sequence of LP problem, defined in the TOMLAB Quick format\n');
fprintf('=====================================================\n');

% Define the LP problem, without slack variables

Name  = 'lptest';
c     = [-7 -5]';      % Coefficients in linear objective function 
A     = [ 1  2
             4  1 ];   % Matrix defining linear constraints
b_U   = [ 6 12 ]';     % Upper bounds on the linear inequalities
x_L   = [ 0  0 ]';     % Lower bounds on x

Prob = lpAssign(c, A, [], b_U, x_L, [], [], 'lpExample');

% Make a test with different randomly disturbed c vectors
% The Prob structure need not be redefined 

n=length(c);

for i=1:5
    cNew =c + 0.1*(rand(n,1)-0.5); % Add random disturbances to vector c

    fprintf('\nc = %f %f\n\n',cNew);
    Prob.QP.c = cNew;

    % Use driver routine, set print level 1 (tomRun calls PrintResult)
    % Save results in a vector of result structures
    Result(i) = tomRun('lpSimplex',Prob,1);
    pause(2)
end

% The objective function value for each of the five runs
for i=1:5
    f_k(i) = Result(i).f_k;
end
xprint(f_k,'f_k(1:5)');


% ---------------------------------------------------------------------
function lp4
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Test of randomly generated LP problems, in the TQ format\n');

[TomV,os,TV] = tomlabVersion;

fprintf('Run with lpSimplex');

if TV(2)
   fprintf(', MINOS, LP-MINOS (special LP interface) and LPOPT\n');
else
   fprintf('\n');
end
if TV(4)
   fprintf('Also try SQOPT and SNOPT\n');
end
fprintf('=====================================================\n');


cases = 5;

NN=480;  % Random number of variables   around NN
MM=220;  % Random number of constraints around MM

% Generate random problems
probs = [];
for i=1:cases
   n = round(NN+50*rand);
   m = round(MM+50*rand);
   A = 10*rand(m,n);
   b = n*10*rand(m,1);
   c = -10*rand(n,1);

   probs(i).n = n;
   probs(i).m = m;
   probs(i).A = A;
   probs(i).b = b;
   probs(i).c = c;

end

% Generate the first problem, to make a Prob structure

n   = probs(1).n;
m   = probs(1).m;
A   = probs(1).A;
b_U = probs(1).b;
c   = probs(1).c;

Prob = lpAssign(c, A, [], b_U, [], [], [], 'Random LP');

% Set same tolerances to be used for all tests
Prob.optParam.MaxIter  =  1000;
Prob.optParam.eps_f    = 1E-6;    
Prob.optParam.eps_Rank = 1E-11; 
Prob.optParam.xTol     = 1E-6;      
Prob.optParam.bTol     = 1E-7;     
Prob = ProbCheck(Prob,'lpSimplex');

% Initial state always the same
rand('state',0);


for i=1:length(probs)
   n   = probs(i).n;
   m   = probs(i).m;
   disp([n m])

   % Change problem
   Prob.Name = ['Random Test ' num2str(i)];
   Prob.P    = i;
   Prob.N    = n;
   Prob.QP.c = probs(i).c; 
   Prob.A    = probs(i).A;      
   Prob.mLin = m;              % Number of linear constraints,size(A,1)
   Prob.b_U  = probs(i).b;    
   Prob.b_L  = -Inf*ones(m,1); % Only upper bounds on linear constraints
   Prob.x_L  = zeros(n,1);     % Lower bounds on variables set to 0
   Prob.x_U  = inf*ones(n,1);  % Infinite upper bounds on variables 

   Result    = tomRun('lpSimplex',Prob,1);

   if i == 1
      ResSave    = Result;
   else
      ResSave(i) = Result;
   end

   if TV(2)
      % Result2   = tomRun('minos',Prob,1);
      Result2   = minosTL(Prob);
      PrintResult(Result2,1);
      % Call the special LP interface to MINOS. Should go faster
      Result3   = tomRun('lp-minos',Prob,1);
      Result4   = tomRun('lpopt',Prob,1);
      if i == 1
         ResSave2    = Result2;
         ResSave3    = Result3;
         ResSave4    = Result4;
      else
         ResSave2(i) = Result2;
         ResSave3(i) = Result3;
         ResSave4(i) = Result4;
      end
   end
   if TV(4)
      Result5   = tomRun('sqopt',Prob,1);
      Result6   = tomRun('snopt',Prob,1);
      if i == 1
         ResSave5    = Result5;
         ResSave6    = Result6;
      else
         ResSave5(i) = Result5;
         ResSave6(i) = Result6;
      end
   end
end

fprintf('\n');
fprintf('lpSimplex  CPU:');
fprintf('%6.3f ',ResSave(1:5).CPUtime);
fprintf('\n');

if TV(2)
   fprintf('MINOS    CPU:');
   fprintf('%6.3f ',ResSave2(1:5).CPUtime);
   fprintf('\n');
   fprintf('MINOS-LP CPU:');
   fprintf('%6.3f ',ResSave3(1:5).CPUtime);
   fprintf('\n');
   fprintf('LPOPT    CPU:');
   fprintf('%6.3f ',ResSave4(1:5).CPUtime);
   fprintf('\n');
end
if TV(4)
   fprintf('SQOPT    CPU:');
   fprintf('%6.3f ',ResSave5(1:5).CPUtime);
   fprintf('\n');
   fprintf('SNOPT    CPU:');
   fprintf('%6.3f ',ResSave6(1:5).CPUtime);
   fprintf('\n');
end

% MINOS runs sparse LP, and LPOPT dense LP
% It is faster to call LP-MINOS then MINOS, but still LPOPT runs much
% faster for these dense LP problems.

% ---------------------------------------------------------------------
function lp5
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Run LP problem defined in Init File format, alternative 1\n');
fprintf('=====================================================\n');

% Run test problem 4 defined in lp_prob.m

% Alternative 1 is to first read the Prob structure, then make any
% changes, and finally call to solver

Prob = probInit('lp_prob',4);

% Now any changes can be made before the call

% For illustration, set the maximal number of iterations to only 2
Prob.optParam.MaxIter = 2;

Result = lpSimplex(Prob);

PrintResult(Result);

% ---------------------------------------------------------------------
function lp6
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Run LP problem defined in Init File format, alternative 2\n');
fprintf('=====================================================\n');

% Run test problem 4 defined in lp_prob.m

% Alternative 2 is to directly call the multi-solver driver routine tomRun
%
% Any changes in the Prob structure should be made by assigning the
% new values before the call to tomRun

% Set the maximal number of iterations to 2 
% lpSimplex will print a conclusion that cycling probably occurred

Prob.optParam.MaxIter = 2;

Result = tomRun('lpSimplex','lp_prob',4,Prob,2);

% MODIFICATION LOG
%
% 041208  med  Prob.N    = n; added to lp4