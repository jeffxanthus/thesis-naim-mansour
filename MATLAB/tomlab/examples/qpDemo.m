function varargout = qpDemo(varargin)
% QPDEMO M-file for qpDemo.fig
%      QPDEMO, by itself, creates a new QPDEMO or raises the existing
%      singleton*.
%
%      H = QPDEMO returns the handle to a new QPDEMO or the handle to
%      the existing singleton*.
%
%      QPDEMO('Property','Value',...) creates a new QPDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to qpDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      QPDEMO('CALLBACK') and QPDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in QPDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qpDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 13:44:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qpDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @qpDemo_OutputFcn, ...
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


% --- Executes just before qpDemo is made visible.
function qpDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for qpDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes qpDemo wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = qpDemo_OutputFcn(hObject, eventdata, handles)
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
        qp1
    case 3
        qp2
    case 4
        qp3
    case 5
        qp4
    case 6
        qp5
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function qpHelp_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
qpHelp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function qpHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file qpDemo in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'the Tomlab Quick (TQ) Format to formulate and solve QP problems.  ' ...
,' ' ...
,'After defining the Tomlab problem structure called Prob, is easy to set'...
,'different variables in Prob before actually solving the problem.' ...
,' ' ...
,'The call to a Tomlab solver is always very simple, in this case'...
,'the call to qpSolve is just:  Result = qpSolve(Prob);' ...
,' ' ...
,'More general is to use the driver routine tomRun.'...
,'Selecting different solvers is easy: just change the name of the solver:'...
,'    Result = tomRun(''qpSolve'',Prob,1);' ...
,'    Result = tomRun(''qpopt'',Prob,1);' ...
,'    Result = tomRun(''sqopt'',Prob,1);' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function qp1
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Run very simple QP defined in the Tomlab Quick format\n');
fprintf('=====================================================\n');

Name  = 'QP Example';  % File qpExample.m
F     = [ 8   2        % Matrix F in 1/2 * x' * F * x + c' * x
          2   8 ];
c     = [ 3  -4 ]';    % Vector c in 1/2 * x' * F * x + c' * x
A     = [ 1   1        % Constraint matrix
          1  -1 ];
b_L   = [-inf  0  ]';  % Lower bounds on the linear constraints
b_U   = [  5   0  ]';  % Upper bounds on the linear constraints
x_L   = [  0   0  ]';  % Lower bounds on the variables
x_U   = [ inf inf ]';  % Upper bounds on the variables
x_0   = [  0   1  ]';  % Starting point

%   x_min and x_max only needed if doing plots
%   x_min = [-1 -1 ];      % Plot region lower bound parameters
%   x_max = [ 6  6 ];      % Plot region upper bound parameters

% Use the Tomlab Quick (TQ) format
%
% Call the Tomlab qpAssign routine that creates a structure with all
% problem information.

Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, 'qpExample');

Result = qpSolve(Prob);

% When the Projected gradient gPr is very small, the minimum is found
% with good accuracy

PrintResult(Result); 

% ---------------------------------------------------------------------
function qp2
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve sequence of QP problem, defined in the Tomlab Quick format\n');
fprintf('=====================================================\n');

% Define the basic QP problem

Name  = 'QP Example';  % File qpExample.m
F     = [ 8   2        % Matrix F in 1/2 * x' * F * x + c' * x
          2   8 ];
c     = [ 3  -4 ]';    % Vector c in 1/2 * x' * F * x + c' * x
A     = [ 1   1        % Constraint matrix
          1  -1 ];
b_L   = [-inf  0  ]';  % Lower bounds on the linear constraints
b_U   = [  5   0  ]';  % Upper bounds on the linear constraints
x_L   = [  0   0  ]';  % Lower bounds on the variables
x_U   = [ inf inf ]';  % Upper bounds on the variables
x_0   = [  0   1  ]';  % Starting point


Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, 'Random QP');

% Make a test with different randomly disturbed c vectors and F matrices
% The Prob structure need not be redefined 

n=length(c);

fprintf('\nc = %f %f\n\n',c);
mPrint(F,'F');

for i=1:5
    cNew =c + 0.1*(rand(n,1)-0.5); % Generate random disturbances to vector c
    FNew =F + 0.1*(rand(n,n)-0.5); % Generate random disturbances to matrix F

    fprintf('\ncNew = %f %f\n\n',cNew);

    mPrint(FNew,'FNew');

    Prob.QP.c = cNew;
    Prob.QP.F = FNew;

    % Use driver routine, set print level 2 (tomRun calls PrintResult)
    % Save results in a vector of result structures
    Result(i) = tomRun('qpSolve',Prob,2);
    pause(2)
end

% All results are saved in the Result structure, and it is easy to
% do postprocessing on the results

% The objective function value for each of the five runs
f_k = zeros(1,5);
for i = 1:5
    f_k(i) = Result(i).f_k;
end
fprintf('\nf_k (optimum) = ');
fprintf('%9.6f ',f_k);
fprintf('\n\n');

% The CPU time for each of the five runs
CPU = zeros(1,5);
for i = 1:5
    CPU(i) = Result(i).CPUtime;
end
fprintf('\nCPU (seconds) = ');
fprintf('%9.3f ',CPU);
fprintf('\n\n');


% ---------------------------------------------------------------------
function qp3
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Test of randomly generated QP problems, in the TQ format\n');
fprintf('=====================================================\n');

[TomV,os,TV] = tomlabVersion;

fprintf('Run with qpSolve');

if TV(2)
   fprintf(', qld, MINOS, QP-MINOS (special QP interface) and QPOPT\n');
else
   fprintf('and qld\n');
end
if TV(4)
   fprintf('Also try SQOPT and SNOPT\n');
end

fprintf('=====================================================\n');


cases = 5;

NN= 100;  % Random number of variables   around NN
MM= 50;  % Random number of constraints around MM

% Generate random problems
probs = [];
for i=1:cases
   n = round(NN+30*rand);
   m = round(MM+30*rand);
   A = 10*rand(m,n);
   b = 3*10*rand(m,1);
   c = -10*rand(n,1);
   F = zeros(n,n);
   for j = 1:n
       % Add to the diagonal of F
       %F(j,j) = F(j,j) + 10 + 20*rand(1,1);
       F(j,j) = 20 + 20*rand(1,1);
       %if j > 1
       %   F(j-1,j) = -5 + 10*rand(1,1);
       %end
       %if j < n
       %   F(j,j+1) = -5 + 10*rand(1,1);
       %end
   end

   probs(i).n = n;
   probs(i).m = m;
   probs(i).A = A;
   probs(i).b = b;
   probs(i).c = c;
   probs(i).F = F;

end

% Generate the first problem, to make a Prob structure

n   = probs(1).n;
m   = probs(1).m;
A   = probs(1).A;
b_U = probs(1).b;
c   = probs(1).c;
F   = probs(1).F;

Prob = qpAssign(F, c, A, [], b_U, [], [], [], 'Random QP');

% Set same tolerances to be used for all tests
Prob.optParam.MaxIter  = 2000;
Prob.optParam.eps_f    = 1E-8;    
Prob.optParam.eps_Rank = 1E-12; 
Prob.optParam.xTol     = 1E-6;      
Prob.optParam.bTol     = 1E-8;     

% Initial state always the same
rand('state',0);

TomV = tomlabVersion;

for i=1:length(probs)
   n   = probs(i).n;
   m   = probs(i).m;
   disp([n m])

   Prob.Name = ['Random Test ' num2str(i)];
   Prob.P    = i;
   Prob.N    = n;
   Prob.QP.c = probs(i).c; 
   Prob.QP.F = probs(i).F; 
   Prob.A    = probs(i).A;      
   Prob.b_L  = -Inf*ones(m,1);    
   Prob.b_U  = probs(i).b;    
   Prob.x_L  = zeros(n,1);               % Lower bounds set to 0
   Prob.x_U  = inf*ones(n,1);            % Infinite upper bounds 
   Prob.x_0  = [];                       % Infinite upper bounds 
   Prob.mLin = size(Prob.A,1);

   Result    = tomRun('qpSolve',Prob,1);
   Result1   = tomRun('qld',Prob,1);

   if i == 1
      ResSave     = Result;
      ResSave1    = Result1;
   else
      ResSave(i)  = Result ;
      ResSave1(i) = Result1;
   end

   if TV(2)
      Result2   = tomRun('minos',Prob,1);
      % Call the special QP interface to MINOS.
      Result3   = tomRun('qp-minos',Prob,1);
      Result4   = tomRun('qpopt',Prob,1);
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
fprintf('qpSolve  CPU:');
fprintf('%6.3f ',ResSave(1:5).CPUtime);
fprintf('\n');
fprintf('qld      CPU:');
fprintf('%6.3f ',ResSave1(1:5).CPUtime);
fprintf('\n');

if TV(2)
   fprintf('MINOS    CPU:');
   fprintf('%6.3f ',ResSave2(1:5).CPUtime);
   fprintf('\n');
   fprintf('QP-MINOS CPU:');
   fprintf('%6.3f ',ResSave3(1:5).CPUtime);
   fprintf('\n');
   fprintf('QPOPT    CPU:');
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

% MINOS runs sparse QP, and QPOPT dense QP
% It is faster to call QPOPT then QP-MINOS for these dense QP problems.
%
% However, note that the time measurements with cputime and etime in Matlab
% does not seem that reliable

% ---------------------------------------------------------------------
function qp4
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Run QP problem defined in Init File format, alternative 1\n');
fprintf('=====================================================\n');

% Run test problem 4 in qp_prob.m

% Alternative 1 is to first read the Prob structure, then make any
% changes, and finally call to solver

Prob = probInit('qp_prob',4);

% Now any changes can be made before the call

% For illustration, change the initial point.
% The point closest in the feasible set will be chosen.
Prob.x_0 = zeros(10,1);

Result = qpSolve(Prob);

PrintResult(Result);

% ---------------------------------------------------------------------
function qp5
% ---------------------------------------------------------------------
format compact
fprintf('=====================================================\n');
fprintf('Run QP problem defined in Init File format, alternative 2\n');
fprintf('=====================================================\n');

% Run test problem 4 in qp_prob.m

% Alternative 2 is to directly call the multi-solver driver routine tomRun
%
% Any changes in the Prob structure should be made by assigning the
% new values before the call to tomRun

% For illustration, change the initial point.
% The point closest in the feasible set will be chosen.

Prob.x_0 = zeros(10,1);

Result = tomRun('qpSolve','qp_prob',4,Prob,2);

% MODIFICATION LOG
%
% 041208  med  Prob.mLin = size(Prob.A,1); added to qp3