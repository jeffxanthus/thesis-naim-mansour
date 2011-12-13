function varargout = sdpDemo(varargin)
% SDPDEMO M-file for sdpDemo.fig
%      SDPDEMO, by itself, creates a new SDPDEMO or raises the existing
%      singleton*.
%
%      H = SDPDEMO returns the handle to a new SDPDEMO or the handle to
%      the existing singleton*.
%
%      SDPDEMO('Property','Value',...) creates a new SDPDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to sdpDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SDPDEMO('CALLBACK') and SDPDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SDPDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sdpDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 14:44:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sdpDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @sdpDemo_OutputFcn, ...
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


% --- Executes just before sdpDemo is made visible.
function sdpDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for sdpDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sdpDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sdpDemo_OutputFcn(hObject, eventdata, handles)
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
        sdp1Demo
    case 3
        sdp2Demo
    case 4
        sdp3Demo
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
function sdpHelp_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sdpHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function sdpHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file sdpDemo.m in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'Tomlab to formulate and solve semidefinite programming (sdp).' ...
,' ' ...
,'The code for the examples are all in this file (sdp1Demo, sdp2Demo etc.)' ...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The first example shows how to solve a problem defined in the SDPLIB '...
,'standard format. It is converted with the routine sdp2pen to a Matlab.' ...
,'structure. The PENSDP is called directly using an interface routine pen' ...
,' ' ...
,'The second example shows how to solve an LP problem formulated in the' ...
,'Tomlab format with PENSDP.' ...
,' ' ...
,'The third example shows how to solve an LP problem formulated in the'...
,'PENSDP format with Tomlab. '...
,'The sdpAssign routine is used to create a Tomlab problem structure Prob'...
,'Then the Tomlab driver routine, tomRun, is called to solve the problem'...
,'and return a standard Tomlab result structure.'...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function sdp1Demo
% ---------------------------------------------------------------------

format compact
fprintf('================================================================\n');
fprintf('Call SDPLIB test problem\n');
fprintf('================================================================\n');

echo on

disp('Read and convert the problem arch0.dat-s from SDPLIB')

sdpFile = 'arch0.dat-s';
p = sdpa2pen(sdpFile);

if isempty(p)
   fprintf('Error! Empty conversion. SDPLIB file %s not found',sdpFile);
else
   disp('The structure that PENSDP is using is shown below.')
   disp('The fields correspond to the pensdp.ps manual in tomlab\docs.')

   p

   fprintf('Default only p.ioptions(1)=0 is set\n')
   fprintf('This implies that the default ioptions are set to:\n')

   p.ioptions = [1 50 100 1 0 0 0 0];

   fprintf('%d ',p.ioptions)
   fprintf('\n');
   fprintf('Setting ioptions(4) = 0; ')
   fprintf('PENSDP would only write a summary output\n')


   p.ioptions(4) = 1; % Output level, 0 or 1

   fprintf('Setting ioptions(5) = 1; ')
   fprintf('PENSDP would not check density of Hessian, ')
   fprintf('necessary for larger problems.\n')

   p.ioptions(5) = 0; % 0 = Check density for Hessian, 1 = do not check

   % If true, PENSDP will do line searches
   p.ioptions(6) = 0; % Line search
   
   fprintf('\nDefault p.foptions is set as empty\n')
   fprintf('This implies that the default foptions are set to:\n')
   p.foptions = [1 0.7 0.1 1E-7 1E-6 1E-14 1E-2 1];
   fprintf('%f ',p.foptions)
   fprintf('\n');


   fprintf('\n\nCall PENSDP:\n')

   [x, fx, uoutput] = pen(p);

   fprintf('\n\n')
   fprintf('Best f(x) from PENSDP = %40.20f\n\n',fx);
   fprintf('PENSDP also returns the optimal x vector and ');
   fprintf('the Lagrange multipliers\n');

end


echo off

% ---------------------------------------------------------------------
function sdp2Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Run PENSDP on the simple LP problem defined in Init File format\n');
fprintf('in file tomlab\\testprob\\lp_prob.m\n');
fprintf('=====================================================\n');


% First read the Prob structure, then make any
% changes, and finally call to solver PENSDP

for P = 1:12
    Prob = probInit('lp_prob',P);

    % Now any changes can be made before the call

    %Prob.PENSDP.ioptions(1) = 100;

    Result = tomRun('pensdp',Prob,2);
end

echo off

% ---------------------------------------------------------------------
function sdp3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Run PENSDP on the simple LP problem defined in PENSDP \n');
fprintf('structure format. The Tomlab driver routine tomRun is used.\n');
fprintf('The sdpAssign routine is used to create a Tomlab problem ');
fprintf('structure Prob\n');
fprintf('Then the Tomlab driver routine, tomRun, is used to solve\n');
fprintf('the problem and return a standard Tomlab result structure.\n');
fprintf('=====================================================\n');

%
% Solve min x_1 + x_2
%
% s/t   -x_1         <=  0
%            -   x_2 <=  0
%       -x_1 - 2 x_2 <= -1
%
% initial x = [ 0 0.1];

% Number of variables
p.vars    = 2;
% Number of linear constraints
p.constr  = 3;
% Number of linear matrix inequalities
p.mconstr = 0;
p.msizes  = [];
p.x0      = [0.,0.1];
p.fobj    = [1,1];
p.ci      = [0,0,-1];
p.bi_dim  = [1,1,2];
p.bi_idx  = [0,1,0,1];
p.bi_val  = [-1,-1,-1,-2];
p.ai_dim  = [0];
p.ai_idx  = [0];
p.ai_nzs  = [0];
p.ai_val  = [0];
p.ai_col  = [0];
p.ai_row  = [0];

p.ioptions = [1 50 100 1 0 0 0 0]; 
p.foptions = [1 0.7 0.1 1.0E-7 1.0E-6 1.0E-14 1.0E-2 1];

Prob = sdpAssign(p);

% If setting the initial value in the Prob.x_0 field, it will be used
% Prob.x_0      = zeros(2,1);

Result = tomRun('pensdp',Prob,2);

echo off
