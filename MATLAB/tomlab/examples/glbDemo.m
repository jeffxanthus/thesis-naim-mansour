function varargout = glbDemo(varargin)
% GLBDEMO M-file for glbDemo.fig
%      GLBDEMO, by itself, creates a new GLBDEMO or raises the existing
%      singleton*.
%
%      H = GLBDEMO returns the handle to a new GLBDEMO or the handle to
%      the existing singleton*.
%
%      GLBDEMO('Property','Value',...) creates a new GLBDEMO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to glbDemo_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GLBDEMO('CALLBACK') and GLBDEMO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GLBDEMO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help glbDemo

% Last Modified by GUIDE v2.5 22-Oct-2003 14:53:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @glbDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @glbDemo_OutputFcn, ...
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


% --- Executes just before glbDemo is made visible.
function glbDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for glbDemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes glbDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = glbDemo_OutputFcn(hObject, eventdata, handles)
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
        glb1Demo
    case 3
        glb2Demo
    case 4
        glb3Demo
    case 5
        glb4Demo
end

% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 

% --------------------------------------------------------------------
function glbHelp_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glbHelp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
function glbHelp
% ---------------------------------------------------------------------

Z = str2mat( ...
 '------------------------------------------- ' ...
,'This is the file glbDemo.m in tomlab\examples ' ...
,'------------------------------------------- ' ...
,' ' ...
,'Here are shown a number of different examples that illustrates the use of'...
,'the TOMLAB Quick (TQ) and Init File (IF) format to solve glb problems.  ' ...
,' ' ...
,'The code for the examples are all in this file (glb1Demo, glb2Demo etc.)'...
,'The names are found in the beginning of the file.)' ...
,' ' ...
,'The solution of the Shekel 5 test problem is used as the main example' ...
,' ' ...
,'Note that the call to a TOMLAB solver is always very simple:'...
,'   Result = glbSolve(Prob);' ...
,' ' ...
,'We may call the TOMLAB driver routine instead. Then the call is'...
,'   Result = tomRun(''glbSolve'',Prob); ' ...
);

disp(Z)
fprintf('\n');
pause(1)

% ---------------------------------------------------------------------
function glb1Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Shekel 5 with TOMLAB Quick format\n');
fprintf('=====================================================\n');

Name   = 'Shekel 5';
x_L    = [ 0  0  0  0]';
x_U    = [10 10 10 10]';

% Generate the problem structure using the TOMLAB Quick format (short call)
Prob   = glcAssign('glb1_f', x_L, x_U, Name);

Result = glbSolve(Prob);

PrintResult(Result);

% ---------------------------------------------------------------------
function glb2Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Shekel 5 with TOMLAB Init File format\n');
fprintf('Shekel 5 is #1 in glb_prob.m\n');
fprintf('=====================================================\n');

% Generate the problem structure using the TOMLAB Init File format
Prob   = probInit('glb_prob', 1);

Result = glbSolve(Prob);

PrintResult(Result);

% ---------------------------------------------------------------------
function glb3Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Shekel 5 with TQ format and 6 restarts\n');
fprintf('Only five iterations per call\n');
fprintf('=====================================================\n');

Name   = 'Shekel 5';
x_L    = [ 0  0  0  0]';
x_U    = [10 10 10 10]';

% Generate the problem structure using the TOMLAB Quick format (short call)
Prob   = glcAssign('glb1_f', x_L, x_U, Name);

% Solve problem with 6 restarts, and only 5 iterations / call

Prob.optParam.MaxIter = 5;
Result                = glbSolve(Prob);

PrintResult(Result);
pause(1)

Prob.WarmStart = 1;
for i = 1:6
    Result = glbSolve(Prob);
    PrintResult(Result);
    pause(1)
end

% ---------------------------------------------------------------------
function glb4Demo
% ---------------------------------------------------------------------

format compact
fprintf('=====================================================\n');
fprintf('Solve Shekel 5 with TQ format, sending info to function\n');
fprintf('=====================================================\n');

Name   = 'Shekel 5';
x_L    = [ 0  0  0  0]';
x_U    = [10 10 10 10]';

% Generate the problem structure using the TOMLAB Quick format (short call)
Prob   = glcAssign('glb4_f', x_L, x_U, Name);

% Add information to be sent to glb4_f. Used in f(x) computation
Prob.user.A = [4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7]';
Prob.user.c = [.1 .2 .2 .4 .4]';

% For illustration, use multi-solver driver routine tomRun, instead of
% direct call to glbSolve. Set print level to 2, gives call to PrintResult

Result = tomRun('glbSolve',Prob,2);
