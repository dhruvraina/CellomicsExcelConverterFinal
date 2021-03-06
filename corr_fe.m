%Correlation Function Front-End (Called from excelconv). Takes in input
%from user, and saves it. corr_bend reads input from file and performs
%plotting. All edits only in OpeningFunction and pbt_Go callback
% Author: Dhruv 
% Last Edit: 030515
% Usage: uiwait(corr_fe)

function varargout = corr_fe(varargin)
% CORR_FE M-file for corr_fe.fig
%      CORR_FE, by itself, creates a new CORR_FE or raises the existing
%      singleton*.
%
%      H = CORR_FE returns the handle to a new CORR_FE or the handle to
%      the existing singleton*.
%
%      CORR_FE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORR_FE.M with the given input arguments.
%
%      CORR_FE('Property','Value',...) creates a new CORR_FE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before corr_fe_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to corr_fe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help corr_fe

% Last Modified by GUIDE v2.5 30-Apr-2015 14:38:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @corr_fe_OpeningFcn, ...
                   'gui_OutputFcn',  @corr_fe_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before corr_fe is made visible.
function corr_fe_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to corr_fe (see VARARGIN)

% Choose default command line output for corr_fe
handles.output = hObject;
if exist('Corr_Vars.mat', 'file')==2
    load('Corr_Vars.mat');
    set(handles.rdbt_minnorm, 'Value', corr_vars(3));
    set(handles.rdbt_saveimg, 'Value', corr_vars(4));
    set(handles.etx_inA, 'String', num2str(corr_vars(1)+1));
    set(handles.etx_inB, 'String', num2str(corr_vars(2)+1));
else
    set(handles.rdbt_minnorm, 'Value', 1);
    set(handles.rdbt_saveimg, 'Value', 1);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes corr_fe wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = corr_fe_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbt_go.
function pbt_go_Callback(hObject, eventdata, handles)
% hObject    handle to pbt_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ColA = (str2double(get(handles.etx_inA,'String'))-1);
ColB = (str2double(get(handles.etx_inB,'String'))-1);
NormFlag = 1;
SaveFlag = get(handles.rdbt_saveimg, 'Value');
corr_vars = [ColA; ColB; NormFlag;SaveFlag];
save([cd '\Corr_Vars.mat'], 'corr_vars')
guidata(hObject, handles)
uiresume


% --- Executes on button press in rdbt_saveimg.
function rdbt_saveimg_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_saveimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_saveimg


% --- Executes on button press in rdbt_minnorm.
function rdbt_minnorm_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_minnorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_minnorm


% --- Executes on button press in rdbt_maxnorm.
function rdbt_maxnorm_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_maxnorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_maxnorm


% --- Executes on button press in rdbt_nonorm.
function rdbt_nonorm_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_nonorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_nonorm



function etx_inA_Callback(hObject, eventdata, handles)
% hObject    handle to etx_inA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_inA as text
%        str2double(get(hObject,'String')) returns contents of etx_inA as a double


% --- Executes during object creation, after setting all properties.
function etx_inA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_inA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_inB_Callback(hObject, eventdata, handles)
% hObject    handle to etx_inB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_inB as text
%        str2double(get(hObject,'String')) returns contents of etx_inB as a double


% --- Executes during object creation, after setting all properties.
function etx_inB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_inB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
