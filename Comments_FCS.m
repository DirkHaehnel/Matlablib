function varargout = Comments_FCS(varargin)
% COMMENTS_FCS M-file for Comments_FCS.fig
%      COMMENTS_FCS, by itself, creates a new COMMENTS_FCS or raises the existing
%      singleton*.
%
%      H = COMMENTS_FCS returns the handle to a new COMMENTS_FCS or the handle to
%      the existing singleton*.
%
%      COMMENTS_FCS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMMENTS_FCS.M with the given input arguments.
%
%      COMMENTS_FCS('Property','Value',...) creates a new COMMENTS_FCS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Comments_FCS_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Comments_FCS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Comments_FCS

% Last Modified by GUIDE v2.5 22-Feb-2006 15:35:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Comments_FCS_OpeningFcn, ...
                   'gui_OutputFcn',  @Comments_FCS_OutputFcn, ...
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


% --- Executes just before Comments_FCS is made visible.
function Comments_FCS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Comments_FCS (see VARARGIN)

% Choose default command line output for Comments_FCS
handles.output = hObject;
handles.res=varargin{1};
handles.head = varargin{2};
handles.t3r_name = varargin{3};
load ('c:\programme\matlab7\work\parameters_ini.mat')
handles.parameters=parameters;
set(handles.text_t3r_name,'String',handles.t3r_name);
set(handles.text3,'String','Pinhole Diameter [µm]:')
set(handles.text_corr,'String','Correction Ring [µm]:')
set(handles.text7,'String','Power [µW/Laser]:')
set(handles.text12,'String',handles.head.FileTime)
set(handles.edit_objective,'String',parameters.Objective);
set(handles.edit_corr,'String',parameters.CorrectionRing);
set(handles.edit_power,'String',parameters.Power);
set(handles.text_acq_time,'String',['',time2str(sec2hr(handles.head.NCounts/handles.head.TTTRCDFRate))],'HorizontalAlignment','left');
set(handles.edit_temperature,'String',parameters.Temperature);
set(handles.edit_pinholesize,'String',parameters.PinholeSize);
set(handles.edit_sample,'String',parameters.Sample);
set(handles.edit_comments,'String',parameters.Comments);
set(handles.edit_reprate,'String',parameters.RepRate);
set(handles.edit_coverslide,'String',parameters.CoverSlideThickness);
set(handles.edit_experiment,'String',parameters.Experiment);
set(handles.edit_piezopos,'String',parameters.PiezoPos);
set(handles.edit_beamwaistradius,'String',handles.parameters.BeamWaistRadius);
set(handles.edit_tubelenspos,'String',handles.parameters.TubeLensPos);
set(handles.edit_focpos,'String',handles.parameters.FocPos);

guidata(hObject, handles);


uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Comments_FCS_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.res;
varargout{2} = handles.head;
varargout{3} = handles.parameters;
delete(handles.figure1)
% --- Executes on button press in pb_objective_01.
function pb_objective_01_Callback(hObject, eventdata, handles)

h=get(handles.pb_objective_01,'String');
set(handles.edit_objective,'String',h);
guidata(hObject, handles);
edit_objective_Callback(hObject, eventdata, handles)

function edit_objective_Callback(hObject, eventdata, handles)

edit_corr_Callback(hObject, eventdata, handles)

function edit_corr_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_objective_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_finishandsave.
function pb_finishandsave_Callback(hObject, eventdata, handles)
handles.parameters.Objective = get(handles.edit_objective,'String');
handles.parameters.CorrectionRing = get(handles.edit_corr,'String');
handles.parameters.AcqTime=get(handles.text_acq_time,'String');
handles.parameters.PinholeSize = get(handles.edit_pinholesize,'String');
handles.parameters.Power = get(handles.edit_power,'String');
handles.parameters.Temperature = get(handles.edit_temperature,'String');
handles.parameters.Sample = get(handles.edit_sample,'String');
handles.parameters.RepRate = get(handles.edit_reprate,'String');
handles.parameters.Comments = get(handles.edit_comments,'String');
handles.parameters.Experiment = get(handles.edit_experiment,'String');
handles.parameters.PiezoPos = get(handles.edit_piezopos,'String');
handles.parameters.BeamWaistRadius = get(handles.edit_beamwaistradius,'String');
handles.parameters.TubeLensPos = get(handles.edit_tubelenspos,'String');
handles.parameters.CoverSlideThickness = get(handles.edit_coverslide,'String');
handles.parameters.T3R_Name = get(handles.text_t3r_name,'String');
handles.parameters.FocPos = get(handles.edit_focpos,'String');
name =handles.parameters.T3R_Name;
name = [name(1:end-3),'mat'];
[FileName,PathName] = uiputfile('*.mat','Save as',name);
if FileName==0  
   h = warndlg('Data not saved!');
   guidata(hObject, handles);
   uiresume
else
    handles.currentfile=FileName;
    handles.currentdir=PathName;
    directory=[handles.currentdir,'\'];
    cd (directory)
    res=handles.res;
    head=handles.head;
    parameters=handles.parameters;
    save (handles.currentfile,'res','head','parameters')
    save ('c:\programme\matlab7\work\parameters_ini.mat','parameters')
end
guidata(hObject, handles);
uiresume

% --- Executes on button press in pb_finish.
function pb_finish_Callback(hObject, eventdata, handles)

button = questdlg('Quit without saving?','Quit','Yes');
switch button
    case 'Yes'
guidata(hObject, handles);
uiresume
    otherwise
guidata(hObject, handles);
end



% --- Executes during object creation, after setting all properties.
function edit_corr_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)


function edit_pinholesize_Callback(hObject, eventdata, handles)

function edit_pinholesize_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_power_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_power_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_reprate_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_reprate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_reprate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sample_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sample as text
%        str2double(get(hObject,'String')) returns contents of edit_sample as a double


% --- Executes during object creation, after setting all properties.
function edit_sample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_comments_Callback(hObject, eventdata, handles)
% hObject    handle to edit_comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_comments as text
%        str2double(get(hObject,'String')) returns contents of edit_comments as a double


% --- Executes during object creation, after setting all properties.
function edit_comments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_temperature_Callback(hObject, eventdata, handles)
% hObject    handle to edit_temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_temperature as text
%        str2double(get(hObject,'String')) returns contents of edit_temperature as a double


% --- Executes during object creation, after setting all properties.
function edit_temperature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_beamwaistradius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beamwaistradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beamwaistradius as text
%        str2double(get(hObject,'String')) returns contents of edit_beamwaistradius as a double


% --- Executes during object creation, after setting all properties.
function edit_beamwaistradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beamwaistradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tubelenspos_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tubelenspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tubelenspos as text
%        str2double(get(hObject,'String')) returns contents of edit_tubelenspos as a double


% --- Executes during object creation, after setting all properties.
function edit_tubelenspos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tubelenspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_piezopos_Callback(hObject, eventdata, handles)
% hObject    handle to edit_piezopos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_piezopos as text
%        str2double(get(hObject,'String')) returns contents of edit_piezopos as a double


% --- Executes during object creation, after setting all properties.
function edit_piezopos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_piezopos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to edit_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_experiment as text
%        str2double(get(hObject,'String')) returns contents of edit_experiment as a double


% --- Executes during object creation, after setting all properties.
function edit_experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_coverslide_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coverslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coverslide as text
%        str2double(get(hObject,'String')) returns contents of edit_coverslide as a double


% --- Executes during object creation, after setting all properties.
function edit_coverslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coverslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_focpos_Callback(hObject, eventdata, handles)
% hObject    handle to edit_focpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_focpos as text
%        str2double(get(hObject,'String')) returns contents of edit_focpos as a double


% --- Executes during object creation, after setting all properties.
function edit_focpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_focpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


