function varargout = QDControl(varargin)
% QDCONTROL M-file for QDControl.fig
%      QDCONTROL, by itself, creates a new QDCONTROL or raises the existing
%      singleton*.
%
%      H = QDCONTROL returns the handle to a new QDCONTROL or the handle to
%      the existing singleton*.
%
%      QDCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QDCONTROL.M with the given input arguments.
%
%      QDCONTROL('Property','Value',...) creates a new QDCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QDControl_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QDControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QDControl

% Last Modified by GUIDE v2.5 22-Feb-2005 12:49:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QDControl_OpeningFcn, ...
                   'gui_OutputFcn',  @QDControl_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before QDControl is made visible.
function QDControl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QDControl (see VARARGIN)

% Choose default command line output for QDControl
handles.output = hObject;

% UIWAIT makes QDControl wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.n0 = 1.52;
handles.n1 = 1.33;
handles.lambda = 0.57;

nn = 20;
handles.pix = 16;%6.45;
[x,y] = meshgrid(-nn:nn,-nn:nn);
handles.p = angle(x+i*y);
handles.r = sqrt(x.^2+y.^2);

handles.kappa = str2double(get(handles.kappa_edit,'String'));
handles.ratio = str2double(get(handles.ratio_edit,'String'));
handles.omega = str2double(get(handles.omega_edit,'String'));
handles.theta = str2double(get(handles.theta_edit,'String'));
handles.phi = 0;
if nargin 
    pos = find(strcmpi('kappa', varargin));
    if ~isempty(pos)
        handles.kappa = varargin{pos+1};
        set(handles.kappa_edit,'String',num2str(handles.kappa));
    end
    pos = find(strcmpi('ratio', varargin));
    if ~isempty(pos)
        handles.ratio = varargin{pos+1};
        set(handles.ratio_edit,'String',num2str(handles.ratio));
    end
    pos = find(strcmpi('omega', varargin));
    if ~isempty(pos)
        handles.omega = varargin{pos+1};
        set(handles.omega_edit,'String',num2str(handles.omega));
    end
    pos = find(strcmpi('theta', varargin));
    if ~isempty(pos)
        handles.theta = varargin{pos+1};
        set(handles.theta_edit,'String',num2str(handles.theta));
    end
    pos = find(strcmpi('phi', varargin));
    if ~isempty(pos)
        handles.phi = varargin{pos+1};
    end
end

set(handles.theta_slider,'sliderstep',1/90*[1 10],'max',90,'min',0,'Value',handles.theta);
set(handles.omega_slider,'sliderstep',1/180*[1 10],'max',180,'min',0,'Value',handles.omega);
set(handles.kappa_slider,'sliderstep',0.01*[1 5],'max',1,'min',-1,'Value',handles.kappa);
set(handles.ratio_slider,'sliderstep',0.02*[1 5],'max',1,'min',0,'Value',handles.ratio);
set(handles.phi_slider,'sliderstep',1/360*[1 10],'max',360,'min',0,'Value',handles.phi);

handles.mag = str2double(get(handles.mag_edit,'String'));
handles.focus = str2double(get(handles.foc_edit,'String'));
handles.na = str2double(get(handles.na_edit,'String'));
if nargin
    pos = find(strcmpi('mag', varargin));
    if ~isempty(pos)
        handles.mag = varargin{pos+1};
        set(handles.mag_edit,'String',num2str(handles.mag));
    end
    pos = find(strcmpi('focus', varargin));
    if ~isempty(pos)
        handles.focus = varargin{pos+1};
        set(handles.foc_edit,'String',num2str(handles.focus));
    end
    pos = find(strcmpi('na', varargin));
    if ~isempty(pos)
        handles.na = varargin{pos+1};
        set(handles.na_edit,'String',num2str(handles.na));
    end
end

[tmp, tmp, tmp, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole([0 ceil(size(handles.r,1)/sqrt(2)/handles.mag*handles.pix)],0,handles.na,handles.n0,handles.n1,handles.n1,[],0,[],handles.lambda,handles.mag,handles.focus);

f00 = fxx0.*conj(byx0);
f01 = -fxx0.*conj(byz);
f02 = -fxx0.*conj(byx2);
f10 = -fxz.*conj(byx0);
f11 = fxz.*conj(byz);
f12 = fxz.*conj(byx2);
f20 = -fxx2.*conj(byx0);
f21 = fxx2.*conj(byz);
f22 = fxx2.*conj(byx2);
handles.rho = rho/handles.pix;
handles.f00 = f00;
handles.f01 = f01;
handles.f02 = f02;
handles.f10 = f10;
handles.f11 = f11;
handles.f12 = f12;
handles.f20 = f20;
handles.f21 = f21;
handles.f22 = f22;

handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));

if nargin
    pos = find(strcmpi('image',varargin));
    if ~isempty(pos)
        axes(handles.imm);
        mim(varargin{pos+1});
    end
end
axes(handles.pic);
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = QDControl_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes during object creation, after setting all properties.
function na_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to na_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject, 'String', '1.2');

function na_edit_Callback(hObject, eventdata, handles)
% hObject    handle to na_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of na_edit as text
%        str2double(get(hObject,'String')) returns contents of na_edit as a double
handles.na = str2double(get(hObject,'String'));
[tmp, tmp, tmp, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole([0 ceil(size(handles.r,1)/sqrt(2)/handles.mag*handles.pix)],0,handles.na,handles.n0,handles.n1,handles.n1,[],0,[],handles.lambda,handles.mag,handles.focus);
f00 = fxx0.*conj(byx0);
f01 = -fxx0.*conj(byz);
f02 = -fxx0.*conj(byx2);
f10 = -fxz.*conj(byx0);
f11 = fxz.*conj(byz);
f12 = fxz.*conj(byx2);
f20 = -fxx2.*conj(byx0);
f21 = fxx2.*conj(byz);
f22 = fxx2.*conj(byx2);
handles.f00 = f00;
handles.f01 = f01;
handles.f02 = f02;
handles.f10 = f10;
handles.f11 = f11;
handles.f12 = f12;
handles.f20 = f20;
handles.f21 = f21;
handles.f22 = f22;
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mag_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mag_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', '110');


function mag_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mag_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mag_edit as text
%        str2double(get(hObject,'String')) returns contents of mag_edit as a double
handles.mag = str2double(get(hObject,'String'));
[tmp, tmp, tmp, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole([0 ceil(size(handles.r,1)/sqrt(2)/handles.mag*handles.pix)],0,handles.na,handles.n0,handles.n1,handles.n1,[],0,[],handles.lambda,handles.mag,handles.focus);
f00 = fxx0.*conj(byx0);
f01 = -fxx0.*conj(byz);
f02 = -fxx0.*conj(byx2);
f10 = -fxz.*conj(byx0);
f11 = fxz.*conj(byz);
f12 = fxz.*conj(byx2);
f20 = -fxx2.*conj(byx0);
f21 = fxx2.*conj(byz);
f22 = fxx2.*conj(byx2);
handles.rho = rho/handles.pix;
handles.f00 = f00;
handles.f01 = f01;
handles.f02 = f02;
handles.f10 = f10;
handles.f11 = f11;
handles.f12 = f12;
handles.f20 = f20;
handles.f21 = f21;
handles.f22 = f22;
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function foc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', '1.2');

function foc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to foc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of foc_edit as text
%        str2double(get(hObject,'String')) returns contents of foc_edit as a double
handles.focus = str2double(get(hObject,'String'));
[tmp, tmp, tmp, rho, tmp, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole([0 ceil(size(handles.r,1)/sqrt(2)/handles.mag*handles.pix)],0,handles.na,handles.n0,handles.n1,handles.n1,[],0,[],handles.lambda,handles.mag,handles.focus);
f00 = fxx0.*conj(byx0);
f01 = -fxx0.*conj(byz);
f02 = -fxx0.*conj(byx2);
f10 = -fxz.*conj(byx0);
f11 = fxz.*conj(byz);
f12 = fxz.*conj(byx2);
f20 = -fxx2.*conj(byx0);
f21 = fxx2.*conj(byz);
f22 = fxx2.*conj(byx2);
handles.f00 = f00;
handles.f01 = f01;
handles.f02 = f02;
handles.f10 = f10;
handles.f11 = f11;
handles.f12 = f12;
handles.f20 = f20;
handles.f21 = f21;
handles.f22 = f22;
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function omega_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', '0');

function omega_edit_Callback(hObject, eventdata, handles)
% hObject    handle to omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of omega_edit as text
%        str2double(get(hObject,'String')) returns contents of omega_edit as a double
handles.omega = str2double(get(hObject,'String'));
set(handles.omega_slider,'Value',handles.omega);
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function theta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', '0');

function theta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to theta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta_edit as text
%        str2double(get(hObject,'String')) returns contents of theta_edit as a double
handles.theta = str2double(get(hObject,'String'));
set(handles.theta_slider,'Value',handles.theta);
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function kappa_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kappa_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', '0');

function kappa_edit_Callback(hObject, eventdata, handles)
% hObject    handle to kappa_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kappa_edit as text
%        str2double(get(hObject,'String')) returns contents of kappa_edit as a double
handles.kappa = str2double(get(hObject,'String'));
set(handles.kappa_slider,'Value',handles.kappa);
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ratio_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject, 'String', '0');

function ratio_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ratio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ratio_edit as text
%        str2double(get(hObject,'String')) returns contents of ratio_edit as a double
handles.ratio = str2double(get(hObject,'String'));
set(handles.ratio_slider,'Value',handles.ratio);
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate pic


% --- Executes during object creation, after setting all properties.
function theta_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function theta_slider_Callback(hObject, eventdata, handles)
% hObject    handle to theta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.theta = get(hObject,'Value');
set(handles.theta_edit,'String',num2str(handles.theta));
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function omega_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function omega_slider_Callback(hObject, eventdata, handles)
% hObject    handle to omega_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.omega = get(hObject,'Value');
set(handles.omega_edit,'String',num2str(handles.omega));
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function kappa_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kappa_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function kappa_slider_Callback(hObject, eventdata, handles)
% hObject    handle to kappa_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.kappa = get(hObject,'Value');
set(handles.kappa_edit,'String',num2str(handles.kappa));
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ratio_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratio_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function ratio_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ratio_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.ratio = get(hObject,'Value');
set(handles.ratio_edit,'String',num2str(handles.ratio));
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function phi_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function phi_slider_Callback(hObject, eventdata, handles)
% hObject    handle to phi_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.phi = get(hObject,'Value');
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, handles.phi, get(handles.radiobutton1,'Value'));
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function radiobutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Value',1);


% --- Executes during object creation, after setting all properties.
function radiobutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Value',0);


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
set(handles.radiobutton2,'Value',0);
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, ...
    handles.phi, 1);
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
set(handles.radiobutton1,'Value',0);
handles.int = int_compute(handles.theta, handles.omega, handles.kappa, handles.ratio, handles.rho, handles.r, handles.p, ...
    handles.f00, handles.f01, handles.f02, handles.f10, handles.f11, handles.f12, handles.f20, handles.f21, handles.f22, ...
    handles.phi, 0);
mim(handles.int); title(int2str(round(handles.phi)),'fontsize',12)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function imm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate imm
return

function int = int_compute(theta, omega, kappa, ratio, rho, r, p, f00, f01, f02, f10, f11, f12, f20, f21, f22, phi, flag);

theta = theta/180*pi;
omega = omega/180*pi;
p = p - phi/180*pi;

if flag==1
    psi = 0;
    intc1(1,:) = 3*f00+3*f22+f11 - f00*cos(2*psi)*(1-cos(2*theta)+kappa*(3+cos(2*theta))*cos(2*omega)) + ...
        (f00-f11+f22)*(cos(2*theta)+2*kappa*cos(2*omega)*sin(theta)^2) + 4*kappa*f00*cos(theta)*sin(2*psi)*sin(2*omega);
    intc1(2,:) = -2*(2*f01-f21+2*f10-f12)*sin(theta)*(cos(psi)*cos(theta)*(1-kappa*cos(2*omega))+kappa*sin(psi)*sin(2*omega));
    intc1(3,:) = -3*f20+f11-3*f02 + (f20+f02)*cos(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)+2*sin(theta)^2) - ...
        (f20+f11+f02)*(cos(2*theta)+2*kappa*cos(2*omega)*sin(theta)^2) - 4*(f20+f02)*kappa*cos(theta)*sin(2*psi)*sin(2*omega);
    intc1(4,:) = 2*(f21+f12)*sin(theta)*(cos(psi)*cos(theta)*(1-kappa*cos(2*omega))+kappa*sin(psi)*sin(2*omega));
    intc1(5,:) = -f22*(cos(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)+2*sin(theta)^2)-4*kappa*cos(theta)*sin(2*psi)*sin(2*omega));
    ints1(1,:) = 2*(f21+f12)*sin(theta)*(cos(theta)*(1-kappa*cos(2*omega))*sin(psi)-kappa*cos(psi)*sin(2*omega));
    ints1(2,:) = (f20+f02)*(sin(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)+2*sin(theta)^2)+4*kappa*cos(2*psi)*cos(theta)*sin(2*omega));
    ints1(3,:) = ints1(1,:);
    ints1(4,:) = -f22*(sin(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)-2*sin(theta)^2)+4*kappa*cos(2*psi)*cos(theta)*sin(2*omega));
    
    psi = pi/2;
    intc2(1,:) = 3*f00+3*f22+f11 - f00*cos(2*psi)*(1-cos(2*theta)+kappa*(3+cos(2*theta))*cos(2*omega)) + ...
        (f00-f11+f22)*(cos(2*theta)+2*kappa*cos(2*omega)*sin(theta)^2) + 4*kappa*f00*cos(theta)*sin(2*psi)*sin(2*omega);
    intc2(2,:) = -2*(2*f01-f21+2*f10-f12)*sin(theta)*(cos(psi)*cos(theta)*(1-kappa*cos(2*omega))+kappa*sin(psi)*sin(2*omega));
    intc2(3,:) = -3*f20+f11-3*f02 + (f20+f02)*cos(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)+2*sin(theta)^2) - ...
        (f20+f11+f02)*(cos(2*theta)+2*kappa*cos(2*omega)*sin(theta)^2) - 4*(f20+f02)*kappa*cos(theta)*sin(2*psi)*sin(2*omega);
    intc2(4,:) = 2*(f21+f12)*sin(theta)*(cos(psi)*cos(theta)*(1-kappa*cos(2*omega))+kappa*sin(psi)*sin(2*omega));
    intc2(5,:) = -f22*(cos(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)+2*sin(theta)^2)-4*kappa*cos(theta)*sin(2*psi)*sin(2*omega));
    ints2(1,:) = 2*(f21+f12)*sin(theta)*(cos(theta)*(1-kappa*cos(2*omega))*sin(psi)-kappa*cos(psi)*sin(2*omega));
    ints2(2,:) = (f20+f02)*(sin(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)+2*sin(theta)^2)+4*kappa*cos(2*psi)*cos(theta)*sin(2*omega));
    ints2(3,:) = ints2(1,:);
    ints2(4,:) = -f22*(sin(2*psi)*(kappa*(3+cos(2*theta))*cos(2*omega)-2*sin(theta)^2)+4*kappa*cos(2*psi)*cos(theta)*sin(2*omega));
    
    % int = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
    int = interp1(rho,real(intc1(1,:)),r,'cubic') + interp1(rho,real(intc2(1,:)),r,'cubic');
    for jj=1:4 
        int = int + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
            sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
        int = int + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
            sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
    end
    
    intc1(1,:) = 4*(f00+f22)*sin(theta)^2+4*f11*cos(theta)^2;
    intc1(2,:) = -2*(f01-f21+f10-f12)*sin(2*theta);
    intc1(3,:) = -4*(f20+f02)*sin(theta)^2;
    
    tmp = interp1(rho,real(intc1(1,:)),r,'cubic');
    for jj=1:2 
        tmp = tmp + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic'); 
    end
else
    intc1(1,:) = (f00+f22)*(cos(omega)^2*cos(theta)^2 + sin(omega)^2) + f11*cos(omega)^2*sin(theta)^2;
    intc1(2,:) = -0.5*(f01-f21+f10-f12)*sin(2*theta)*cos(omega)^2;
    intc1(3,:) = 0.25*(f20+f02)*(1-3*cos(2*omega)-2*cos(omega)^2*cos(2*theta));
    ints1(1,:) = -0.5*(f01-f21+f10-f12)*sin(2*omega)*sin(theta);
    ints1(2,:) = -(f20+f02)*sin(2*omega)*cos(theta);
    
    omega = omega + pi/2;
    intc2(1,:) = (f00+f22)*(cos(omega)^2*cos(theta)^2 + sin(omega)^2) + f11*cos(omega)^2*sin(theta)^2;
    intc2(2,:) = -0.5*(f01-f21+f10-f12)*sin(2*theta)*cos(omega)^2;
    intc2(3,:) = 0.25*(f20+f02)*(1-3*cos(2*omega)-2*cos(omega)^2*cos(2*theta));
    ints2(1,:) = -0.5*(f01-f21+f10-f12)*sin(2*omega)*sin(theta);
    ints2(2,:) = -(f20+f02)*sin(2*omega)*cos(theta);

    int = (1-kappa)/2*interp1(rho,real(intc1(1,:)),r,'cubic') + (1+kappa)/2*interp1(rho,real(intc2(1,:)),r,'cubic');
    for jj=1:2
        int = int + (1-kappa)/2*cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
            (1-kappa)/2*sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
        int = int + (1+kappa)/2*cos(jj*p).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
            (1+kappa)/2*sin(jj*p).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
    end
    
    intc1(1,:) = (f00+f22)*sin(theta)^2+f11*cos(theta)^2;
    intc1(2,:) = -0.5*(f01-f21+f10-f12)*sin(2*theta);
    intc1(3,:) = -(f20+f02)*sin(theta)^2;
    
    tmp = interp1(rho,real(intc1(1,:)),r,'cubic');
    for jj=1:2 
        tmp = tmp + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic'); 
    end
    
end

int = (1-ratio)*int + ratio*tmp;
return

% --- Executes during object deletion, before destroying properties.
function pic_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to pic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function imm_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to imm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


