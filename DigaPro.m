function varargout = DigaPro(varargin)
% DIGAPRO M-file for DigaPro.fig
%      DIGAPRO, by itself, creates a new DIGAPRO or raises the existing
%      singleton*.
%
%      H = DIGAPRO returns the handle to a new DIGAPRO or the handle to
%      the existing singleton*.
%
%      DIGAPRO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGAPRO.M with the given input arguments.
%
%      DIGAPRO('Property','Value',...) creates a new DIGAPRO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DigaPro_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DigaPro_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DigaPro

% Last Modified by GUIDE v2.5 09-Mar-2004 08:13:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DigaPro_OpeningFcn, ...
                   'gui_OutputFcn',  @DigaPro_OutputFcn, ...
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


% --- Executes just before DigaPro is made visible.
function DigaPro_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DigaPro (see VARARGIN)

% Choose default command line output for DigaPro
handles.output = hObject;
handles.action = 'Find Molecules';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DigaPro wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DigaPro_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function tsh_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsh_slider (see GCBO)
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
function tsh_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tsh_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

switch handles.action
    case 'Find Molecules'
        handles.tsh = get(hObject,'Value');
        handles.act_im = handles.cim>handles.tsh*sqrt(handles.err);
        mim(handles.act_im);
        guidata(hObject, handles);
    case 'Threshold Reset'
        set(hObject,'Value',0);
    case 'Continue'
        handles.tsh = sum(handles.len<round(get(hObject,'Value')*max(handles.len)));
        handles.act_im = handles.cl>handles.tsh;
        mim(handles.act_im);
        guidata(hObject, handles);
end

% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if exist('DigaPro.ini','file')==2
    load -MAT 'DigaPro.ini' path;
else
    path = '';
end

[name, path] = uigetfile('*.tif','',path);

if ~(name==0)
    
    if exist('DigaPro.ini','file')==2
        save 'DigaPro.ini' path -MAT -APPEND
    else
        save 'DigaPro.ini' path -MAT
    end
    
    name = [path name];
    handles.imlen = length(imfinfo(name));
    im = double(imread(name,1)); 
    for j=2:handles.imlen 
        tmp = double(imread(name,j));
        im = im + tmp;
    end
    handles.im = im;
    handles.act_im = im;
    
    mim(handles.act_im);
    
end
guidata(hObject, handles);
    

% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1)

% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mim(handles.im)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mim(handles.act_im)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'String')
    case 'Find Molecules'
        im = handles.im;
        nn = 10;
        t = 2:nn;
        
        a2 = autocorr(im);
        a2 = sum(a2');
        a2 = simplex('exp2fun',1,0,[],[],[],2:nn,a2(2:nn));
        
        a1 = autocorr(im');
        a1 = sum(a1');
        a1 = simplex('exp2fun',1,0,[],[],[],2:nn,a1(2:nn));
        
        nn = ceil(mean([a1,a2]));
        n2 = 2*nn;
        [x,y] = meshgrid(-n2:n2,-n2:n2);
        mask = exp(-x.^2/2/a1^2-y.^2/2/a2^2);
        
        bck = disk(n2); 
        mask = mask.*bck;
        
        n1 = size(mask,1);
        n2 = size(mask,2);
        bck = bck/sqrt(sum(sum(bck.^2)));
        for j=1:size(mask,3) 
            mask(:,:,j) = mask(:,:,j)./sqrt(sum(sum(mask(:,:,j).^2)));
        end
        im0 = mconv2(im,bck);
        im02 = mconv2(im.^2,bck/max(bck(:)));
        err = inf*im;
        cim = im;
        sim = 1+0*im;
        for s = 1:size(mask,3)
            crs = sum(sum(bck.*mask(:,:,s)));
            crs = inv([1 crs; crs 1]);
            im1 = mconv2(im,squeeze(mask(:,:,s)));
            tmperr = im02 - crs(1,1)*im0.*im0 - 2*crs(1,2)*im0.*im1 - crs(2,2)*im1.*im1;
            tmpcim = crs(1,2)*im0 + crs(2,2)*im1;
            cim(tmperr<err) = tmpcim(tmperr<err);
            sim(tmperr<err) = s;
            err(tmperr<err) = tmperr(tmperr<err);
        end
        
        handles.tsh = 1;
        handles.sim = sim;
        handles.cim = cim;
        handles.err = err;
        handles.act_im = cim>handles.tsh*sqrt(err);
        mim(handles.act_im);
        set(hObject,'String','Continue')
        guidata(hObject, handles);
        
    case 'Continue'
        handles.action = 'Threshold Reset';
        DigaPro('tsh_slider_Callback',gcbo,[],guidata(gcbo));
        handles.action = 'Continue';
        [handles.cl,handles.len]=Cluster(handles.act_im);
        set(hObject,'String','Finish')
        handles.act_im = handles.cl>0;
        mim(handles.act_im);
        guidata(hObject, handles);
        
    case 'Finish'
        handles.action = 'Finish';
        [a,b] = meshgrid(1:size(handles.im,2),1:size(handles.im,1));
        cnt = 1;
        handles.len(1:handles.tsh) = [];
        handles.cl = handles.cl-handles.tsh; 
        handles.cl(handles.cl<0) = 0;
        for j=1:length(handles.len)
            ind = handles.cl==j;
            tmpxc = a(ind);
            tmpyc = b(ind);
            tmperr = handles.err(ind);
            tmps = handles.sim(ind);
            ind = tmperr==min(tmperr);
            sc(cnt) = tmps(ind);
            xc(cnt) = tmpxc(ind);
            yc(cnt) = tmpyc(ind);    
            cnt = cnt+1;
        end
        xc = round(xc);
        yc = round(yc);
%         ind = xc<n2+1 | xc>size(im,2)-n2 | yc<n1+1 | yc>size(im,1)-n1; 
%         xc(ind)=[]; yc(ind)=[]; sc(ind) = []; len(ind) = [];
        set(hObject,'String','Save Results')
        guidata(hObject, handles);
        
end


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


