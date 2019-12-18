function varargout = AGIAS(varargin)
% AGIAS MATLAB code for AGIAS.fig
%      AGIAS, by itself, creates a new AGIAS or raises the existing
%      singleton*.
%
%      H = AGIAS returns the handle to a new AGIAS or the handle to
%      the existing singleton*.
%
%      AGIAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGIAS.M with the given input arguments.
%
%      AGIAS('Property','Value',...) creates a new AGIAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AGIAS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AGIAS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE'editarea Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AGIAS

% Last Modified by GUIDE v2.5 18-May-2018 17:42:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AGIAS_OpeningFcn, ...
    'gui_OutputFcn',  @AGIAS_OutputFcn, ...
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


% --- Executes just before AGIAS is made visible.
function AGIAS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AGIAS (see VARARGIN)

% Choose default command line output for AGIAS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(hObject,'toolbar','figure')
% UIWAIT makes AGIAS wait for user response (see UIRESUME)
% uiwait(handles.figure1);

axes(handles.logo)
logo=imread('logo.png');
imshow(logo)

try
    datain=evalin('base','AGIAS_SAVE');
    dataout=datain;
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
    view([-30 10])
catch
    
end

% --- Outputs from this function are returned to the command line.
function varargout = AGIAS_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in surfacebox.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       AGIAS Aircraft Generation, Iteration, and Analysis Systems       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIST BOXES %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function surfacebox_Callback(hObject, eventdata, handles)
% hObject    handle to surfacebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns surfacebox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from surfacebox
i=get(handles.surfacebox,'Value');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;
set(handles.editname,'String',surfaces{i}.name)
set(handles.editarea,'String',num2str(surfaces{i}.s))
set(handles.editspan,'String',num2str(surfaces{i}.b))
set(handles.edittaper,'String',num2str(surfaces{i}.taper))
set(handles.editx,'String',num2str(surfaces{i}.coord(1)))
set(handles.edity,'String',num2str(surfaces{i}.coord(2)))
set(handles.editz,'String',num2str(surfaces{i}.coord(3)))
set(handles.editsweep,'String',num2str(surfaces{i}.sweep))
if ~isempty(surfaces{i}.surfaces)
   for j=1:length(surfaces{i}.surfaces)
      names{j}=surfaces{i}.surfaces(j).name;
   end
   set(handles.controlbox,'String',names)
   set(handles.controlbox,'Value',1)
   controlbox_Callback(handles.controlbox,eventdata,handles)
else
   set(handles.controlbox,'String','')
   set(handles.cnameedit,'String','N/A')
set(handles.ccordedit,'String','N/A')
set(handles.cbiedit,'String','N/A')
set(handles.cbfedit,'String','N/A')
set(handles.csedit,'String','N/A')
end

if strcmp(class(surfaces{i}),'Main_Wing')
    set(handles.editdihedral,'String',num2str(surfaces{i}.dihedral))
    set(handles.editairfoilout,'String',surfaces{i}.airfoil_out)
    set(handles.editairfoilin,'String',surfaces{i}.airfoil_in)
elseif strcmp(class(surfaces{i}),'Hor_Stab')
    set(handles.editdihedral,'String',num2str(surfaces{i}.dihedral))
    set(handles.editairfoilout,'String','N/A')
    set(handles.editairfoilin,'String',surfaces{i}.airfoil)
else
    set(handles.editdihedral,'String','N/A')
    set(handles.editairfoilout,'String','N/A')
    set(handles.editairfoilin,'String',surfaces{i}.airfoil)
end

function surfacebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to surfacebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

try
    datain=evalin('base','AGIAS_SAVE');
    surfaces=datain.surfaces;
    for i=1:numel(datain.surfaces)
        names{i}=datain.surfaces{i}.name;
    end
    set(hObject,'String',names)
catch
    list={'wing','htail','vtail'};
    set(hObject,'String',list)
    surfaces{1}=Main_Wing(5,5,1,0,[],[],0,[],[0,0,0],'wing');
    surfaces{2}=Hor_Stab(2,1,1,0,[],0,[],[2.5,0,0],'htail');
    surfaces{3}=Ver_Stab(1,.5,1,0,[],[],[2.5,0,0],'vtail');
    param.v=10;
    param.g=9.81;
    param.d=1.225;
    param.beta=0;
    param.alpha=0;
    param.use_alpha=0;
    aircraft=Aircraft(surfaces{1},surfaces{2},surfaces{3},[],[0,0,0],1.5,[1 1 1],'plane',param);
    dataout.surfaces=surfaces;
    dataout.aircraft=aircraft;
    dataout.parameters.v=10;
    dataout.parameters.g=9.81;
    dataout.parameters.d=1.225;
    assignin('base','AGIAS_SAVE',dataout)
end

function sectionbox_Callback(hObject, eventdata, handles)

function sectionbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sectionbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function controlbox_Callback(hObject, eventdata, handles)

j=get(hObject,'Value');
i=get(handles.surfacebox,'Value');
names=get(hObject,'String');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;
control=surfaces{i}.surfaces(j);
set(handles.cnameedit,'String',control.name)
set(handles.ccordedit,'String',num2str(control.chinge))
set(handles.cbiedit,'String',num2str(control.bi))
set(handles.cbfedit,'String',num2str(control.be))
set(handles.csedit,'String',num2str(control.sign))

function controlbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controlbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Edit Surfaces %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

function editarea_Callback(hObject, eventdata, handles)
% hObject    handle to editarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
text=get(hObject,'String');
try
    s=str2num(text);
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    surfaces=datain.surfaces;
    current=surfaces{i};
    current.s=s;
    surfaces{i}=current;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
catch
    error('Enter a Number')
end

function editarea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.s))
catch
    set(hObject,'String',num2str(1))
end

function editspan_Callback(hObject, eventdata, handles)
% hObject    handle to editspan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editspan as text
%        str2num(get(hObject,'String')) returns contents of editspan as a double
text=get(hObject,'String');
b=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
current=surfaces{i};
current.b=b;
surfaces{i}=current;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editspan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editspan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.b))
catch
    set(hObject,'String',num2str(1))
end

function edittaper_Callback(hObject, eventdata, handles)
% hObject    handle to edittaper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittaper as text
%        str2num(get(hObject,'String')) returns contents of edittaper as a double
text=get(hObject,'String');
taper=str2num(text);
if taper < 0
    errordlg('taper ratio must be positive')
    set(hObject,'String',1)
else
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    surfaces=datain.surfaces;
    current=surfaces{i};
    current.taper=taper;
    surfaces{i}=current;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
end

function edittaper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittaper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.taper))
catch
    set(hObject,'String',num2str(1))
end

function editdihedral_Callback(hObject, eventdata, handles)
% hObject    handle to editdihedral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdihedral as text
%        str2num(get(hObject,'String')) returns contents of editdihedral as a double

text=get(hObject,'String');

dihedral=str2num(text);
if dihedral < -pi/2 ||  dihedral > pi/2
    errordlg('dihedral must be between -pi/2 and pi/2 radians')
    set(hObject,'String',0)
else
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    surfaces=datain.surfaces;
    
    if strcmp(class(surfaces{i}),'Ver_Stab')
        disp('hello')
        set(hObject,'String','N\A')
    else
        
        surfaces=datain.surfaces;
        current=surfaces{i};
        current.dihedral=dihedral;
        surfaces{i}=current;
        dataout=datain;
        dataout.surfaces=surfaces;
        wings=[];
        htails=[];
        vtails=[];
        for j=1:numel(surfaces)
            if strcmp(class(surfaces{1,j}),'Main_Wing')
                wings=[wings,surfaces{j}];
            elseif class(surfaces{j})=='Hor_Stab'
                htails=[htails,surfaces{j}];
            elseif class(surfaces{j})=='Ver_Stab'
                vtails=[vtails,surfaces{j}];
            end
        end
        dataout.aircraft.wings=wings;
        dataout.aircraft.htails=htails;
        dataout.aircraft.vtails=vtails;
        assignin('base','AGIAS_SAVE',dataout)
        axes(handles.plotplane)
        cla
        aircraft=dataout.aircraft;
        aircraft.plot_all;
    end
end

function editdihedral_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editdihedral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    if strcmp(class(datain.surfaces{1}),'Main_Wing')
        set(hObject,'String',num2str(datain.surfaces{1}.dihedral))
    else
        set(hObject,'String','N/A')
    end
catch
    set(hObject,'String',num2str(0))
end

function editsweep_Callback(hObject, eventdata, handles)
% hObject    handle to editsweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editsweep as text
%        str2double(get(hObject,'String')) returns contents of editsweep as a double
text=get(hObject,'String');
sweep=str2num(text);
if sweep > pi/2 || sweep < - pi/2
    errordlg('sweep must be between -pi/2 and pi/2 radians')
    set(hObject,'String',0)
else
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    surfaces=datain.surfaces;
    current=surfaces{i};
    current.sweep=sweep;
    surfaces{i}=current;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
end

function editsweep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editsweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.sweep))
catch
    set(hObject,'String',num2str(0))
end

function editx_Callback(hObject, eventdata, handles)
% hObject    handle to editx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editx as text
%        str2num(get(hObject,'String')) returns contents of editx as a double
text=get(hObject,'String');

xcord=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
current=surfaces{i};
current.coord(1)=xcord;
surfaces{i}=current;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.coord(1)))
catch
    set(hObject,'String',num2str(0))
end

function edity_Callback(hObject, eventdata, handles)
% hObject    handle to edity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edity as text
%        str2num(get(hObject,'String')) returns contents of edity as a double
text=get(hObject,'String');

ycord=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
current=surfaces{i};
current.coord(2)=ycord;
surfaces{i}=current;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function edity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.coord(2)))
catch
    set(hObject,'String',num2str(0))
end

function editz_Callback(hObject, eventdata, handles)
% hObject    handle to editz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editz as text
%        str2num(get(hObject,'String')) returns contents of editz as a double
text=get(hObject,'String');
zcord=str2num(text); %#ok<*ST2NM>
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
current=surfaces{i};
current.coord(3)=zcord;
surfaces{i}=current;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.surfaces{1}.coord(3)))
catch
    set(hObject,'String',num2str(0))
end

function editairfoilin_Callback(hObject, eventdata, handles)
% hObject    handle to editairfoilin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editairfoilin as text
%        str2num(get(hObject,'String')) returns contents of editairfoilin as a double
text=get(hObject,'String');
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
current=surfaces{i};
current.airfoil_in=text;
surfaces{i}=current;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editairfoilin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editairfoilin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.surfaces{1}.airfoil_in)
catch
    set(hObject,'String',[])
end

function editairfoilout_Callback(hObject, eventdata, handles)
% hObject    handle to editairfoilout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editairfoilout as text
%        str2num(get(hObject,'String')) returns contents of editairfoilout as a double
text=get(hObject,'String');
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
if ~strcmp(class(surfaces{i}),'Main_Wing')
    set(hObject,'String','N/A')
else
    surfaces=datain.surfaces;
    current=surfaces{i};
    current.airfoil_out=text;
    surfaces{i}=current;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
end

function editairfoilout_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to editairfoilout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    if strcmp(class(datain.surfaces{1}),'Main_Wing')
        set(hObject,'String',datain.surfaces{1}.airfoil_out)
    else
        set(hObject,'String','N/A')
    end
catch
    set(hObject,'String',[])
end

function editname_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to editname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editname as text
%        str2num(get(hObject,'String')) returns contents of editname as a double
text=get(hObject,'String');
if isempty(text)
    text='xxx';
end
set(hObject,'String',text)
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
surfaces=datain.surfaces;
current=surfaces{i};
current.name=text;
surfaces{i}=current;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    names{j}=surfaces{j}.name;  %#ok<*AGROW>
end
set(handles.surfacebox,'String',names)
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing') %#ok<*STISA>
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.surfaces{1}.name)
catch
    set(hObject,'String','wing')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Edit Aircraft %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

function editcgx_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
cgx=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
dataout=datain;
dataout.aircraft.cg(1)=cgx;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editcgx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcgx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.aircraft.cg(1)))
catch
    set(hObject,'String',num2str(1))
end

function editcgy_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
cgy=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
dataout=datain;
dataout.aircraft.cg(2)=cgy;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editcgy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcgy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.aircraft.cg(2)))
catch
    set(hObject,'String',num2str(1))
end

function editcgz_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
cgz=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
dataout=datain;
dataout.aircraft.cg(3)=cgz;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editcgz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcgz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.aircraft.cg(3)))
catch
    set(hObject,'String',num2str(1))
end

function editixx_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
ixx=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
dataout=datain;
dataout.aircraft.inert(1)=ixx;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editixx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editixx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.aircraft.inert(1)))
catch
    set(hObject,'String',num2str(1))
end

function editiyy_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
iyy=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
dataout=datain;
dataout.aircraft.inert(2)=iyy;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editiyy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editiyy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.aircraft.inert(2)))
catch
    set(hObject,'String',num2str(1))
end

function editiz_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
izz=str2num(text);
datain=evalin('base','AGIAS_SAVE');
i=get(handles.surfacebox,'Value');
dataout=datain;
dataout.aircraft.inert(3)=izz;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editiz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editiz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',num2str(datain.aircraft.inert(3)))
catch
    set(hObject,'String',num2str(1))
end

function editv_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
v=str2num(text);
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
dataout.aircraft.param.v=v;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.aircraft.param.v)
catch
    set(hObject,'String',num2str(10))
end

function editd_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
d=str2num(text);
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
dataout.aircraft.param.d=d;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.aircraft.param.d)
catch
    set(hObject,'String',num2str(1.225))
end

function editg_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
g=str2num(text);
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
dataout.aircraft.param.g=g;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.aircraft.param.g)
catch
    set(hObject,'String',num2str(9.81))
end

function editm_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
m=str2num(text);
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
dataout.aircraft.m=m;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.aircraft.m)
catch
    set(hObject,'String',num2str(1.5))
end

function editalpha_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
alpha=str2num(text);
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
dataout.aircraft.param.alpha=alpha;
dataout.aircraft.param.use_alpha=1;
set(handles.setalphabutton,'Value',1)
set(handles.getalphabutton,'Value',0)
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editalpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editalpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    if datain.aircraft.param.use_alpha==1
        set(hObject,'String',datain.aircraft.param.alpha)
    else
        set(hObject,'String','N/A')
    end
catch
    set(hObject,'String','N/A')
end

function editbeta_Callback(hObject, eventdata, handles)
text=get(hObject,'String');
beta=str2num(text);
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
dataout.aircraft.param.beta=beta;
assignin('base','AGIAS_SAVE',dataout)
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function editbeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editbeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    datain=evalin('base','AGIAS_SAVE');
    set(hObject,'String',datain.aircraft.param.beta)
catch
    set(hObject,'String',num2str(0))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit Control %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function cnameedit_Callback(hObject, eventdata, handles)

    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    j=get(handles.controlbox,'Value');
    name=get(hObject,'String');
    surfaces=datain.surfaces;
    controls=surfaces{i}.surfaces;
    current=controls(j);
    current.name=name;
    controls(j)=current;
    for k=1:length(controls)
       names{k}=controls(k).name; 
    end
    set(handles.controlbox,'String',names)
    surfaces{i}.surfaces=controls;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;

function cnameedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cnameedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ccordedit_Callback(hObject, eventdata, handles)
    str=get(hObject,'String');
    chord=str2num(str);
    if chord >= 0 && chord <= 1
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    j=get(handles.controlbox,'Value');
    surfaces=datain.surfaces;
    controls=surfaces{i}.surfaces;
    current=controls(j);
    current.chinge=chord;
    controls(j)=current;
    surfaces{i}.surfaces=controls;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
    else
        set(hObject,'String','0.25')
        errordlg('Chord must be between 0 and 1')
    end

function ccordedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ccordedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cbiedit_Callback(hObject, eventdata, handles)

    str=get(hObject,'String');
    bi=str2num(str);
    if bi >= 0 && bi <= 1
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    j=get(handles.controlbox,'Value');
    surfaces=datain.surfaces;
    controls=surfaces{i}.surfaces;
    current=controls(j);
    current.bi=bi;
    controls(j)=current;
    surfaces{i}.surfaces=controls;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
    else
        set(hObject,'String','0.1')
        errordlg('bi must be between 0 and 1')
    end

function cbiedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cbiedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cbfedit_Callback(hObject, eventdata, handles)

    str=get(hObject,'String');
    bf=str2num(str);
    if bf >= 0 && bf <= 1
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    j=get(handles.controlbox,'Value');
    surfaces=datain.surfaces;
    controls=surfaces{i}.surfaces;
    current=controls(j);
    current.be=bf;
    controls(j)=current;
    surfaces{i}.surfaces=controls;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
    else
        set(hObject,'String','0.9')
        errordlg('bf must be between 0 and 1')
    end

function cbfedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cbfedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function csedit_Callback(hObject, eventdata, handles)

    str=get(hObject,'String');
    s=str2num(str);
    if s==1 || s == -1
    datain=evalin('base','AGIAS_SAVE');
    i=get(handles.surfacebox,'Value');
    j=get(handles.controlbox,'Value');
    surfaces=datain.surfaces;
    controls=surfaces{i}.surfaces;
    current=controls(j);
    current.sign=s;
    controls(j)=current;
    surfaces{i}.surfaces=controls;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing')
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
    axes(handles.plotplane)
    cla
    aircraft=dataout.aircraft;
    aircraft.plot_all;
    else
        set(hObject,'String','0.25')
        errordlg('sign must be 1 or -1')
    end

function csedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to csedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PUSH BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

function pushairfoilin_Callback(hObject, eventdata, handles)
% hObject    handle to pushairfoilin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=get(handles.surfacebox,'Value');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;
[FileName,~]=uigetfile('*.dat','Airfoil File');
set(handles.editairfoilin,'String',FileName);
surfaces{i}.airfoil_in=FileName;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing') %#ok<*STISA>
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)

function pushairfoilout_Callback(hObject, eventdata, handles)
% hObject    handle to pushairfoilout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=get(handles.surfacebox,'Value');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;

if strcmp(class(surfaces{i}),'Main_Wing') %#ok<*STISA>
    [FileName,~]=uigetfile('*.dat','Airfoil File');
    set(handles.editairfoilout,'String',FileName);
    surfaces{i}.airfoil_out=FileName;
    dataout=datain;
    dataout.surfaces=surfaces;
    wings=[];
    htails=[];
    vtails=[];
    for j=1:numel(surfaces)
        if strcmp(class(surfaces{1,j}),'Main_Wing') %#ok<*STISA>
            wings=[wings,surfaces{j}];
        elseif class(surfaces{j})=='Hor_Stab'
            htails=[htails,surfaces{j}];
        elseif class(surfaces{j})=='Ver_Stab'
            vtails=[vtails,surfaces{j}];
        end
    end
    dataout.aircraft.wings=wings;
    dataout.aircraft.htails=htails;
    dataout.aircraft.vtails=vtails;
    assignin('base','AGIAS_SAVE',dataout)
else
    errordlg('Stabalizers Can only Have 1 Airfoil')
end

function runavlbutton_Callback(hObject, eventdata, handles)
datain=evalin('base','AGIAS_SAVE');
aircraft=datain.aircraft;
set(handles.statusbutton,'String','Building')
aircraft.build_file
aircraft.build_mass
set(handles.statusbutton,'String','Running')
aircraft=aircraft.run_avl;
axes(handles.plotplane)
cla
aircraft.plot_all
dataout=datain;
dataout.aircraft=aircraft;
assignin('base','AGIAS_SAVE',dataout)
CL=aircraft.run.CLtot;
CLstring=sprintf('CL=%1.4f',CL);
set(handles.cltext,'String',CLstring);
CD=aircraft.run.CDtot;
CDstring=sprintf('CD=%1.4f',CD);
set(handles.cdtext,'String',CDstring);
roll=sprintf('z=%.3f   w=%.3f hz',aircraft.eig.roll(1),aircraft.eig.roll(2)/(2*pi));
dutch=sprintf('z=%.3f   w=%.3f hz',aircraft.eig.dutch(1),aircraft.eig.dutch(2)/(2*pi));
short=sprintf('z=%.3f   w=%.3f hz',aircraft.eig.short(1),aircraft.eig.short(2)/(2*pi));
spiral=sprintf('z=%.3f   w=%.3f hz',aircraft.eig.spiral(1),aircraft.eig.spiral(2)/(2*pi));
phugoid=sprintf('z=%.3f   w=%.3f hz',aircraft.eig.phugoid(1),aircraft.eig.phugoid(2)/(2*pi));
set(handles.rolltext,'String',roll)
set(handles.dutchrolltext,'String',dutch)
set(handles.shorttext,'String',short)
set(handles.spiraltext,'String',spiral)
set(handles.phugoidtext,'String',phugoid)
alpha=sprintf('alpha=%.2f',aircraft.run.alpha);
Cla=sprintf('CLa=%.2f',aircraft.st.CLa);
clb=sprintf('Clb=%.2f',aircraft.st.Clb);
cma=sprintf('Cma=%.2f',aircraft.st.Cma);
cnb=sprintf('Cnb=%.2f',aircraft.st.Cnb);
clp=sprintf('Clp=%.2f',aircraft.st.Clp);
cmq=sprintf('Cmq=%.2f',aircraft.st.Cmq);
cnr=sprintf('Cnr=%.2f',aircraft.st.Cnr);
i=1;
while ~strcmp('Main_Wing',class(datain.surfaces{i}))
    i=i+1;
end
if i > numel(datain.surfaces)
    i=1;
end
chord=datain.surfaces{i}.mean_chord;
static_margin=(-aircraft.cg(1)+aircraft.st.NP)/chord*100;
static_text=sprintf('Static Margin=%.2f %%',static_margin);
set(handles.atext,'String',alpha)
set(handles.clatext,'String',Cla)
set(handles.clbtext,'String',clb)
set(handles.cmatext,'String',cma)
set(handles.cnbtext,'String',cnb)
set(handles.clptext,'String',clp)
set(handles.cmqtext,'String',cmq)
set(handles.cnrtext,'String',cnr)
set(handles.staticmargintext,'String',static_text)
set(handles.statusbutton,'String','Done!')

function pushnewsurface_Callback(hObject, eventdata, handles)
% hObject    handle to pushnewsurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=get(handles.surfacebox,'Value');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;
type=menu('What type of surface?','Main Wing','Horizontal Stabalizer','Vertical Stabalizer');
l=numel(surfaces);
if type==1
    surfaces{l+1}=Main_Wing(1,1,0,0,[],[],0,[],[0,0,0],'new_wing');
elseif type==2
    surfaces{l+1}=Hor_Stab(1,1,0,0,[],0,[],[1,0,0],'new_htail');
elseif type==3
    surfaces{l+1}=Ver_Stab(1,1,0,0,[],[],[1,0,0],'new_vtail');
end

for i=1:numel(surfaces)
    names{i}=surfaces{i}.name;
end
set(handles.surfacebox,'String',names)

dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)

function pushdeletesurface_Callback(hObject, eventdata, handles)
% hObject    handle to pushdeletesurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=get(handles.surfacebox,'Value');
datain=evalin('base','AGIAS_SAVE');
len=numel(datain.surfaces);
surfaces=datain.surfaces([1:i-1,i+1:len]);
for i=1:numel(surfaces)
    names{i}=surfaces{i}.name;
end
set(handles.surfacebox,'String',names)
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
set(handles.surfacebox,'Value',1')
axes(handles.plotplane)
cla
aircraft=dataout.aircraft;
aircraft.plot_all;

function setalphabutton_Callback(hObject, eventdata, handles)
i=get(hObject,'Value');
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
if i==1
    set(hObject,'Value',1)
    set(handles.getalphabutton,'Value',0)
    datain=evalin('base','AGIAS_SAVE');
    set(handles.editalpha,'String',datain.aircraft.param.alpha)
    dataout.aircraft.param.use_alpha=1;
else
    set(hObject,'Value',0)
    set(handles.getalphabutton,'Value',1)
    set(handles.editalpha,'String','N/A')
    dataout.aircraft.param.use_alpha=0;
end
assignin('base','AGIAS_SAVE',dataout)

function setalphabutton_CreateFcn(hObject, eventdata, handles)
try
    datain=evalin('base','AGIAS_SAVE');
    if datain.aircraft.param.use_alpha==1
        set(hObject,'Value',1)
    else
        set(hObject,'Value',0)
    end
catch
    set(hObject,'Value',0)
end

function getalphabutton_Callback(hObject, eventdata, handles)
i=get(hObject,'Value');
datain=evalin('base','AGIAS_SAVE');
dataout=datain;
if i==1
    set(hObject,'Value',1)
    set(handles.setalphabutton,'Value',0)
    set(handles.editalpha,'String','N/A')
    dataout.aircraft.param.use_alpha=0;
else
    set(hObject,'Value',0)
    set(handles.setalphabutton,'Value',1)
    set(handles.editalpha,'String',dataout.aircraft.param.alpha)
    dataout.aircraft.param.use_alpha=1;
end
assignin('base','AGIAS_SAVE',dataout)

function getalphabutton_CreateFcn(hObject, eventdata, handles)
try
    datain=evalin('base','AGIAS_SAVE');
    if datain.aircraft.param.use_alpha==1
        set(hObject,'Value',0)
    else
        set(hObject,'Value',1)
    end
catch
    set(hObject,'Value',1)
end

function newcontrolbutton_Callback(hObject, eventdata, handles)

i=get(handles.surfacebox,'Value');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;
controls=surfaces{i}.surfaces;
if isempty(controls)
controls=Surface(.25,1,.2,.8,'new_surface');
else
controls(length(controls)+1)=Surface(.25,1,.2,.8,'new_surface');
end
for k=1:length(controls)
    names{k}=controls(k).name;
end
set(handles.controlbox,'String',names)
surfaces{i}.surfaces=controls;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)

function deletecontrolbutton_Callback(hObject, eventdata, handles)
i=get(handles.surfacebox,'Value');
j=get(handles.controlbox,'Value');
datain=evalin('base','AGIAS_SAVE');
surfaces=datain.surfaces;
controls=surfaces{i}.surfaces;
len=length(controls);
controls=controls([1:j-1,j+1:len]);
if len-1 > 0
for k=1:length(controls)
    names{k}=controls(k).name;
end
else 
    names='';
end
set(handles.controlbox,'String',names)
surfaces{i}.surfaces=controls;
dataout=datain;
dataout.surfaces=surfaces;
wings=[];
htails=[];
vtails=[];
for j=1:numel(surfaces)
    if strcmp(class(surfaces{1,j}),'Main_Wing')
        wings=[wings,surfaces{j}];
    elseif class(surfaces{j})=='Hor_Stab'
        htails=[htails,surfaces{j}];
    elseif class(surfaces{j})=='Ver_Stab'
        vtails=[vtails,surfaces{j}];
    end
end
dataout.aircraft.wings=wings;
dataout.aircraft.htails=htails;
dataout.aircraft.vtails=vtails;
assignin('base','AGIAS_SAVE',dataout)
set(handles.controlbox,'Value',1)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE WINDOWS %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function logo_CreateFcn(hObject, eventdata, handles)

function plotplane_CreateFcn(hObject, eventdata, handles)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEXT BOXES %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function cltext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','CL=x.xxx')

function cdtext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','CD=x.xxx')

function rolltext_CreateFcn(hObject, eventdata, handles)
text=sprintf('z=x.xxx   w=x.xxx hz');
set(hObject,'String',text)

function phugoidtext_CreateFcn(hObject, eventdata, handles)
text=sprintf('z=x.xxx   w=x.xxx hz');
set(hObject,'String',text)

function spiraltext_CreateFcn(hObject, eventdata, handles)
text=sprintf('z=x.xxx   w=x.xxx hz');
set(hObject,'String',text)

function shorttext_CreateFcn(hObject, eventdata, handles)
text=sprintf('z=x.xxx   w=x.xxx hz');
set(hObject,'String',text)

function dutchrolltext_CreateFcn(hObject, eventdata, handles)
text=sprintf('z=x.xxx   w=x.xxx hz');
set(hObject,'String',text)

function atext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','alpha=xx.xx')

function clatext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','CLa=xx.xxx')

function clbtext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Clb=xx.xxx')

function cmatext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Cma=xx.xxx')

function cnbtext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Cnb=xx.xxx')

function clptext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Clp=xx.xxx')

function cmqtext_CreateFcn(hObject, eventdata, handles)

function cnrtext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Cnr=xx.xxx')

function staticmargintext_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Static Margin=xx.xx%%')

function statusbutton_Callback(hObject, eventdata, handles)

function statusbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String','...')
