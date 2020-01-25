%翻转--------------------------------------------------------------------
function varargout = mflip(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mflip_OpeningFcn, ...
                   'gui_OutputFcn',  @mflip_OutputFcn, ...
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
function mflip_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = mflip_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
global im;
axes(handles.axes1);
imshow(im);

function pushbutton1_Callback(hObject, eventdata, handles)
global im;
axes(handles.axes1);
I=im; 
[ROW,COL,DIM] = size(I);
Ih = uint8(zeros(ROW, COL,DIM));%Horizontal mirroring
%水平镜像
for i =1:ROW
    for j=1:COL
        for k=1:DIM
        x = i;
        y = COL-j+1;
        z = k;
        Ih(x,y,z) =I(i,j,k);
        end
    end
end
imshow(Ih);
function pushbutton3_Callback(hObject, eventdata, handles)
global im;
axes(handles.axes1);
I=im; 
imshow(I);
[ROW,COL,DIM] = size(I);
Ihv = uint8(zeros(ROW, COL,DIM));
%水平垂直镜像
for i =1:ROW
    for j=1:COL
        for k=1:DIM
        x = ROW-i+1;
        y = COL-j+1;
        z = k;
        Ihv(x,y,z) =I(i,j,k);
        end
    end
end
imshow(Ihv);

function pushbutton4_Callback(hObject, eventdata, handles)
global im;
axes(handles.axes1);
I=im; 
[ROW,COL,DIM] = size(I);
Iv = uint8(zeros(ROW, COL,DIM));%Vertical mirroring
%垂直镜像
for i =1:ROW
    for j=1:COL
        for k=1:DIM
        x = ROW-i+1;
        y = j;
        z = k;
        Iv(x,y,z) =I(i,j,k);
        end
    end
end
imshow(Iv);
