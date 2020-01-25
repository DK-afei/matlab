%缩放--------------------------------------------------------------------
function varargout = enlarge_narrow(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @enlarge_narrow_OpeningFcn, ...
                   'gui_OutputFcn',  @enlarge_narrow_OutputFcn, ...
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
function enlarge_narrow_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = enlarge_narrow_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
axes(handles.axes1);
imshow(I1);
[m,n]=size(I1);
set(handles.axes1,'units','pixels');
pos=get(handles.axes1,'pos');
pos(3:4)=[m n];
set(handles.axes1,'pos',pos);
%缩放
function pushbutton1_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
axes(handles.axes1);
data1=get(handles.edit1,'string');
data1_num=str2num(data1);
a=I1;
mul=data1_num;
type=1;
%输入图像灰度值
%mul:缩放倍数
%type:1表示最邻近法 2位双线性插值法
%画出缩放后图像
[m,n]=size(a);
m1=m*mul;n1=n*mul;
if type==1
for i=1:m1
    for j=1:n1;
    b(i,j)=a(round(i/mul),round(j/mul));
    end
end
elseif type==2
    for i=1:m1-1
      for j=1:n1-1;
      u0=i/mul;v0=j/mul;
      u=round(u0);v=round(v0);
       s=u0-u;t=v0-v;
       b(i,j)=(a(u+1,v)-a(u,v))*s+(a(u,v+1)-a(u,v))*t+(a(u+1,v+1)+a(u,v)-a(u,v+1)-a(u+1,v))*s*t+a(u,v);
      end
    end
end 
%*****************************************************
b=uint8(b);
imshow(b);
title('处理后图像');
[m,n]=size(b);
set(handles.axes1,'units','pixels');
pos=get(handles.axes1,'pos');
pos(3:4)=[m n];
set(handles.axes1,'pos',pos);
function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
