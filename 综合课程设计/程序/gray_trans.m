%灰度变换--------------------------------------------------------------------
function varargout = gray_trans(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gray_trans_OpeningFcn, ...
                   'gui_OutputFcn',  @gray_trans_OutputFcn, ...
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
function gray_trans_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = gray_trans_OutputFcn(hObject, eventdata, handles) 
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
%灰度变换------------------------------
function pushbutton1_Callback(hObject, eventdata, handles)
%线性变换
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
axes(handles.axes1);
data1=get(handles.edit1,'string');
data2=get(handles.edit2,'string');
data3=get(handles.edit3,'string');
data4=get(handles.edit4,'string');
data5=get(handles.edit5,'string');
data6=get(handles.edit6,'string');
data1_num=str2num(data1);
data2_num=str2num(data2);
data3_num=str2num(data3);
data4_num=str2num(data4);
data5_num=str2num(data5);
data6_num=str2num(data6);
I2=data1_num*I1+data2_num;
imshow(I2);
function text2_DeleteFcn(hObject, eventdata, handles)
function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton2_Callback(hObject, eventdata, handles)
%非线性变换log
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
axes(handles.axes1);
I2=log(double(I1)+1)*10;
imshow(uint8(I2));
title('log变换—log(I1+1)*10');
function pushbutton3_Callback(hObject, eventdata, handles)
%非线性变换幂律
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
axes(handles.axes1);
I_D=double(I1);
C=I_D/255;
I2=uint8(255*(C.^1.5));
imshow(I2);
title('幂律变换—γ=1.5');

function pushbutton4_Callback(hObject, eventdata, handles)
%灰度值映射
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
axes(handles.axes1);
data3=get(handles.edit3,'string');
data4=get(handles.edit4,'string');
data5=get(handles.edit5,'string');
data6=get(handles.edit6,'string');
data3_num=str2num(data3);
data4_num=str2num(data4);
data5_num=str2num(data5);
data6_num=str2num(data6);
I=I1;
Y=double(I1); %将参数I转为双精度浮点类型
[M,N]=size(Y);
for i=1:M 
    for j=1:N 
        if I(i,j)<data3_num 
            Y(i,j)=I(i,j); 
        else
            if I(i,j)<=data4_num
                Y(i,j)=(data6_num-data5_num)/(data4_num-data3_num)*(I(i,j)-data3_num)+data3_num; 
            else
                Y(i,j)=I(i,j);
            end
        end
    end
end
imshow(uint8(Y));

function edit3_Callback(hObject, eventdata, handles)
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit4_Callback(hObject, eventdata, handles)
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit5_Callback(hObject, eventdata, handles)
function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit6_Callback(hObject, eventdata, handles)
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

    

