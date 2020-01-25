%平移--------------------------------------------------------------------
function varargout = ti(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ti_OpeningFcn, ...
                   'gui_OutputFcn',  @ti_OutputFcn, ...
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
function ti_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = ti_OutputFcn(hObject, eventdata, handles) 
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
data2=get(handles.edit2,'string');
data1_num=str2num(data1);
data2_num=str2num(data2);
%%
%平移
init = I1; % 读取图像
[R, C] = size(init); % 获取图像大小
res = zeros(R, C); % 构造结果矩阵。每个像素点默认初始化为0（黑色）
delX = data1_num; % 平移量X
delY = data2_num; % 平移量Y
tras = [1 0 delX; 0 1 delY; 0 0 1]; % 平移的变换矩阵 
for i = 1 : R
    for j = 1 : C
        temp = [i; j; 1];
        temp = tras * temp; % 矩阵乘法
        x = temp(1, 1);
        y = temp(2, 1);
        % 变换后的位置判断是否越界
        if (x <= R) & (y <= C) & (x >= 1) & (y >= 1)
            res(x, y) = init(i, j);
        end
    end
end;
imshow(uint8(res)); % 显示图像
