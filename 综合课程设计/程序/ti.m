%ƽ��--------------------------------------------------------------------
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
%ƽ��
init = I1; % ��ȡͼ��
[R, C] = size(init); % ��ȡͼ���С
res = zeros(R, C); % ����������ÿ�����ص�Ĭ�ϳ�ʼ��Ϊ0����ɫ��
delX = data1_num; % ƽ����X
delY = data2_num; % ƽ����Y
tras = [1 0 delX; 0 1 delY; 0 0 1]; % ƽ�Ƶı任���� 
for i = 1 : R
    for j = 1 : C
        temp = [i; j; 1];
        temp = tras * temp; % ����˷�
        x = temp(1, 1);
        y = temp(2, 1);
        % �任���λ���ж��Ƿ�Խ��
        if (x <= R) & (y <= C) & (x >= 1) & (y >= 1)
            res(x, y) = init(i, j);
        end
    end
end;
imshow(uint8(res)); % ��ʾͼ��
