%��ת--------------------------------------------------------------------
function varargout = rotate(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rotate_OpeningFcn, ...
                   'gui_OutputFcn',  @rotate_OutputFcn, ...
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
function rotate_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = rotate_OutputFcn(hObject, eventdata, handles) 
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
%��ת--------------------------------
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
img=I1;
degree=data1_num;
%��ȡͼƬ��Ϣ ע����ͨ����ȡ�� ��������������
[m,n,dep]=size(img);
%�������ת֮���γ�һ������εĳ��� ���Կ�Ч��ͼ
rm=round(m*abs(cosd(degree))+n*abs(sind(degree)));
rn=round(m*abs(sind(degree))+n*abs(cosd(degree)));
%����һ���¾�����ͨ���ģ��洢��ͼƬ����Ϣ
newimage=zeros(rm,rn,dep);
%����任 ������ 
m1=[1,0,0;0,1,0;-0.5*rm,-0.5*rn,1];
m2=[cosd(degree),sind(degree),0;-sind(degree),cosd(degree),0;0,0,1];
m3=[1,0,0;0,1,0;0.5*m,0.5*n,1];
%����ѭ������ÿһ�����ص���б任
for i=1:rm
    for j=1:rn
        tem=[i j 1];
        tem=tem*m1*m2*m3;
        x=tem(1,1);
        y=tem(1,2);
        x=round(x);
        y=round(y);
        if(x>0&&x<=m)&&(y>0&&y<=n)
        newimage(i,j,:)=img(x,y,:);
        end
        end
end
imshow(uint8(newimage));
