%����ͼ����������--------------------------------------------------------------------
function varargout = ImageProcess(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageProcess_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageProcess_OutputFcn, ...
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
function ImageProcess_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.axes1,'visible','off');
set(handles.axes2,'visible','off');
guidata(hObject, handles);
function varargout = ImageProcess_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
function file_Callback(hObject, eventdata, handles)
function edit_Callback(hObject, eventdata, handles)
function filt_Callback(hObject, eventdata, handles)

%һ��ͼ����뱣�� ----------------------------------------------------------
% 1.ͼ���*****************************
function file_open_Callback(hObject, eventdata, handles)
axes(handles.axes1);
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'ѡ��ͼ��');
str=[pathname filename];
global im
im = imread(str);  
imshow(im);	
mysize=size(im);
%disp(mysize);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
[m,n]=size(I1);
str=[num2str(n),'x',num2str(m)];
set(handles.edit1,'String',str);
% 2.ͼ�񱣴�*****************************
function file_save_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.gif','GIF(*.gif)';...
                                 '*.*',  'All Files (*.*)'},...
                                 'Save Picture','Untitled');
if FileName==0
    return;
else
    h=getframe(handles.axes2);
    imwrite(h.cdata,[PathName,FileName]);
end;

%�������α任 ------------------------------------------------------------
function Geometric_transformation_Callback(hObject, eventdata, handles)
% 1.����*****************************
function enlarge_narrow_Callback(hObject, eventdata, handles)
enlarge_narrow;
% 2.ƽ��*****************************
function translation_Callback(hObject, eventdata, handles)
ti;
% 3.��ת*****************************
function rotate_Callback(hObject, eventdata, handles)
rotate;
% 4.��ת*****************************
function flip_Callback(hObject, eventdata, handles)
mflip;

%����ͼ��任 --------------------------------------------------------------
% 1.ͼ��ϳ�***********************
function synthesis_Callback(hObject, eventdata, handles)
%       aͼ���ں�
function blend_Callback(hObject, eventdata, handles)
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'ѡ��ͼ��');
str=[pathname filename];
I1 = imread(str);  
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'ѡ��ͼ��');
str=[pathname filename];
I2 = imread(str); 
[m n]=size(I1);
I3=zeros(m,n);
if(size(I1)==size(I2))
    I3=I1+I2;
    figure,
    subplot(131),imshow(I1);  
    subplot(132),imshow(I2);
    subplot(133),imshow(I3);
else
    sprintf('������ѡ��ͼƬ');
end
%       bֱ��ƴ��
function joint_Callback(hObject, eventdata, handles)
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'ѡ��ͼ��');
str=[pathname filename];
I1 = imread(str);  
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'ѡ��ͼ��');
str=[pathname filename];
I2 = imread(str); 
 img1=I1;
 img2=I2;
[H,W,k]=size(img1);
l_r=5;%�ص���ȣ�W-�� �� W��---�����������ƥ������ֱ��д�غ�����
L=W+1-l_r;%������
R=W;%�ұ�β��
n=R-L+1;%�ص���ȣ�����l_r
%ֱ��ƴ��ͼ
im=[img1,img2(:,n:W,:)];%1ȫͼ+2�ĺ��沿��
figure;imshow(im);title('ֱ��ƴ��ͼ');
% 2.�Ҷȱ任***********************
%       aͼ��ת
function reversion_Callback(hObject, eventdata, handles)
global im;
J=double(im);
J=-J+(256-1);                 %���Ա任
H=uint8(J);
axes(handles.axes2);
imshow(H);	
%       b�Ҷȱ任
function Linear_transform_Callback(hObject, eventdata, handles)
gray_trans;
% 3.ֱ��ͼ���⻯*******************
function histeq_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
f=I1;
[m,n]=size(f);
nf = m*n;%�����ظ���
h=zeros(1,256);%ÿ���Ҷȼ����ظ���
for i=1:m;
   for j=1:n;
       for k=1:256;
            if(f(i,j)==k)
              h(k)=h(k)+1;
            end
       end
   end;
end;
hs=zeros(1,256);
for l=1:256;
    hs(l)=h(l)/nf;%ÿ���Ҷȼ����ظ���ռ�ٷֱ�
end
%ֱ��ͼ���⻯
hp=zeros(1,256);
for i=1:256;
    if i==1
        hp(i)=hs(i);
    else
        hp(i)=hp(i-1)+hs(i);
    end
end
g=zeros(1,256);
for i=1:256;
    g(i)=hp(i)*255;
end
%�Ҷ�ֵӳ����g(i)�滻ԭͼ ���лҶ�ֵΪi�����ص�
for i = 1:m
    for j = 1: n
        f(i,j) = g(f(i,j));
    end
end
axes(handles.axes2);
imshow(f);%��ʾ���⻯���ͼ
title('���⻯���ͼ');

%�ġ�ͼ��ȥ�� --------------------------------------------------------------
% 1.�����˲�***********************
%      a��ֵ�˲�
function medfilt_Callback(hObject, eventdata, handles)
global im;
axes(handles.axes2);
data1=get(handles.edit7,'string');
data1_num=str2num(data1);
mysize=size(im);
if numel(mysize)>2
    I=rgb2gray(im);
else
    I=im;
end
Img=I;
masksize=data1_num;
exsize=floor(masksize/2);   %��������չ��С
Imgex=padarray(Img,[exsize,exsize],'replicate','both'); %��չͼƬ
[m,n]=size(Img);
Img_out=Img;    %��Img_out׼��Ϊ��Img��ͬ��size
for i=1:m
    for j=1:n
        neighbor=Imgex(i:i+masksize-1,j:j+masksize-1);  %��ȡ����
        Img_out(i,j)=median(neighbor(:));   %��ֵ�˲�
    end
end
axes(handles.axes2);
imshow(Img_out);
%          ��ֵ�˲�ģ������
function pushbutton2_Callback(hObject, eventdata, handles)
medfilt_Callback(hObject, eventdata, handles);
%      b��ֵ�˲�
function Mean_Callback(hObject, eventdata, handles)
global im;
data1=get(handles.edit6,'string');
data1_num=str2num(data1);
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
x=uint8(I1);
%x����Ҫ�˲���ͼ��,n��ģ���С(��n��n)  
n=data1_num;
a(1:n,1:n)=1;   %a��n��nģ��,Ԫ��ȫ��1  
[height, width]=size(x);   %����ͼ����hightxwidth��,��hight>n,width>n  
x1=double(x);  
x2=x1;  
for i=1:height-n+1  
	    for j=1:width-n+1  
	        c=x1(i:i+(n-1),j:j+(n-1)).*a; %ȡ��x1�д�(i,j)��ʼ��n��n��Ԫ����ģ�����  
        s=sum(sum(c));                 %��c�����и�Ԫ��֮��  
        x2(i+(n-1)/2,j+(n-1)/2)=s/(n*n); %����ģ�������ĸ�Ԫ�صľ�ֵ����ģ������λ�õ�Ԫ��  
    end  
end  
%δ����ֵ��Ԫ��ȡԭֵ  
imshow(uint8(x2));
%          ��ֵ�˲�ģ������
function pushbutton3_Callback(hObject, eventdata, handles)
Mean_Callback(hObject, eventdata, handles);
% 2.Ƶ���˲�***********************
%      a��ͨ�˲�
function high_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
data1=get(handles.edit8,'string');
data1_num=str2num(data1);
%% ��ͨ�˲�����ͼ��
moon=I1;
[m,n]=size(moon);
I=im2double(moon);
I_spectrum=fft2(I);
I_spectrum=fftshift(I_spectrum);
%% highpass filter
H=zeros(m,n);
centerx=m/2;
centery=n/2;
D0=data1_num;  % �ɵ���ͨ���뾶������ͨ���ĸ�Ƶ����
for x=1:m
    for y=1:n
        H(x,y)=exp(-((x-centerx)^2+(y-centery)^2)/(2*D0^2));  %�����˹�˲�ģ��
    end
end
H=1-H;

g1=H.*I_spectrum;% ��Ƶ����ͼ�񣬱�Եͼ��
g2=g1+I_spectrum;
% �ȷ����Ļ�����ת������ͼ��
g3=ifftshift(g2);
I2=real(ifft2(g3));
axes(handles.axes2);
imshow(real(ifft2(ifftshift(g1))),[]);
%          ��ͨ�˲���ֵ����
function pushbutton6_Callback(hObject, eventdata, handles)
high_Callback(hObject, eventdata, handles);
%      b��ͨ�˲�
function low_Callback(hObject, eventdata, handles)
%����Butterworth��ͨ�˲�
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
data1=get(handles.edit2,'string');
data1_num=str2num(data1);
J=double(I1);
f=fft2(J);
g=fftshift(f);
[M,N]=size(f);
n=3;d0=data1_num;
n1=floor(M/2);n2=floor(N/2);
for i=1:M
    for j=1:N
        d=sqrt((i-n1)^2+(j-n2)^2);
        h=1/(1+0.414*(d/d0)^(2*n));
        g(i,j)=h*g(i,j);
    end
end
g=ifftshift(g);
g=uint8(real(ifft2(g)));
axes(handles.axes2);
imshow(g);
%          ��ͨ�˲���ֵ���� 
function pushbutton5_Callback(hObject, eventdata, handles)
low_Callback(hObject, eventdata, handles);


% �塢ͼ���Ե���------------------------------------------------------------
function edge_detection_Callback(hObject, eventdata, handles)
% 1.Laplacian����ģ��
function Laplace_Callback(hObject, eventdata, handles)
%
l=[0 1 0 ;1 -4 1; 0 1 0];
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
[m n]=size(I1);
x=I1;
y=I1;
I1=im2double(I1);
P0=conv2(I1,l);
for i=1:m
    for j=1:n
        if(P0(i,j)>0.15)
            P0(i,j)=255;
        else 
            P0(i,j)=0;
        end
    end
end
axes(handles.axes2);
imshow(P0);
% 2.Prewitt����ģ��
function Prewitt_Callback(hObject, eventdata, handles)
p1=[-1 -1 -1;0 0 0;1 1 1];
p2=[-1 0 1; -1 0 1;-1 0 1];
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
[m n]=size(I1);
I1=im2double(I1);
P1=conv2(I1,p1);
P2=conv2(I1,p2);
P3=P1+P2;
for i=1:m
    for j=1:n
        if(P3(i,j)>0.30)
            P3(i,j)=255;
        else 
            P3(i,j)=0;
        end
    end
end
axes(handles.axes2);
imshow(P3);
% 3.sobel����ģ��
function Sobel_Callback(hObject, eventdata, handles)
s1=[-1 -2 -1;0 0 0;1 2 1];
s2=[-1 0 1;-2 0 2;-1 0 1];
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
I3=I1;
[m n]=size(I1);
I1=im2double(I1);
P3=conv2(I1,s1);
P4=conv2(I1,s2);
P5=P4+P3;
for i=1:m
    for j=1:n
        if(P5(i,j)>0.30)
            P5(i,j)=255;
        else 
            P5(i,j)=0;
        end
    end
end
axes(handles.axes2);
imshow(P5);
% 4.Roberts����ģ��
function Roberts_Callback(hObject, eventdata, handles)
r1=[-1 0;0 1];
r2=[0 -1;1 0];
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
[m n]=size(I1);
I1=im2double(I1);
P6=conv2(I1,r1);
P7=conv2(I1,r2);
P8=P6+P7;
for i=1:m
    for j=1:n
        if(P8(i,j)>0.10)
            P8(i,j)=255;
        else 
            P8(i,j)=0;
        end
    end
end
axes(handles.axes2);
imshow(P8);
% 5.�����������ģ��
function Log_Callback(hObject, eventdata, handles)
e=[-1 -1 -1;-1 8 -1;-1 -1 -1];
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
[m n]=size(I1);
I1=im2double(I1);
P9=conv2(I1,e);
for i=1:m
    for j=1:n
        if(P9(i,j)>0.30)
            P9(i,j)=255;
        else 
            P9(i,j)=0;
        end
    end
end
axes(handles.axes2);
imshow(P9);

%����ͼ����----------------------------------------------------------------
function sharpen_Callback(hObject, eventdata, handles)
% 1.Laplacian������
function laplace_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
%I1Ϊ�Ҷ�ͼ
I=im2double(I1);
KernelType=4;
c=1;
%��չ�����������
KernelSize=3;
len=floor(KernelSize/2);
%��ԭʼͼ�������չ���˴������˾�����չ��Ŀ���ǽ����Ե���������
f_pad=padarray(I,[len,len],'symmetric');
[M,N]=size(f_pad);
switch KernelType
    case -4
        L=[0 1 0;
            1 -4 1;
            0 1 0];
    case -8
        L=[1 1 1;
            1 -8 1;
            1 1 1];
    case 4
        L=[0 -1 0;
            -1 4 -1;
            0 -1 0];
    case 8 
        L=[-1 -1 -1;
            -1 8 -1;
            -1 -1 -1];
    %�����������Ǻϳ�������˹����
    case 5
        L=[0 -1 0;
            -1 5 -1;
            0 -1 0];     
    case 9 
        L=[-1 -1 -1;
            -1 9 -1;
            -1 -1 -1];        
end
if KernelType>0
    a=1;
else 
    a=-1;
end
for i=1+len:M-len
    for j=1+len:N-len
        %����չͼ���У�ȡ���ֲ�ͼ��
        Block=f_pad(i-len:i+len,j-len:j+len);
        %��������˹���ӵĽ��������ԭʼͼ�񣬵õ����ͼ��       
        g(i-len,j-len)=I(i-len,j-len)+ a*sum(sum(Block.*L));
        %����������˹���ӵ�������
        edge(i-len,j-len)=a*sum(sum(Block.*L));
    end
end
figure,imshow(edge);
axes(handles.axes2);
imshow(g);
% 2.������
function other_sharpen_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
size_1 = size(I1);
h = size_1(1);
w = size_1(2);
img_2 = repmat(uint8(0),h-2, w-2);
for i = 1:h-2
    for j = 1:w-2
        x = I1(i:i+2, j:j+2);
        a = sum(x(:));
        a_1 = int16(I1(i+1,j+1));
        a_1 = 10*a_1;
        a_2 = -a+a_1;
        img_2(i,j) = a_2;
    end
end
axes(handles.axes2);
imshow(img_2);

% �ߡ�ͼ��ָ�----------------------------------------------------------------
function image_segmentation_Callback(hObject, eventdata, handles)
% 1.��ֵ��********************
%        a�ֶ���ֵ
function Threshold_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
I=I1;
data1=get(handles.edit5,'string');
data1_num=str2num(data1);
%�˹�ѡ����ֵ���зָѡ����ֵΪ data1_num
[width,height]=size(I);
T1=data1_num;
for i=1:width
    for j=1:height
        if(I(i,j)<T1)
            BW1(i,j)=0;
        else 
            BW1(i,j)=1;
        end
    end
end
% figure;imshow(BW1),title('�˹���ֵ���зָ�');
%�Զ�ѡ����ֵ
% T2=graythresh(I);
% BW2=im2bw(I,T2);%Otus��ֵ���зָ�
% figure;imshow(BW2),title('Otus��ֵ���зָ�');
axes(handles.axes2);
imshow(BW1);
function pushbutton4_Callback(hObject, eventdata, handles)
Threshold_Callback(hObject, eventdata, handles);
%        b�Զ���ֵostu
function ostu_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
I=I1;
%�Զ�ѡ����ֵ
T2=graythresh(I);
BW2=im2bw(I,T2);%Otus��ֵ���зָ�
axes(handles.axes2);
imshow(BW2);
title('Otus��ֵ���зָ�');
% 2.��������
function Regional_growth_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
data3=get(handles.edit12,'string');
data3_num=str2num(data3);
I=I1;
if isinteger(I)
    I=im2double(I);
end
axes(handles.axes1);
imshow(I);
[M,N]=size(I);
 [y,x]=getpts;             %�������������ʼ��
  if x<0
      x=-x;
  end
  
  if y<0
      y=-y;
  end

x1=round(x);            %������ȡ��
y1=round(y);            %������ȡ��
try
seed=I(x1,y1);           %��������ʼ��Ҷ�ֵ����seed��
catch
end
J=zeros(M,N);          %��һ��ȫ����ԭͼ��ȴ��ͼ�����J����Ϊ���ͼ�����
J(x1,y1)=1;             %��J������ȡ�����Ӧλ�õĵ�����Ϊ��
sum=seed;              %��������������������ĵ�ĻҶ�ֵ�ĺ�
suit=1;                 %��������������������ĵ�ĸ���
count=1;               %��¼ÿ���ж�һ����Χ�˵�����������µ����Ŀ
threshold=data3_num;         %��ֵ��ע����Ҫ��double���ʹ洢��ͼ�������
while count>0
    s=0;                   %��¼�ж�һ����Χ�˵�ʱ�������������µ�ĻҶ�ֵ֮��
     count=0;
     for i=1:M
       for j=1:N
         if J(i,j)==1
          if (i-1)>0 && (i+1)<(M+1) && (j-1)>0 && (j+1)<(N+1)  %�жϴ˵��Ƿ�Ϊͼ��߽��ϵĵ�
           for u= -1:1                               %�жϵ���Χ�˵��Ƿ������ֵ����
            for v= -1:1
              if  J(i+u,j+v)==0 && abs(I(i+u,j+v)-seed)<=threshold&& 1/(1+1/15*abs(I(i+u,j+v)-seed))>0.8
                           J(i+u,j+v)=1;
                    %�ж��Ƿ���δ��ǣ�����Ϊ������ֵ�����ĵ�
                    %����������������������J����֮λ�ö�Ӧ�ĵ�����Ϊ��
                 count=count+1;
                 s=s+I(i+u,j+v);                      %�˵�ĻҶ�֮����s��
              end
            end
           end
          end
         end
       end
     end
    suit=suit+count;                                   %��n������ϵ�����������
    sum=sum+s;                                     %��s������ϵ�ĻҶ�ֵ�ܺ���
    seed=sum/suit;                                    %�����µĻҶ�ƽ��ֵ
end
axes(handles.axes2);
imshow(J);
%        ����������������ֵ
function pushbutton8_Callback(hObject, eventdata, handles)
Regional_growth_Callback(hObject, eventdata, handles);
% 3.���������ϲ�
function Regional_division_and_merger_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
% ʹ���������ͺϲ���ͼ��ָ�
f = I1;
g=splitmerge(f,2,@predicate);%2����ָ���������С�Ŀ飬predicate��������1��˵����Ҫ�ٷ��ѣ�����0˵������Ҫ��������
% figure,imshow(g);
% title('mindimΪ2ʱ�ķָ�ͼ��');
se=ones(8,8);
gdilate=imdilate(g,se);%������Ϊ�����ն�
% figure;imshow(gdilate);
% title('���ͺ��ͼ')
gerode=imerode(gdilate,se);%��ʴ��Ϊ������ԭ����С
% figure;imshow(gerode);
% title('��ʴ���ͼ')
axes(handles.axes2);
imshow(gerode);
%%%%%%%%%%
function g=splitmerge(f,mindim,fun)%f�Ǵ��ָ��ԭͼ��mindim�Ƕ���ֽ������������С�Ŀ飬������2������������
Q=2^nextpow2(max(size(f)));
[M,N]=size(f);
f=padarray(f,[Q-M,Q-N],'post');%�����ͼ���������顣f������ͼ������������ͼ���Ƚ�ͼ����䵽2���ݴ���ʹ����ķֽ����
%Ȼ��������������������post����ʾ��ÿһά�����һ��Ԫ�غ����,B = padarray(A,padsize,padval,direction)
%����padval����0���,Q��������ͼ��Ĵ�С��
S=qtdecomp(f,@split_test,mindim,fun);%S����split_test��qtdecomp divides a square image into four
% different sizes.S�ǰ����Ĳ����ṹ��ϡ����󣬴洢��ֵ�ǿ�Ĵ�С�����꣬��ϡ�������ʽ�洢
Lmax=full(max(S(:)));%����ϡ�����洢��ʽ�洢�ľ���任������ͨ����full matrix����ʽ�洢��full��sparseֻ�Ǵ洢��ʽ�Ĳ�ͬ
g=zeros(size(f));
MARKER=zeros(size(f));
for k=1:Lmax
    [vals,r,c]=qtgetblk(f,S,k);%vals��һ�����飬����f���Ĳ����ֽ��д�СΪk*k�Ŀ��ֵ����һ��k*k*�����ľ���
%������ָS���ж��ٸ�������С�Ŀ飬f�Ǳ��Ĳ����ֵ�ԭͼ��r��c�Ƕ�Ӧ�����Ͻǿ��������2*2�飬����������Ͻǿ�ʼ�������
        if ~isempty(vals)
            for I=1:length(r)
                    xlow=r(I);
                    ylow=c(I);
                    xhigh=xlow+k-1;
                    yhigh=ylow+k-1;
                    region=f(xlow:xhigh,ylow:yhigh);%�ҵ���Ӧ������
                    flag=feval(fun,region);%evaluates the function handle, fhandle,using arguments x1 through xn.ִ�к���fun��region�ǲ���
                    if flag%������ص���1������б��
                            g(xlow:xhigh,ylow:yhigh)=1;%Ȼ�󽫶�Ӧ��������1
                            MARKER(xlow,ylow)=1;%MARKER�����Ӧ�����Ͻ�������1
                    end
            end
        end
end	
	g=bwlabel(imreconstruct(MARKER,g));%imreconstructĬ��2Dͼ��8��ͨ���������������ϵ�����
	g=g(1:M,1:N);%����ԭͼ��Ĵ�С
%%%%%%%%%%
function v=split_test(B,mindim,fun)
 K=size(B,3);%B����qtdecomp�����������ģ�����ǰsize(B,3)���ص���B�Ĳ���������B�Ǽ�ά�ģ�����ʵ���Ͼ����м���B������С��ͼ���
%���������˼�Ǵ�qtdecomp������������B���ǵ�ǰ�ֽ�ɵ�K���m*m��ͼ��飬K��ʾ�ж��ٸ�������С��ͼ���
 v(1:K)=false;
   for I=1:K
        quadregion=B(:,:,I);
        if size(quadregion,1)<=mindim%����ֵĿ�Ĵ�СС��mindim��ֱ�ӽ���
                v(I)=false;
                continue
        end
        flag=feval(fun,quadregion);%quadregion��fun�����Ĳ���
        if flag%���flag��1��������Ҫ�ٷ�
            v(I)=true;%������൱��split_test����һ������predicate�����ã����صľ���ppredicate��ֵ
        end
   end
function flag=predicate(region)
sd=std2(region);
m=mean2(region);
flag=(sd>20)&(m>26)&(m<255);
%predicate��������Ŀ�ģ���split_test�б�����ʱ���ж��Ƿ�ñ��֣����������������ͷ���1����Ҫ�ٷ֣�����ͷ���0�������ٱ���
%�ڿ�ʼ�ϲ�ʱ�������������ж��Ƿ�ý��кϲ���ǣ�����1����ʾͨ�����ԣ��Ĳ�������1��䣬����0����û��ͨ�����ԣ���0���
%�������������ϲ��ķָ���У�predicate�������û��Զ���ģ��Ķ��������ؼ����غ��ָ�Ч���ĺû�����������Ŀ��ʱ���Էֱ��岻ͬ��
%����׼���Դﵽ��õ�Ч����

% �ˡ�ͼ��ѹ��----------------------------------------------------------------
function compression_Callback(hObject, eventdata, handles)
% 1.����ѹ�� 
function lossless_Callback(hObject, eventdata, handles)
%  ����������ѹ��
function Huffman_Callback(hObject, eventdata, handles)
Huffen_Compression;
% 2.����ѹ�� 
function lossy_Callback(hObject, eventdata, handles)
%  ͼ��DCTѹ�� 
function dct_Callback(hObject, eventdata, handles)
DCT_Compression(handles)
function pushbutton7_Callback(hObject, eventdata, handles)
DCT_Compression(handles);



%�������
function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit1_ButtonDownFcn(hObject, eventdata, handles)
function spatial_Callback(hObject, eventdata, handles)
function spectrum_Callback(hObject, eventdata, handles)
function Untitled_8_Callback(hObject, eventdata, handles)
function Untitled_9_Callback(hObject, eventdata, handles)
function Untitled_12_Callback(hObject, eventdata, handles)
function Untitled_13_Callback(hObject, eventdata, handles)
function Untitled_5_Callback(hObject, eventdata, handles)
function Untitled_6_Callback(hObject, eventdata, handles)
function Untitled_7_Callback(hObject, eventdata, handles)
function Untitled_10_Callback(hObject, eventdata, handles)
function Untitled_11_Callback(hObject, eventdata, handles)
function Grayscale_Callback(hObject, eventdata, handles)
function non_Linear_transform_Callback(hObject, eventdata, handles)
function subtraction_Callback(hObject, eventdata, handles)
function fusion_Callback(hObject, eventdata, handles)
function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Untitled_2_Callback(hObject, eventdata, handles)
function Untitled_3_Callback(hObject, eventdata, handles)
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

function edit7_Callback(hObject, eventdata, handles)
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit8_Callback(hObject, eventdata, handles)
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit9_Callback(hObject, eventdata, handles)
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit10_Callback(hObject, eventdata, handles)
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit11_Callback(hObject, eventdata, handles)
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit12_Callback(hObject, eventdata, handles)
function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
