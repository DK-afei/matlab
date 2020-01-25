%数字图像处理主函数--------------------------------------------------------------------
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

%一、图像打开与保存 ----------------------------------------------------------
% 1.图像打开*****************************
function file_open_Callback(hObject, eventdata, handles)
axes(handles.axes1);
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'选择图像');
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
% 2.图像保存*****************************
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

%二、几何变换 ------------------------------------------------------------
function Geometric_transformation_Callback(hObject, eventdata, handles)
% 1.缩放*****************************
function enlarge_narrow_Callback(hObject, eventdata, handles)
enlarge_narrow;
% 2.平移*****************************
function translation_Callback(hObject, eventdata, handles)
ti;
% 3.旋转*****************************
function rotate_Callback(hObject, eventdata, handles)
rotate;
% 4.翻转*****************************
function flip_Callback(hObject, eventdata, handles)
mflip;

%三、图像变换 --------------------------------------------------------------
% 1.图像合成***********************
function synthesis_Callback(hObject, eventdata, handles)
%       a图像融合
function blend_Callback(hObject, eventdata, handles)
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'选择图像');
str=[pathname filename];
I1 = imread(str);  
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'选择图像');
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
    sprintf('请重新选择图片');
end
%       b直接拼接
function joint_Callback(hObject, eventdata, handles)
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'选择图像');
str=[pathname filename];
I1 = imread(str);  
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.tif';'*.jpg'},'选择图像');
str=[pathname filename];
I2 = imread(str); 
 img1=I1;
 img2=I2;
[H,W,k]=size(img1);
l_r=5;%重叠宽度（W-宽 至 W）---如果不用特征匹配这里直接写重合区宽
L=W+1-l_r;%左边起点
R=W;%右边尾点
n=R-L+1;%重叠宽度：就是l_r
%直接拼接图
im=[img1,img2(:,n:W,:)];%1全图+2的后面部分
figure;imshow(im);title('直接拼接图');
% 2.灰度变换***********************
%       a图像反转
function reversion_Callback(hObject, eventdata, handles)
global im;
J=double(im);
J=-J+(256-1);                 %线性变换
H=uint8(J);
axes(handles.axes2);
imshow(H);	
%       b灰度变换
function Linear_transform_Callback(hObject, eventdata, handles)
gray_trans;
% 3.直方图均衡化*******************
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
nf = m*n;%总像素个数
h=zeros(1,256);%每个灰度级像素个数
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
    hs(l)=h(l)/nf;%每个灰度级像素个数占百分比
end
%直方图均衡化
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
%灰度值映射用g(i)替换原图 所有灰度值为i的像素点
for i = 1:m
    for j = 1: n
        f(i,j) = g(f(i,j));
    end
end
axes(handles.axes2);
imshow(f);%显示均衡化后的图
title('均衡化后的图');

%四、图像去噪 --------------------------------------------------------------
% 1.空域滤波***********************
%      a中值滤波
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
exsize=floor(masksize/2);   %各方向扩展大小
Imgex=padarray(Img,[exsize,exsize],'replicate','both'); %扩展图片
[m,n]=size(Img);
Img_out=Img;    %将Img_out准备为和Img相同的size
for i=1:m
    for j=1:n
        neighbor=Imgex(i:i+masksize-1,j:j+masksize-1);  %截取邻域
        Img_out(i,j)=median(neighbor(:));   %中值滤波
    end
end
axes(handles.axes2);
imshow(Img_out);
%          中值滤波模板设置
function pushbutton2_Callback(hObject, eventdata, handles)
medfilt_Callback(hObject, eventdata, handles);
%      b均值滤波
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
%x是需要滤波的图像,n是模板大小(即n×n)  
n=data1_num;
a(1:n,1:n)=1;   %a即n×n模板,元素全是1  
[height, width]=size(x);   %输入图像是hightxwidth的,且hight>n,width>n  
x1=double(x);  
x2=x1;  
for i=1:height-n+1  
	    for j=1:width-n+1  
	        c=x1(i:i+(n-1),j:j+(n-1)).*a; %取出x1中从(i,j)开始的n行n列元素与模板相乘  
        s=sum(sum(c));                 %求c矩阵中各元素之和  
        x2(i+(n-1)/2,j+(n-1)/2)=s/(n*n); %将与模板运算后的各元素的均值赋给模板中心位置的元素  
    end  
end  
%未被赋值的元素取原值  
imshow(uint8(x2));
%          均值滤波模板设置
function pushbutton3_Callback(hObject, eventdata, handles)
Mean_Callback(hObject, eventdata, handles);
% 2.频域滤波***********************
%      a高通滤波
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
%% 高通滤波，锐化图像
moon=I1;
[m,n]=size(moon);
I=im2double(moon);
I_spectrum=fft2(I);
I_spectrum=fftshift(I_spectrum);
%% highpass filter
H=zeros(m,n);
centerx=m/2;
centery=n/2;
D0=data1_num;  % 可调节通带半径来控制通过的高频分量
for x=1:m
    for y=1:n
        H(x,y)=exp(-((x-centerx)^2+(y-centery)^2)/(2*D0^2));  %计算高斯滤波模板
    end
end
H=1-H;

g1=H.*I_spectrum;% 高频部分图像，边缘图像
g2=g1+I_spectrum;
% 先反中心化，在转到空域图像
g3=ifftshift(g2);
I2=real(ifft2(g3));
axes(handles.axes2);
imshow(real(ifft2(ifftshift(g1))),[]);
%          高通滤波阈值设置
function pushbutton6_Callback(hObject, eventdata, handles)
high_Callback(hObject, eventdata, handles);
%      b低通滤波
function low_Callback(hObject, eventdata, handles)
%三阶Butterworth低通滤波
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
%          低通滤波阈值设置 
function pushbutton5_Callback(hObject, eventdata, handles)
low_Callback(hObject, eventdata, handles);


% 五、图像边缘检测------------------------------------------------------------
function edge_detection_Callback(hObject, eventdata, handles)
% 1.Laplacian算子模板
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
% 2.Prewitt算子模板
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
% 3.sobel算子模板
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
% 4.Roberts算子模板
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
% 5.对数检测算子模板
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

%六、图像锐化----------------------------------------------------------------
function sharpen_Callback(hObject, eventdata, handles)
% 1.Laplacian算子锐化
function laplace_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
%I1为灰度图
I=im2double(I1);
KernelType=4;
c=1;
%扩展区域的行列数
KernelSize=3;
len=floor(KernelSize/2);
%对原始图像进行扩展，此处采用了镜像扩展，目的是解决边缘计算的问题
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
    %接下来两个是合成拉普拉斯算子
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
        %从扩展图像中，取出局部图像
        Block=f_pad(i-len:i+len,j-len:j+len);
        %将拉普拉斯算子的结果作用于原始图像，得到输出图像       
        g(i-len,j-len)=I(i-len,j-len)+ a*sum(sum(Block.*L));
        %保留拉普拉斯算子的运算结果
        edge(i-len,j-len)=a*sum(sum(Block.*L));
    end
end
figure,imshow(edge);
axes(handles.axes2);
imshow(g);
% 2.其他锐化
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

% 七、图像分割----------------------------------------------------------------
function image_segmentation_Callback(hObject, eventdata, handles)
% 1.阈值法********************
%        a手动阈值
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
%人工选定阈值进行分割，选择阈值为 data1_num
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
% figure;imshow(BW1),title('人工阈值进行分割');
%自动选择阈值
% T2=graythresh(I);
% BW2=im2bw(I,T2);%Otus阈值进行分割
% figure;imshow(BW2),title('Otus阈值进行分割');
axes(handles.axes2);
imshow(BW1);
function pushbutton4_Callback(hObject, eventdata, handles)
Threshold_Callback(hObject, eventdata, handles);
%        b自动阈值ostu
function ostu_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
I=I1;
%自动选择阈值
T2=graythresh(I);
BW2=im2bw(I,T2);%Otus阈值进行分割
axes(handles.axes2);
imshow(BW2);
title('Otus阈值进行分割');
% 2.区域生长
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
 [y,x]=getpts;             %获得区域生长起始点
  if x<0
      x=-x;
  end
  
  if y<0
      y=-y;
  end

x1=round(x);            %横坐标取整
y1=round(y);            %纵坐标取整
try
seed=I(x1,y1);           %将生长起始点灰度值存入seed中
catch
end
J=zeros(M,N);          %作一个全零与原图像等大的图像矩阵J，作为输出图像矩阵
J(x1,y1)=1;             %将J中与所取点相对应位置的点设置为白
sum=seed;              %储存符合区域生长条件的点的灰度值的和
suit=1;                 %储存符合区域生长条件的点的个数
count=1;               %记录每次判断一点周围八点符合条件的新点的数目
threshold=data3_num;         %阈值，注意需要和double类型存储的图像相符合
while count>0
    s=0;                   %记录判断一点周围八点时，符合条件的新点的灰度值之和
     count=0;
     for i=1:M
       for j=1:N
         if J(i,j)==1
          if (i-1)>0 && (i+1)<(M+1) && (j-1)>0 && (j+1)<(N+1)  %判断此点是否为图像边界上的点
           for u= -1:1                               %判断点周围八点是否符合阈值条件
            for v= -1:1
              if  J(i+u,j+v)==0 && abs(I(i+u,j+v)-seed)<=threshold&& 1/(1+1/15*abs(I(i+u,j+v)-seed))>0.8
                           J(i+u,j+v)=1;
                    %判断是否尚未标记，并且为符合阈值条件的点
                    %符合以上两条件即将其在J中与之位置对应的点设置为白
                 count=count+1;
                 s=s+I(i+u,j+v);                      %此点的灰度之加入s中
              end
            end
           end
          end
         end
       end
     end
    suit=suit+count;                                   %将n加入符合点数计数器中
    sum=sum+s;                                     %将s加入符合点的灰度值总合中
    seed=sum/suit;                                    %计算新的灰度平均值
end
axes(handles.axes2);
imshow(J);
%        区域生长法设置阈值
function pushbutton8_Callback(hObject, eventdata, handles)
Regional_growth_Callback(hObject, eventdata, handles);
% 3.区域分裂与合并
function Regional_division_and_merger_Callback(hObject, eventdata, handles)
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
% 使用区域分离和合并的图像分割
f = I1;
g=splitmerge(f,2,@predicate);%2代表分割中允许最小的块，predicate函数返回1，说明需要再分裂，返回0说明不需要继续分裂
% figure,imshow(g);
% title('mindim为2时的分割图像');
se=ones(8,8);
gdilate=imdilate(g,se);%膨胀是为了填充空洞
% figure;imshow(gdilate);
% title('膨胀后的图')
gerode=imerode(gdilate,se);%腐蚀是为了缩回原来大小
% figure;imshow(gerode);
% title('腐蚀后的图')
axes(handles.axes2);
imshow(gerode);
%%%%%%%%%%
function g=splitmerge(f,mindim,fun)%f是待分割的原图，mindim是定义分解中所允许的最小的块，必须是2的正整数次幂
Q=2^nextpow2(max(size(f)));
[M,N]=size(f);
f=padarray(f,[Q-M,Q-N],'post');%：填充图像或填充数组。f是输入图像，输出是填充后的图像，先将图像填充到2的幂次以使后面的分解可行
%然后是填充的行数和列数，post：表示在每一维的最后一个元素后填充,B = padarray(A,padsize,padval,direction)
%不含padval就用0填充,Q代表填充后图像的大小。
S=qtdecomp(f,@split_test,mindim,fun);%S传给split_test，qtdecomp divides a square image into four
% different sizes.S是包含四叉树结构的稀疏矩阵，存储的值是块的大小及坐标，以稀疏矩阵形式存储
Lmax=full(max(S(:)));%将以稀疏矩阵存储形式存储的矩阵变换成以普通矩阵（full matrix）形式存储，full，sparse只是存储形式的不同
g=zeros(size(f));
MARKER=zeros(size(f));
for k=1:Lmax
    [vals,r,c]=qtgetblk(f,S,k);%vals是一个数组，包含f的四叉树分解中大小为k*k的块的值，是一个k*k*个数的矩阵，
%个数是指S中有多少个这样大小的块，f是被四叉树分的原图像，r，c是对应的左上角块的坐标如2*2块，代表的是左上角开始块的坐标
        if ~isempty(vals)
            for I=1:length(r)
                    xlow=r(I);
                    ylow=c(I);
                    xhigh=xlow+k-1;
                    yhigh=ylow+k-1;
                    region=f(xlow:xhigh,ylow:yhigh);%找到对应的区域
                    flag=feval(fun,region);%evaluates the function handle, fhandle,using arguments x1 through xn.执行函数fun，region是参数
                    if flag%如果返回的是1，则进行标记
                            g(xlow:xhigh,ylow:yhigh)=1;%然后将对应的区域置1
                            MARKER(xlow,ylow)=1;%MARKER矩阵对应的左上角坐标置1
                    end
            end
        end
end	
	g=bwlabel(imreconstruct(MARKER,g));%imreconstruct默认2D图像8连通，这个函数就是起合的作用
	g=g(1:M,1:N);%返回原图像的大小
%%%%%%%%%%
function v=split_test(B,mindim,fun)
 K=size(B,3);%B就是qtdecomp函数传过来的，代表当前size(B,3)返回的是B的层数，就是B是几维的，这里实际上就是有几个B这样大小的图像块
%这句代码的意思是从qtdecomp函数传过来的B，是当前分解成的K块的m*m的图像块，K表示有多少个这样大小的图像块
 v(1:K)=false;
   for I=1:K
        quadregion=B(:,:,I);
        if size(quadregion,1)<=mindim%如果分的块的大小小于mindim就直接结束
                v(I)=false;
                continue
        end
        flag=feval(fun,quadregion);%quadregion是fun函数的参数
        if flag%如果flag是1，代表需要再分
            v(I)=true;%这里就相当于split_test是起一个调用predicate的作用，返回的就是ppredicate的值
        end
   end
function flag=predicate(region)
sd=std2(region);
m=mean2(region);
flag=(sd>20)&(m>26)&(m<255);
%predicate用于两个目的，在split_test中被调用时，判断是否该被分，如果满足这个条件就返回1，需要再分，否则就返回0，不能再被分
%在开始合并时调用它，用于判断是否该进行合并标记，返回1，表示通过测试，四叉区域都用1填充，返回0代表没有通过测试，用0填充
%所以区域分裂与合并的分割方法中，predicate函数是用户自定义的，改动性最大，最关键，关乎分割效果的好坏，用于两个目的时可以分别定义不同的
%函数准则，以达到最好的效果。

% 八、图像压缩----------------------------------------------------------------
function compression_Callback(hObject, eventdata, handles)
% 1.无损压缩 
function lossless_Callback(hObject, eventdata, handles)
%  霍夫曼编码压缩
function Huffman_Callback(hObject, eventdata, handles)
Huffen_Compression;
% 2.有损压缩 
function lossy_Callback(hObject, eventdata, handles)
%  图像DCT压缩 
function dct_Callback(hObject, eventdata, handles)
DCT_Compression(handles)
function pushbutton7_Callback(hObject, eventdata, handles)
DCT_Compression(handles);



%各个组件
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
