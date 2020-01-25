close all
clc
[fn,pn,fi]=uigetfile('III/粤A3Y347.jpg')
YuanShi=imread([pn fn]);%输入原始图像
figure(1);subplot(3,2,1),imshow(YuanShi),title('原始图像');
%1、图像预处理
YuanShiHuiDu=rgb2gray(YuanShi);%转化为灰度图像
subplot(3,2,2),imshow(YuanShiHuiDu),title('灰度图像');

BianYuan=edge(YuanShiHuiDu,'canny',0.5);%Canny算子边缘检测
subplot(3,2,3),imshow(BianYuan),title('Canny算子边缘检测后图像');

se1=[1;1;1]; %线型结构元素 
FuShi=imerode(BianYuan,se1);    %腐蚀图像
subplot(3,2,4),imshow(FuShi),title('腐蚀后边缘图像');

se2=strel('rectangle',[25,25]); %矩形结构元素
TianChong=imclose(FuShi,se2);%图像聚类、填充图像
subplot(3,2,5),imshow(TianChong),title('填充后图像');

YuanShiLvBo=bwareaopen(TianChong,2000);%从对象中移除面积小于2000的小对象
figure(2);
subplot(2,2,1),imshow(YuanShiLvBo),title('形态滤波后图像');
%2、车牌定位
[y,x]=size(YuanShiLvBo);%size函数将数组的行数返回到第一个输出变量，将数组的列数返回到第二个输出变量
YuCuDingWei=double(YuanShiLvBo);
%2.1、车牌粗定位之一确定行的起始位置和终止位置
%找出行累计和最大行，然后从累计和最大行开始索引
%自增自减来找到行起始/终止位置
Y1=zeros(y,1);%产生y行1列全零数组
for i=1:y
    for j=1:x
        if(YuCuDingWei(i,j)==1)
            Y1(i,1)= Y1(i,1)+1;%白色像素点统计
        end
    end
end
[temp,MaxY]=max(Y1);%Y方向车牌区域确定。返回行向量temp和MaxY，temp向量记录Y1的每列的最大值，MaxY向量记录Y1每列最大值的行号
subplot(2,2,2),plot(0:y-1,Y1),title('原图行方向像素点值累计和'),xlabel('行值'),ylabel('像素'); 
PY1=MaxY;
while ((Y1(PY1,1)>=50)&&(PY1>1))
        PY1=PY1-1;
end
PY2=MaxY;
while ((Y1(PY2,1)>=50)&&(PY2<y))
        PY2=PY2+1;
end
IY=YuanShi(PY1:PY2,:,:);
%2.2、车牌粗定位之二确定列的起始位置和终止位置
%找出列累计和最大列，然后从累计和最大列开始索引
%自增自减来找到列起始/终止位置
X1=zeros(1,x);%产生1行x列全零数组
for j=1:x
    for i=PY1:PY2
        if(YuCuDingWei(i,j,1)==1)
                X1(1,j)= X1(1,j)+1;               
         end  
    end       
end
subplot(2,2,4),plot(0:x-1,X1),title('原图列方向像素点值累计和'),xlabel('列值'),ylabel('像数');
PX1=1;
while ((X1(1,PX1)<3)&&(PX1<x))
       PX1=PX1+1;
end    
PX3=x;
while ((X1(1,PX3)<3)&&(PX3>PX1))
        PX3=PX3-1;
end
CuDingWei=YuanShi(PY1:PY2,PX1:PX3,:);
subplot(2,2,3),imshow(CuDingWei),title('粗定位后的彩色车牌图像')
%2.3、车牌精定位之一预处理
CuDingWei(:, :, 3) = CuDingWei(:, :, 1);
CuDingWei(:, :, 2) = CuDingWei(:, :, 1);
figure
imshow(CuDingWei)
CuDingWeiHuiDu=rgb2gray(CuDingWei); %将RGB图像转化为灰度图像

CuDingWeiErZhi = im2bw(CuDingWeiHuiDu,0);
[W, H] = size(CuDingWeiErZhi);
CuDingWeiErZhi = zeros([W + 20, H + 20]);
delta = 0;
for i = 1 : W 
    for j = 1 : H
        if CuDingWei(i, j, 1) + CuDingWei(i, j, 2) >= 180
            CuDingWeiErZhi(i + delta, j + delta) = 0;
        else
            CuDingWeiErZhi(i + delta, j + delta) = 1;
        end
    end
end
figure(3);
subplot(2,2,1),imshow(CuDingWeiErZhi),title('粗定位的二值车牌图像')%DingWei
%2.4、车牌精定位之二去除边框干扰
[r,s]=size(CuDingWeiErZhi);%size函数将数组的行数返回到第一个输出变量，将数组的列数返回到第二个输出变量
YuJingDingWei=double(CuDingWeiErZhi);%;CuDingWeiErZhi
X2=zeros(1,s);%产生1行s列全零数组
for i=1:r
    for j=1:s
        if(YuJingDingWei(i,j)==1)
            X2(1,j)= X2(1,j)+1;%白色像素点统计
        end
    end
end
[temp,MaxX]=max(X2);
subplot(2,2,2),plot(0:s-1,X2),title('粗定位车牌图像列方向像素点值累计和'),xlabel('列值'),ylabel('像素');
%2.4.1、去除左侧边框干扰
[g,h]=size(YuJingDingWei);
ZuoKuanDu=0;YouKuanDu=0;KuanDuYuZhi=5;
while sum(YuJingDingWei(:,ZuoKuanDu+1))~=0
    ZuoKuanDu=ZuoKuanDu+1;
end
if ZuoKuanDu<KuanDuYuZhi   % 认为是左侧干扰
    YuJingDingWei(:,[1:ZuoKuanDu])=0;%给图像d中1到KuanDu宽度间的点赋值为零
    YuJingDingWei=QieGe(YuJingDingWei); %值为零的点会被切割
end
subplot(2,2,3),imshow(YuJingDingWei),title('去除左侧边框的二值车牌图像')
%2.4.1、去除右侧边框干扰
[e,f]=size(YuJingDingWei);%上一步裁剪了一次，所以需要再次获取图像大小
d=f;
while sum(YuJingDingWei(:,d-1))~=0
    YouKuanDu=YouKuanDu+1;
    d=d-1;
end
if YouKuanDu<KuanDuYuZhi   % 认为是右侧干扰
    YuJingDingWei(:,[(f-YouKuanDu):f])=0;%
    YuJingDingWei=QieGe(YuJingDingWei); %值为零的点会被切割
end
subplot(2,2,4),imshow(YuJingDingWei),title('精确定位的车牌二值图像')


%%
ChePaiErZhi=YuJingDingWei;%logical()
ChePaiLvBo=bwareaopen(ChePaiErZhi,20);
subplot(1,2,1),imshow(ChePaiLvBo),title('形态学滤波后的车牌二值图像')
ChePaiYuFenGe=double(ChePaiLvBo);
figure
imshow(ChePaiLvBo)
[p,q] = size(ChePaiLvBo);
XX = zeros(1,q);
for j=1:q
    for i=1:p
       if(ChePaiLvBo(i,j)==0) 
           XX(1,j)=XX(1,j)+1;
       end
    end
end
i = 1;
while ((XX(1,i)==0)) 
    i = i + 1;
end
l = i;
count = 0;
figure
while (i<q)
    while ((XX(1,i)~=0)&&(i<q))
        i = i + 1;
    end
    r = i;
    part = ChePaiLvBo(:,l:r);
    count = count + 1;
    subplot(1,8,count)
    imshow(part)
    while ((XX(1,i)==0)&&(i<q))
        i = i + 1;
    end
    l = i;
    hold on;
end
subplot(1,2,2),plot(0:q-1,XX),title('车牌列方向像素点灰度值累计和'),xlabel('列值'),ylabel('累计像素量');



