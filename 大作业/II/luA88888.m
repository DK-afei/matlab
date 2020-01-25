close all
clc
CuDingWei=imread('II/鲁A88888.jpg')
imshow(CuDingWei)
%%%%%%%%%%车牌预处理%%%%%%%%%%%
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
% return ;
subplot(2,2,1),imshow(CuDingWeiErZhi),title('二值车牌图像')%DingWei
%%%%%%%%%%去除边框干扰%%%%%%%%%%%
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

%%%%%%%%%去除左侧边框干扰%%%%%%%%%%%
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
%%%%%%%%去除右侧边框干扰%%%%%%%%%%%
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
subplot(2,2,4),imshow(YuJingDingWei),title('车牌二值图像')


%%
ChePaiErZhi=YuJingDingWei;
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



