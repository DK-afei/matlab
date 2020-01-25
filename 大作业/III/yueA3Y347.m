close all
clc
[fn,pn,fi]=uigetfile('III/��A3Y347.jpg')
YuanShi=imread([pn fn]);%����ԭʼͼ��
figure(1);subplot(3,2,1),imshow(YuanShi),title('ԭʼͼ��');
%1��ͼ��Ԥ����
YuanShiHuiDu=rgb2gray(YuanShi);%ת��Ϊ�Ҷ�ͼ��
subplot(3,2,2),imshow(YuanShiHuiDu),title('�Ҷ�ͼ��');

BianYuan=edge(YuanShiHuiDu,'canny',0.5);%Canny���ӱ�Ե���
subplot(3,2,3),imshow(BianYuan),title('Canny���ӱ�Ե����ͼ��');

se1=[1;1;1]; %���ͽṹԪ�� 
FuShi=imerode(BianYuan,se1);    %��ʴͼ��
subplot(3,2,4),imshow(FuShi),title('��ʴ���Եͼ��');

se2=strel('rectangle',[25,25]); %���νṹԪ��
TianChong=imclose(FuShi,se2);%ͼ����ࡢ���ͼ��
subplot(3,2,5),imshow(TianChong),title('����ͼ��');

YuanShiLvBo=bwareaopen(TianChong,2000);%�Ӷ������Ƴ����С��2000��С����
figure(2);
subplot(2,2,1),imshow(YuanShiLvBo),title('��̬�˲���ͼ��');
%2�����ƶ�λ
[y,x]=size(YuanShiLvBo);%size������������������ص���һ�������������������������ص��ڶ����������
YuCuDingWei=double(YuanShiLvBo);
%2.1�����ƴֶ�λ֮һȷ���е���ʼλ�ú���ֹλ��
%�ҳ����ۼƺ�����У�Ȼ����ۼƺ�����п�ʼ����
%�����Լ����ҵ�����ʼ/��ֹλ��
Y1=zeros(y,1);%����y��1��ȫ������
for i=1:y
    for j=1:x
        if(YuCuDingWei(i,j)==1)
            Y1(i,1)= Y1(i,1)+1;%��ɫ���ص�ͳ��
        end
    end
end
[temp,MaxY]=max(Y1);%Y����������ȷ��������������temp��MaxY��temp������¼Y1��ÿ�е����ֵ��MaxY������¼Y1ÿ�����ֵ���к�
subplot(2,2,2),plot(0:y-1,Y1),title('ԭͼ�з������ص�ֵ�ۼƺ�'),xlabel('��ֵ'),ylabel('����'); 
PY1=MaxY;
while ((Y1(PY1,1)>=50)&&(PY1>1))
        PY1=PY1-1;
end
PY2=MaxY;
while ((Y1(PY2,1)>=50)&&(PY2<y))
        PY2=PY2+1;
end
IY=YuanShi(PY1:PY2,:,:);
%2.2�����ƴֶ�λ֮��ȷ���е���ʼλ�ú���ֹλ��
%�ҳ����ۼƺ�����У�Ȼ����ۼƺ�����п�ʼ����
%�����Լ����ҵ�����ʼ/��ֹλ��
X1=zeros(1,x);%����1��x��ȫ������
for j=1:x
    for i=PY1:PY2
        if(YuCuDingWei(i,j,1)==1)
                X1(1,j)= X1(1,j)+1;               
         end  
    end       
end
subplot(2,2,4),plot(0:x-1,X1),title('ԭͼ�з������ص�ֵ�ۼƺ�'),xlabel('��ֵ'),ylabel('����');
PX1=1;
while ((X1(1,PX1)<3)&&(PX1<x))
       PX1=PX1+1;
end    
PX3=x;
while ((X1(1,PX3)<3)&&(PX3>PX1))
        PX3=PX3-1;
end
CuDingWei=YuanShi(PY1:PY2,PX1:PX3,:);
subplot(2,2,3),imshow(CuDingWei),title('�ֶ�λ��Ĳ�ɫ����ͼ��')
%2.3�����ƾ���λ֮һԤ����
CuDingWei(:, :, 3) = CuDingWei(:, :, 1);
CuDingWei(:, :, 2) = CuDingWei(:, :, 1);
figure
imshow(CuDingWei)
CuDingWeiHuiDu=rgb2gray(CuDingWei); %��RGBͼ��ת��Ϊ�Ҷ�ͼ��

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
subplot(2,2,1),imshow(CuDingWeiErZhi),title('�ֶ�λ�Ķ�ֵ����ͼ��')%DingWei
%2.4�����ƾ���λ֮��ȥ���߿����
[r,s]=size(CuDingWeiErZhi);%size������������������ص���һ�������������������������ص��ڶ����������
YuJingDingWei=double(CuDingWeiErZhi);%;CuDingWeiErZhi
X2=zeros(1,s);%����1��s��ȫ������
for i=1:r
    for j=1:s
        if(YuJingDingWei(i,j)==1)
            X2(1,j)= X2(1,j)+1;%��ɫ���ص�ͳ��
        end
    end
end
[temp,MaxX]=max(X2);
subplot(2,2,2),plot(0:s-1,X2),title('�ֶ�λ����ͼ���з������ص�ֵ�ۼƺ�'),xlabel('��ֵ'),ylabel('����');
%2.4.1��ȥ�����߿����
[g,h]=size(YuJingDingWei);
ZuoKuanDu=0;YouKuanDu=0;KuanDuYuZhi=5;
while sum(YuJingDingWei(:,ZuoKuanDu+1))~=0
    ZuoKuanDu=ZuoKuanDu+1;
end
if ZuoKuanDu<KuanDuYuZhi   % ��Ϊ��������
    YuJingDingWei(:,[1:ZuoKuanDu])=0;%��ͼ��d��1��KuanDu��ȼ�ĵ㸳ֵΪ��
    YuJingDingWei=QieGe(YuJingDingWei); %ֵΪ��ĵ�ᱻ�и�
end
subplot(2,2,3),imshow(YuJingDingWei),title('ȥ�����߿�Ķ�ֵ����ͼ��')
%2.4.1��ȥ���Ҳ�߿����
[e,f]=size(YuJingDingWei);%��һ���ü���һ�Σ�������Ҫ�ٴλ�ȡͼ���С
d=f;
while sum(YuJingDingWei(:,d-1))~=0
    YouKuanDu=YouKuanDu+1;
    d=d-1;
end
if YouKuanDu<KuanDuYuZhi   % ��Ϊ���Ҳ����
    YuJingDingWei(:,[(f-YouKuanDu):f])=0;%
    YuJingDingWei=QieGe(YuJingDingWei); %ֵΪ��ĵ�ᱻ�и�
end
subplot(2,2,4),imshow(YuJingDingWei),title('��ȷ��λ�ĳ��ƶ�ֵͼ��')


%%
ChePaiErZhi=YuJingDingWei;%logical()
ChePaiLvBo=bwareaopen(ChePaiErZhi,20);
subplot(1,2,1),imshow(ChePaiLvBo),title('��̬ѧ�˲���ĳ��ƶ�ֵͼ��')
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
subplot(1,2,2),plot(0:q-1,XX),title('�����з������ص�Ҷ�ֵ�ۼƺ�'),xlabel('��ֵ'),ylabel('�ۼ�������');



