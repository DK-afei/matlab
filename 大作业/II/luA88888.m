close all
clc
CuDingWei=imread('II/³A88888.jpg')
imshow(CuDingWei)
%%%%%%%%%%����Ԥ����%%%%%%%%%%%
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
% return ;
subplot(2,2,1),imshow(CuDingWeiErZhi),title('��ֵ����ͼ��')%DingWei
%%%%%%%%%%ȥ���߿����%%%%%%%%%%%
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

%%%%%%%%%ȥ�����߿����%%%%%%%%%%%
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
%%%%%%%%ȥ���Ҳ�߿����%%%%%%%%%%%
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
subplot(2,2,4),imshow(YuJingDingWei),title('���ƶ�ֵͼ��')


%%
ChePaiErZhi=YuJingDingWei;
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



