%小波变换压缩图像
clear
load wbarb;
%装入图像
subplot(221);
X=imread('lena.bmp')
[m n]=size(X);
image(X);
%显示图像
colormap(map)
title('原始图像');
axis square
disp('压缩前图像X的大小');
whos('X')
[c,s]=wavedec2(X,2,'bior3.7');
%对图像用小波进行层分解
cal=appcoef2(c,s,'bior3.7',1);
%提取小波分解结构中的一层的低频系数和高频系数
ch1=detcoef2('h',c,s,1);
%水平方向
cv1=detcoef2('v',c,s,1);
%垂直方向
cd1=detcoef2('d',c,s,1);
%斜线方向
a1=wrcoef2('a',c,s,'bior3.7',1);
h1=wrcoef2('h',c,s,'bior3.7',1);
v1=wrcoef2('v',c,s,'bior3.7',1);
d1=wrcoef2('d',c,s,'bior3.7',1);
%各频率成份重构
c1=[a1,h1;v1,d1];
subplot(222);
image(c1);
%显示分频信息
axis square;
title ('分解后低频和高频信息');
%进行图像压缩
%保留小波分解第一层低频信息
%首先对第一层信息进行量化编码
ca1=appcoef2(c,s,'bior3.7',1);
ca1=wcodemat(ca1,440,'mat',0);
ca1=0.5*ca1;
subplot(223);
image(ca1);
%改变图像高度并显示
colormap(map);
axis square;
title('第一次压缩图像');
disp('第一次压缩图像的大小为：');
whos('ca1')
disp('图像的压缩比：')
disp(m*n/length(ca1))
ca2=appcoef2(c,s,'bior3.7',2);
%保留小波分解第二层低频信息进行压缩
ca2=wcodemat(ca2,440,'mat',0);
%首先对第二层信息进行量化编码
ca2=0.25*ca2;
%改变图像高度并显示
subplot(224);
image(ca2);
colormap(map);
axis square;
title('第二次压缩图像');
disp('第二次压缩图像的大小为：');
whos('ca2')
disp('图像的压缩比：')
disp(m*n/length(ca2))