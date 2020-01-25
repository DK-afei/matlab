t=imread('lena.bmp');
[m,n,z]=size(t);
y=0+0.1*randn(m,n);%二维高斯分布矩阵 0是均值 0.1是标准差

%先将其double化，再除以255 便于后面计算
t1=double(t)/255;

%加上噪声
t1=t1+y;

%将像素范围扩大至0--255
t1=t1*255;

%转换为uint8类型
t1=uint8(t1);

subplot(3,2,1),imshow(t),title('原图');
subplot(3,2,3),imshow(t1),title('加入均值为0，标准差为0.1的高斯噪声后');

K = medfilt2(t1);
subplot(3,2,4),imshow(K)
title('中值滤波图像')

K1= filter2(fspecial('average',3),t1)/255;
%应用3×3邻域窗口法
subplot(3,2,5),imshow(K1)
title('3×3窗的邻域平均滤波图像')
K2= filter2(fspecial('average',7),t1)/255;
%应用7×7邻域窗口法
subplot(3,2,6),imshow(K2)
title('7×7窗的邻域平均滤波图像')
