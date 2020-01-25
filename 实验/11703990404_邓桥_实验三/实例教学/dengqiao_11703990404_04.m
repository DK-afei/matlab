%中值滤波处理
I = imread('eight.tif');
J = imnoise(I,'salt & pepper',0.02);
%图像添加盐椒噪声
K = medfilt2(J);
%缺省3×3的邻域窗的中值滤波
subplot(121),imshow(J)
title('添加盐椒噪声图像')
subplot(122),imshow(K)
title('缺省3×3邻域窗的中值滤波图像')
