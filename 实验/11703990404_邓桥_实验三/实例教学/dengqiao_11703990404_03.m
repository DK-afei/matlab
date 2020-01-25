%均值滤波处理
I = imread('eight.tif');
J = imnoise(I,'salt & pepper',0.02);
%添加椒盐噪声
subplot(221),imshow(I)
title('原图像')
subplot(222),imshow(J)
title('添加椒盐噪声图像')
K1= filter2(fspecial('average',3),J)/255;
%应用3×3邻域窗口法
subplot(223),imshow(K1)
title('3×3窗的邻域平均滤波图像')
K2= filter2(fspecial('average',7),J)/255;
%应用7×7邻域窗口法
subplot(224),imshow(K2)
title('7×7窗的邻域平均滤波图像')
