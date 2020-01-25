%通过函数变化对图像亮度进行调整
[X,map] = imread('forest.tif');
I = ind2gray(X,map);
%把索引图像转换成灰度图像
imshow(I)
title('原图像')
figure;subplot(121)
plot(0:0.01:1,sqrt(0:0.01:1))
axis square 
title('平方根灰度变换函数')
subplot(122)
J=sqrt(double(I));
%平方根变换
imshow(uint8(J))
title('平方根变换后图像')
figure;subplot(121)
plot(0:0.01:1,1-(0:0.01:1))
axis square 
title('图像反转曲线')
subplot(122)
K=255-J;
%对平方根变换后的图像反转变换
imshow(uint8(K))
title('平方根变换后图像反转变换')
