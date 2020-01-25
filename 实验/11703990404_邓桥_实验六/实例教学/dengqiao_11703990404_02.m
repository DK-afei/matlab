%Log算子与零交叉检测
I = imread('lena.bmp');
BW1 = edge(I,'log');
% 以自动阈值选择法对图像进行Log算子检测
[BW1,thresh1] = edge(I,'log');
% 返回当前Log算子边缘检测的阈值
disp('Log算子自动选择的阈值为：')
disp(thresh1)
subplot(121),imshow(BW1);
title('自动阈值的Log算子边缘检测')
BW1 = edge(I,'log',0.005);
% 以阈值为0.005对图像进行Log算子检测
subplot(122),imshow(BW1);
title('阈值为0.005的Log算子边缘检测')
h=fspecial('gaussian',5);
% 设计高斯滤波器
[BW2,thresh2] = edge(I,'zerocross',[],h);
% 返回当前零交叉检测边缘检测的阈值
disp('零交叉检测自动选择的阈值为：')
disp(thresh2)
figure,subplot(121),imshow(BW2);
title('自动阈值的零交叉边缘检测')
BW2 = edge(I,'zerocross',0.03,h);
% 以阈值为0.03对图像进行零交叉检测
subplot(122),imshow(BW2);
title('阈值为0.03的零交叉边缘检测')
