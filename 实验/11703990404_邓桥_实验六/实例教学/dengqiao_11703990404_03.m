%Canny算子检测
I = imread('lena.bmp');
BW = edge(I,'canny');
% 以自动阈值选择法对图像进行Canny算子检测
[BW,thresh] = edge(I,'canny');
% 返回当前Canny算子边缘检测的阈值
disp('Canny算子自动选择的阈值为：')
disp(thresh)
subplot(121),imshow(BW);
title('自动阈值的Canny算子边缘检测')
BW = edge(I,'Canny',[0.2 0.5]);
% 以阈值为[0.1 0.5]对图像进行Canny算子检测
subplot(122),imshow(BW);
title('阈值为[0.1 0.5]的Canny算子边缘检测')
