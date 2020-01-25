I = imread('lena.bmp');
imshow(I);
BW1 = edge(I,'roberts');
% 以自动阈值选择法对图像进行Roberts算子检测
[BW1,thresh1] = edge(I,'roberts');
% 返回当前Roberts算子边缘检测的阈值
disp('Roberts算子自动选择的阈值为：')
disp(thresh1)
subplot(121),imshow(BW1);
title('自动阈值的Roberts算子边缘检测')
BW1 = edge(I,'roberts',0.05);
% 以阈值为0.05对图像进行Roberts算子检测
subplot(122),imshow(BW1);
title('阈值为0.05的Roberts算子边缘检测')
BW2 = edge(I,'sobel');
% 以自动阈值选择法对图像进行Sobel算子检测
figure,subplot(131),imshow(BW2);
title('自动阈值的Sobel算子边缘检测')
[BW2,thresh2] = edge(I,'sobel');
% 返回当前Sobel算子边缘检测的阈值
disp('Sobel算子自动选择的阈值为：')
disp(thresh2)
BW2 = edge(I,'sobel',0.05,'horizontal');
% 以阈值为0.05水平方向对图像进行Sobel算子检测
subplot(132),imshow(BW2);
title('阈值0.05水平方向Sobel算子')
BW2 = edge(I,'sobel',0.05,'vertical');
% 以阈值为0.05垂直方向对图像进行Sobel算子检测
subplot(133),imshow(BW2);
title('阈值0.05垂直方向Sobel算子')
BW3 = edge(I,'prewitt');
% 以自动阈值选择法对图像进行Prewitt算子检测
figure,subplot(131),imshow(BW3);
title('自动阈值的Prewitt算子边缘检测')
[BW3,thresh3] = edge(I,'prewitt');
% 返回当前Prewitt算子边缘检测的阈值
disp('Prewitt算子自动选择的阈值为：')
disp(thresh3)
BW3 = edge(I,'prewitt',0.05,'horizontal');
% 以阈值为0.05水平方向对图像进行Prewitt算子检测
subplot(132),imshow(BW3);
title('阈值0.05水平方向Prewitt算子')
BW3 = edge(I,'prewitt',0.05,'vertical');
% 以阈值为0.05垂直方向对图像进行Prewitt算子检测
subplot(133),imshow(BW3);
title('阈值0.05垂直方向Prewitt算子')
