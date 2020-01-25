%图像亮度的调整
I = imread('lax.tif');
J = imadjust(I,[0 0.2],[0.5 1]);
%变换[0 0.2]--->[0.5 1]
imshow(I)
figure, imshow(J)
