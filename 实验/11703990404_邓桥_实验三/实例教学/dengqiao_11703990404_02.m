%ͼ�����ȵĵ���
I = imread('lax.tif');
J = imadjust(I,[0 0.2],[0.5 1]);
%�任[0 0.2]--->[0.5 1]
imshow(I)
figure, imshow(J)
