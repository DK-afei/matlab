%��ֵ�˲�����
I = imread('eight.tif');
J = imnoise(I,'salt & pepper',0.02);
%ͼ������ν�����
K = medfilt2(J);
%ȱʡ3��3�����򴰵���ֵ�˲�
subplot(121),imshow(J)
title('����ν�����ͼ��')
subplot(122),imshow(K)
title('ȱʡ3��3���򴰵���ֵ�˲�ͼ��')
