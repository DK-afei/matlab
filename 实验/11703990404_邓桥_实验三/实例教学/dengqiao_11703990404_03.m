%��ֵ�˲�����
I = imread('eight.tif');
J = imnoise(I,'salt & pepper',0.02);
%��ӽ�������
subplot(221),imshow(I)
title('ԭͼ��')
subplot(222),imshow(J)
title('��ӽ�������ͼ��')
K1= filter2(fspecial('average',3),J)/255;
%Ӧ��3��3���򴰿ڷ�
subplot(223),imshow(K1)
title('3��3��������ƽ���˲�ͼ��')
K2= filter2(fspecial('average',7),J)/255;
%Ӧ��7��7���򴰿ڷ�
subplot(224),imshow(K2)
title('7��7��������ƽ���˲�ͼ��')
