t=imread('lena.bmp');
[m,n,z]=size(t);
y=0+0.1*randn(m,n);%��ά��˹�ֲ����� 0�Ǿ�ֵ 0.1�Ǳ�׼��

%�Ƚ���double�����ٳ���255 ���ں������
t1=double(t)/255;

%��������
t1=t1+y;

%�����ط�Χ������0--255
t1=t1*255;

%ת��Ϊuint8����
t1=uint8(t1);

subplot(3,2,1),imshow(t),title('ԭͼ');
subplot(3,2,3),imshow(t1),title('�����ֵΪ0����׼��Ϊ0.1�ĸ�˹������');

K = medfilt2(t1);
subplot(3,2,4),imshow(K)
title('��ֵ�˲�ͼ��')

K1= filter2(fspecial('average',3),t1)/255;
%Ӧ��3��3���򴰿ڷ�
subplot(3,2,5),imshow(K1)
title('3��3��������ƽ���˲�ͼ��')
K2= filter2(fspecial('average',7),t1)/255;
%Ӧ��7��7���򴰿ڷ�
subplot(3,2,6),imshow(K2)
title('7��7��������ƽ���˲�ͼ��')
