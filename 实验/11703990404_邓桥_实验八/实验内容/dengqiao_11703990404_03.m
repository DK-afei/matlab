%��ɢ���ұ任����ѹ��ͼ��
I = imread('lena.bmp');
[m n]=size(I);
I = im2double(I);
%ͼ��洢����ת��
T = dctmtx(8);
%��ɢ���ұ任����
B = blkproc(I,[8 8],'P1*x*P2',T,T');
mask = [1   1   1   1   0   0   0   0
       1   1   1   0   0   0   0   0
       1   1   0   0   0   0   0   0
       1   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0];
B2 = blkproc(B,[8 8],'P1.*x',mask);
I2 = blkproc(B2,[8 8],'P1*x*P2',T',T);
imshow(I)
title('ԭʼͼ��')
figure;
imshow(I2)
title('ѹ�����ͼ��')
disp('ͼ���ѹ���ȣ�')
disp(m*n/length(I2))

