%ͨ�������仯��ͼ�����Ƚ��е���
[X,map] = imread('forest.tif');
I = ind2gray(X,map);
%������ͼ��ת���ɻҶ�ͼ��
imshow(I)
title('ԭͼ��')
figure;subplot(121)
plot(0:0.01:1,sqrt(0:0.01:1))
axis square 
title('ƽ�����Ҷȱ任����')
subplot(122)
J=sqrt(double(I));
%ƽ�����任
imshow(uint8(J))
title('ƽ�����任��ͼ��')
figure;subplot(121)
plot(0:0.01:1,1-(0:0.01:1))
axis square 
title('ͼ��ת����')
subplot(122)
K=255-J;
%��ƽ�����任���ͼ��ת�任
imshow(uint8(K))
title('ƽ�����任��ͼ��ת�任')
