%Canny���Ӽ��
I = imread('lena.bmp');
BW = edge(I,'canny');
% ���Զ���ֵѡ�񷨶�ͼ�����Canny���Ӽ��
[BW,thresh] = edge(I,'canny');
% ���ص�ǰCanny���ӱ�Ե������ֵ
disp('Canny�����Զ�ѡ�����ֵΪ��')
disp(thresh)
subplot(121),imshow(BW);
title('�Զ���ֵ��Canny���ӱ�Ե���')
BW = edge(I,'Canny',[0.2 0.5]);
% ����ֵΪ[0.1 0.5]��ͼ�����Canny���Ӽ��
subplot(122),imshow(BW);
title('��ֵΪ[0.1 0.5]��Canny���ӱ�Ե���')
