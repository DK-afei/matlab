I = imread('lena.bmp');
imshow(I);
BW1 = edge(I,'roberts');
% ���Զ���ֵѡ�񷨶�ͼ�����Roberts���Ӽ��
[BW1,thresh1] = edge(I,'roberts');
% ���ص�ǰRoberts���ӱ�Ե������ֵ
disp('Roberts�����Զ�ѡ�����ֵΪ��')
disp(thresh1)
subplot(121),imshow(BW1);
title('�Զ���ֵ��Roberts���ӱ�Ե���')
BW1 = edge(I,'roberts',0.05);
% ����ֵΪ0.05��ͼ�����Roberts���Ӽ��
subplot(122),imshow(BW1);
title('��ֵΪ0.05��Roberts���ӱ�Ե���')
BW2 = edge(I,'sobel');
% ���Զ���ֵѡ�񷨶�ͼ�����Sobel���Ӽ��
figure,subplot(131),imshow(BW2);
title('�Զ���ֵ��Sobel���ӱ�Ե���')
[BW2,thresh2] = edge(I,'sobel');
% ���ص�ǰSobel���ӱ�Ե������ֵ
disp('Sobel�����Զ�ѡ�����ֵΪ��')
disp(thresh2)
BW2 = edge(I,'sobel',0.05,'horizontal');
% ����ֵΪ0.05ˮƽ�����ͼ�����Sobel���Ӽ��
subplot(132),imshow(BW2);
title('��ֵ0.05ˮƽ����Sobel����')
BW2 = edge(I,'sobel',0.05,'vertical');
% ����ֵΪ0.05��ֱ�����ͼ�����Sobel���Ӽ��
subplot(133),imshow(BW2);
title('��ֵ0.05��ֱ����Sobel����')
BW3 = edge(I,'prewitt');
% ���Զ���ֵѡ�񷨶�ͼ�����Prewitt���Ӽ��
figure,subplot(131),imshow(BW3);
title('�Զ���ֵ��Prewitt���ӱ�Ե���')
[BW3,thresh3] = edge(I,'prewitt');
% ���ص�ǰPrewitt���ӱ�Ե������ֵ
disp('Prewitt�����Զ�ѡ�����ֵΪ��')
disp(thresh3)
BW3 = edge(I,'prewitt',0.05,'horizontal');
% ����ֵΪ0.05ˮƽ�����ͼ�����Prewitt���Ӽ��
subplot(132),imshow(BW3);
title('��ֵ0.05ˮƽ����Prewitt����')
BW3 = edge(I,'prewitt',0.05,'vertical');
% ����ֵΪ0.05��ֱ�����ͼ�����Prewitt���Ӽ��
subplot(133),imshow(BW3);
title('��ֵ0.05��ֱ����Prewitt����')
