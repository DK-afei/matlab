%Log�������㽻����
I = imread('lena.bmp');
BW1 = edge(I,'log');
% ���Զ���ֵѡ�񷨶�ͼ�����Log���Ӽ��
[BW1,thresh1] = edge(I,'log');
% ���ص�ǰLog���ӱ�Ե������ֵ
disp('Log�����Զ�ѡ�����ֵΪ��')
disp(thresh1)
subplot(121),imshow(BW1);
title('�Զ���ֵ��Log���ӱ�Ե���')
BW1 = edge(I,'log',0.005);
% ����ֵΪ0.005��ͼ�����Log���Ӽ��
subplot(122),imshow(BW1);
title('��ֵΪ0.005��Log���ӱ�Ե���')
h=fspecial('gaussian',5);
% ��Ƹ�˹�˲���
[BW2,thresh2] = edge(I,'zerocross',[],h);
% ���ص�ǰ�㽻�����Ե������ֵ
disp('�㽻�����Զ�ѡ�����ֵΪ��')
disp(thresh2)
figure,subplot(121),imshow(BW2);
title('�Զ���ֵ���㽻���Ե���')
BW2 = edge(I,'zerocross',0.03,h);
% ����ֵΪ0.03��ͼ������㽻����
subplot(122),imshow(BW2);
title('��ֵΪ0.03���㽻���Ե���')
