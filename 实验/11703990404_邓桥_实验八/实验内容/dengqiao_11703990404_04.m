%С���任ѹ��ͼ��
clear
load wbarb;
%װ��ͼ��
subplot(221);
X=imread('lena.bmp')
[m n]=size(X);
image(X);
%��ʾͼ��
colormap(map)
title('ԭʼͼ��');
axis square
disp('ѹ��ǰͼ��X�Ĵ�С');
whos('X')
[c,s]=wavedec2(X,2,'bior3.7');
%��ͼ����С�����в�ֽ�
cal=appcoef2(c,s,'bior3.7',1);
%��ȡС���ֽ�ṹ�е�һ��ĵ�Ƶϵ���͸�Ƶϵ��
ch1=detcoef2('h',c,s,1);
%ˮƽ����
cv1=detcoef2('v',c,s,1);
%��ֱ����
cd1=detcoef2('d',c,s,1);
%б�߷���
a1=wrcoef2('a',c,s,'bior3.7',1);
h1=wrcoef2('h',c,s,'bior3.7',1);
v1=wrcoef2('v',c,s,'bior3.7',1);
d1=wrcoef2('d',c,s,'bior3.7',1);
%��Ƶ�ʳɷ��ع�
c1=[a1,h1;v1,d1];
subplot(222);
image(c1);
%��ʾ��Ƶ��Ϣ
axis square;
title ('�ֽ���Ƶ�͸�Ƶ��Ϣ');
%����ͼ��ѹ��
%����С���ֽ��һ���Ƶ��Ϣ
%���ȶԵ�һ����Ϣ������������
ca1=appcoef2(c,s,'bior3.7',1);
ca1=wcodemat(ca1,440,'mat',0);
ca1=0.5*ca1;
subplot(223);
image(ca1);
%�ı�ͼ��߶Ȳ���ʾ
colormap(map);
axis square;
title('��һ��ѹ��ͼ��');
disp('��һ��ѹ��ͼ��Ĵ�СΪ��');
whos('ca1')
disp('ͼ���ѹ���ȣ�')
disp(m*n/length(ca1))
ca2=appcoef2(c,s,'bior3.7',2);
%����С���ֽ�ڶ����Ƶ��Ϣ����ѹ��
ca2=wcodemat(ca2,440,'mat',0);
%���ȶԵڶ�����Ϣ������������
ca2=0.25*ca2;
%�ı�ͼ��߶Ȳ���ʾ
subplot(224);
image(ca2);
colormap(map);
axis square;
title('�ڶ���ѹ��ͼ��');
disp('�ڶ���ѹ��ͼ��Ĵ�СΪ��');
whos('ca2')
disp('ͼ���ѹ���ȣ�')
disp(m*n/length(ca2))