%ͼ�������ж�
[X,map]=imread('canoe.tif');
Image_Type_YN1=isind(X);
%����ͼ��Ķ�ȡ���ж�
I = imread('moon.tif'); 
Image_Type_YN2=isgray(I);
%�Ҷ�ͼ��Ķ�ȡ���ж�
RGB=imread('airplane.tif'); 
Image_Type_YN3=isrgb(RGB);
%���ɫͼ��Ķ�ȡ���ж�
BW = imread('couple.tif');
Image_Type_YN4=isbw(BW);
%��ֵͼ��Ķ�ȡ���ж�
disp([Image_Type_YN1   Image_Type_YN2   Image_Type_YN3   Image_Type_YN4]);
%����ֵ�����MATLAB��������ʾ
