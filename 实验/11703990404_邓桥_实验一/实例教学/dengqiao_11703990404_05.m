%图像类型判断
[X,map]=imread('canoe.tif');
Image_Type_YN1=isind(X);
%索引图像的读取与判断
I = imread('moon.tif'); 
Image_Type_YN2=isgray(I);
%灰度图像的读取与判断
RGB=imread('airplane.tif'); 
Image_Type_YN3=isrgb(RGB);
%真彩色图像的读取与判断
BW = imread('couple.tif');
Image_Type_YN4=isbw(BW);
%二值图像的读取与判断
disp([Image_Type_YN1   Image_Type_YN2   Image_Type_YN3   Image_Type_YN4]);
%返回值结果在MATLAB工作间显示
