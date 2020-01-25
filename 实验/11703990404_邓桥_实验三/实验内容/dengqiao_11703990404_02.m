I=imread('lena.bmp');
I1=I*2;
I2=I/2;
figure;
subplot(3,3,1);
imshow(I);
title('原图');
subplot(3,3,2);
imshow(I1);
title('线性2倍');
subplot(3,3,3);
imshow(I2);
title('线性1/2倍');
A=double(I);
B=40*(log(A+1));
I3=uint8(B);
subplot(3,3,4);
imshow(I3);
title('对数');
I_D=double(I);
C=I_D/255;
I4=uint8(255*(C.^0.7));
subplot(3,3,5);
imshow(I4);
title('幂律变换—γ=0.5');
I5=uint8(255*(C.^1.3));
subplot(3,3,6);
imshow(I5);
title('幂律变换—γ=1.5');

%将图像每一个像素增亮或变暗
I=imread('lena.bmp');
subplot(337);       %在2行2列的第1个位置
imshow(I);         %显示rich.tif灰度图像
title('原灰度图像')
J=I;              %重新定义一副图像，该图像暂时与I相同
K=I;
k=0;
add=100;           %图像亮度的改变量
for i=1:512;
   for j=1:512;
      k=double(I(i,j));
  	  if(k+add>255)J(i,j)=255;
      else   J(i,j)=uint8(k+add);
      end;
      if(k-add<0)K(i,j)=0;
      else   K(i,j)=uint8(k-add);
      end;
   end;
end;
subplot(338);
imshow(J);
title('增亮图像')
subplot(339);
imshow(K);
title('变暗图像')
