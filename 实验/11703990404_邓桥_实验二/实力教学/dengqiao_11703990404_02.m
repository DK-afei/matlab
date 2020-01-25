%将图像每一个像素增亮或变暗
I=imread('lax.tif');
subplot(221);       %在2行2列的第1个位置
imshow(I);         %显示rich.tif灰度图像
title('原灰度图像')
J=I;              %重新定义一副图像，该图像暂时与I相同
K=I;
k=0;
add=50;           %图像亮度的改变量
for i=1:256;
   for j=1:256;
      k=double(I(i,j));
  	  if(k+add>255)J(i,j)=255;
      else   J(i,j)=uint8(k+add);
      end;
      if(k-add<0)K(i,j)=0;
      else   K(i,j)=uint8(k-add);
      end;
   end;
end;
subplot(222);
imshow(J);
title('增亮图像')
subplot(223);
imshow(K);
title('变暗图像')
