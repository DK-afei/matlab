%��ͼ��ÿһ������������䰵
I=imread('lax.tif');
subplot(221);       %��2��2�еĵ�1��λ��
imshow(I);         %��ʾrich.tif�Ҷ�ͼ��
title('ԭ�Ҷ�ͼ��')
J=I;              %���¶���һ��ͼ�񣬸�ͼ����ʱ��I��ͬ
K=I;
k=0;
add=50;           %ͼ�����ȵĸı���
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
title('����ͼ��')
subplot(223);
imshow(K);
title('�䰵ͼ��')
