I=imread('lena.bmp');
I1=I*2;
I2=I/2;
figure;
subplot(3,3,1);
imshow(I);
title('ԭͼ');
subplot(3,3,2);
imshow(I1);
title('����2��');
subplot(3,3,3);
imshow(I2);
title('����1/2��');
A=double(I);
B=40*(log(A+1));
I3=uint8(B);
subplot(3,3,4);
imshow(I3);
title('����');
I_D=double(I);
C=I_D/255;
I4=uint8(255*(C.^0.7));
subplot(3,3,5);
imshow(I4);
title('���ɱ任����=0.5');
I5=uint8(255*(C.^1.3));
subplot(3,3,6);
imshow(I5);
title('���ɱ任����=1.5');

%��ͼ��ÿһ������������䰵
I=imread('lena.bmp');
subplot(337);       %��2��2�еĵ�1��λ��
imshow(I);         %��ʾrich.tif�Ҷ�ͼ��
title('ԭ�Ҷ�ͼ��')
J=I;              %���¶���һ��ͼ�񣬸�ͼ����ʱ��I��ͬ
K=I;
k=0;
add=100;           %ͼ�����ȵĸı���
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
title('����ͼ��')
subplot(339);
imshow(K);
title('�䰵ͼ��')
