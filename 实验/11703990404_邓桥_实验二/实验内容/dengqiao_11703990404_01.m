I=imread('airplane.tif');
figure(1);
subplot(2,2,1);
imshow(uint8(I));
title('(a)原始图像')
subplot(2,2,2);
imshow(uint8(I));
title('(b)原始图像')
I=double(I);
h=size(I);
I_fliplr(1:h(1),1:h(2),1:h(3))=I(1:h(1),h(2):-1:1,1:h(3));
I1=uint8(I_fliplr);
subplot(2,2,3);
imshow(I1);
title('(c)水平镜像变化')
I_flipud(1:h(1),1:h(2),1:h(3))=I(h(1):-1:1,1:h(2),1:h(3));
I2=uint8(I_flipud);
subplot(2,2,4);
imshow(I2);
title('(d)垂直镜像变化')