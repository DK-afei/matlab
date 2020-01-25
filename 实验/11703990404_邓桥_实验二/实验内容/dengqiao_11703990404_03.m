I= imread('lax.tif');
J = histeq(I, 256);
subplot(2,2,1);
imshow(uint8(I));
title('(a)原始图像')
subplot(2,2,2);
imshow(uint8(J));
title('(b)均衡后的图像')