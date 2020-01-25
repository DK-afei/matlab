%ÀëÉ¢ÓàÏÒ±ä»»±àÂëÑ¹ËõÍ¼Ïñ
I = imread('lena.bmp');
[m n]=size(I);
I = im2double(I);
%Í¼Ïñ´æ´¢ÀàĞÍ×ª»»
T = dctmtx(8);
%ÀëÉ¢ÓàÏÒ±ä»»¾ØÕó
B = blkproc(I,[8 8],'P1*x*P2',T,T');
mask = [1   1   1   1   0   0   0   0
       1   1   1   0   0   0   0   0
       1   1   0   0   0   0   0   0
       1   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0
       0   0   0   0   0   0   0   0];
B2 = blkproc(B,[8 8],'P1.*x',mask);
I2 = blkproc(B2,[8 8],'P1*x*P2',T',T);
imshow(I)
title('Ô­Ê¼Í¼Ïñ')
figure;
imshow(I2)
title('Ñ¹ËõºóµÄÍ¼Ïñ')
disp('Í¼ÏñµÄÑ¹Ëõ±È£º')
disp(m*n/length(I2))

