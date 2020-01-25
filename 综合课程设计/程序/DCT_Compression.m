%DCT压缩--------------------------------------------------------------------
function DCT_Compression(handles)
% 设置压缩比，cr=0.5为2:1压缩；cr=0.1250为8:1压缩
% fprintf('请设置压缩比为：1、0.5、0.25或0.125\n')
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
data1=get(handles.edit9,'string');
data1_num=str2num(data1);
cr=data1_num;
initialimage = I1;
initialimage = double(initialimage)/255;
[mm,nn]=size(initialimage);
figure;
subplot(121);
imshow(initialimage);
title('原始图像');
disp('原图像大小：')
whos('initialimage');
tic;
%对图像进行DCT变换
t = dctmtx(8);
dctcoe = blkproc(initialimage, [8 8], 'P1*x*P2', t, t');
%将DCT变换后的矩阵转换成列，并按升序排列
coevar = im2col(dctcoe, [8 8], 'distinct');
coe = coevar;
[y, ind] = sort(coevar);
[m, n] = size(coevar);
%舍去不重要的系数
snum = 64-64 * cr;
for i = 1:n
    coe(ind(1:snum), i) = 0;
end
%把列变换为二维矩阵
b2 = col2im(coe, [8 8], [mm nn], 'distinct');
t1=toc;
xl1=sprintf('图像编码用时%4.2f秒\n图像压缩比为：%4.1f：1\n',t1,1/cr);
xlabel(xl1);

%逆DCT变换

tic;
i2 = blkproc(b2, [8 8], 'P1*x*P2', t', t);
t2=toc;

PSNR = psnr(initialimage, i2);
xl2=sprintf('图像解码用时%4.2f秒\n两幅图像的峰值信噪比为：%4.2f dB\n\n',t2,PSNR);
subplot(122);
imshow(i2);
title('解码图像');xlabel(xl2);
disp('解码图像大小：')
whos('i2');
end

