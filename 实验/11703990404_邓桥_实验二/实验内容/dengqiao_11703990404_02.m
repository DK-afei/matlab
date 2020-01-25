 f= imread('lax.tif');

[m,n] = size(f);

f1 = im2uint8(ones(m,n));

%得到输入图像的直方图h

h = imhist(f);

l = length(h);

%概率密度PDF

PDF = h/numel(f);

%分布函数CDF

CDF = cumsum(PDF);

%取整扩展，得到均衡化之后的灰度分布直方图

j = CDF.*256;

%由于灰度级数为1-256之间的整数，故需对扩展之后的灰度灰度级数取整才有意义，

%得到的J矩阵为1X256大小，表示扩展之前的灰度级数，其中每个级数对应元素的

%值为该灰度级数扩展后的灰度级数值。如J（3）=24，表示原始灰度直方图为3灰度值

%的地方经灰度扩展后其灰度值为24

J = round(j);

%将扩展后的灰度级数对应映射到图片中

for i=1:l%l=256

    nn = find(J==i);%找出扩展后的灰度级数对应的扩展前的灰度级数

    L = length(nn);

    for k=1:L

    nn1 = find(f==(nn(k)-1));%再找到扩展前的灰度级数对应的像素点，

    f1(nn1)=i;              %并将像素点对应灰度值值置为扩展后的灰度值

    end

end
subplot(2,2,1);
imshow(uint8(f));
title('(a)原始图像')
subplot(2,2,2);
imshow(uint8(f1));
title('(b)均衡后的图像')

