%已知一个30×30大小的二值图像，在图像中间有个长为5高为20的白色区%域，其它为黑色。对这幅图进行傅立叶变换分析（主要应用FFT算法）
f = zeros(30,30);
f(5:24,13:17) = 1;
%定义图像数组
imshow(f,'notruesize')
F = fft2(f);
%二维傅立叶变换（FFT算法）
mesh(fftshift(abs(F)));
%绘制频谱图
F2 = fftshift(log(abs(F)));
imshow(F2,[-1 5],'notruesize'); 
%显示频谱图像，频谱的零频率系数被移到频谱中间
colormap(jet); colorbar
%在上面的变换前的矩阵没有被填充，下面比较一下填充矩阵填充后的情况。
F = fft2(f,256,256);
%在变换前f被用0填充成256×256的矩阵，变换后的矩阵大小也是256×256
imshow(fftshift(log(abs(F))),[-1 5]); 
colormap(jet); colorbar
