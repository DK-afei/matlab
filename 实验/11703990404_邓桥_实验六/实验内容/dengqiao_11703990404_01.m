%Laplacian算子模板
l=[0 1 0 ;1 -4 1; 0 1 0]
%Prewitt算子模板
p1=[-1 -1 -1;0 0 0;1 1 1];
p2=[-1 0 1; -1 0 1;-1 0 1];
%sobel算子模板
s1=[-1 -2 -1;0 0 0;1 2 1];
s2=[-1 0 1;-2 0 2;-1 0 1];
% Roberts算子模板
r1=[-1 0;0 1]
r2=[0 -1;1 0]
%点检测算子模板
e=[-1 -1 -1;-1 8 -1;-1 -1 -1]


I = imread('lena.bmp');
x=I;
y=I;

subplot(241)
imshow(I)
title('原始图像')
I=im2double(I); 
P0=conv2(I,l)
for i=2:511
    for j=2:511
        if(P0(i,j)>0.15)
            P0(i,j)=255;
        else 
            P0(i,j)=0;
        end
    end
end
subplot(242)
imshow(P0)
title('Laplace')

P1=conv2(I,p1)
P2=conv2(I,p2)
P3=P2+P2
for i=2:511
    for j=2:511
        if(P3(i,j)>0.30)
            P3(i,j)=255;
        else 
            P3(i,j)=0;
        end
    end
end

subplot(243);
imshow(P3);
title('Prewitt');

P3=conv2(I,s1)
P4=conv2(I,s2)
P5=P4+P3
for i=2:511
    for j=2:511
        if(P5(i,j)>0.30)
            P5(i,j)=255;
        else 
            P5(i,j)=0;
        end
    end
end
subplot(244);
imshow(P5);
title('Sobel')

P6=conv2(I,r1)
P7=conv2(I,r2)
P8=P6+P7
for i=2:511
    for j=2:511
        if(P8(i,j)>0.10)
            P8(i,j)=255;
        else 
            P8(i,j)=0;
        end
    end
end
subplot(245);
imshow(P8);
title('Roberts')


P9=conv2(I,e)
for i=2:511
    for j=2:511
        if(P9(i,j)>0.30)
            P9(i,j)=255;
        else 
            P9(i,j)=0;
        end
    end
end
subplot(246);
imshow(P9);
title('log')


P10=edge(y,'canny')
subplot(248);
imshow(P10);
title('canny')

P11=edge(x,'zerocross')
subplot(247);
imshow(P11);
title('零交叉')



