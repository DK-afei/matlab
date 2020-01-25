global im;
I1 = imresize(im,0.5);
imshow(I1);	

%%
%放大*****************************************************************************
I=I1;
Igray=rgb2gray(I); %rgb2gray函数是将彩色图像转为灰度图像
figure,imshow(Igray);
[a,b]=size(Igray);
m=1.5;
J=zeros(a*m,b*m); %建立一个放大后尺寸的矩阵J
for row1=1:a*m
    for col1=1:b*m  %遍历J中所有的点
        row=round(row1/m);  %round函数为四舍五入的意思，对于近邻插值来说，他的原理就是J中的点对应原始图像Igray中的点，哪个近就取哪个。
        col=round(col1/m);
        if row<1
            row=1;  
        end   %原因是row1小于m时，row=0，但对于Igray矩阵来说行数不可能为0，因此人为地定义row1小于m时，row=1。
        if row>a
            row=a;
        end   %其实我在运行程序时实验了一下，发现如果将row大于a和col大于b两个if都去掉，是不影响输出结果的，因此我也不太清楚是否这个是毫无用处的。但小于1的两种情况是必要的，否则将显示【下标索引必须为正整数类型或逻辑类型。】。
        if col<1
            col=1;
        end
        if col>b
            col=b
        end
        J(row1,col1)=Igray(row,col);  %在规定完row、col的特殊情况后，就可以将Igray的值赋给J。
    end
end
%J=uint8(J);  %注意是uint！我之前总以为是unit。。。必须要改一下J的类型，否则J是double类型，将无法显示。原因是J图像矩阵是经过了运算的，他的数据类型会从uint变为double型。如果直接运行imshow(I)，我们会发现显示的是一个白色的图像。这是因为imshow()显示图像时对double型是认为在0~1范围内，即大于1时都是显示为白色，而imshow显示uint8型时是0~255范围。而经过运算的范围在0-255之间的double型数据就被不正常得显示为白色图像了。
figure,imshow(J,[]);
%如果有J=uint8(J)程序，就可以直接imshow(J)。imshow(I,[])是为了自动调整数据的范围以便于显示 
%%

%%缩小**********************************************************************

I2 = rgb2gray(I1);
a=I2;
mul=0.5;
type=1;
%****************************************************
%a:输入图像灰度值
%mul:缩放倍数
%type:1表示最邻近法，2表示双极性插值法
%画出缩放后图像并返回其灰度值
%****************************************************
[m,n]=size(a);
m1=m*mul;n1=n*mul;
%****************************************************
if type==1
for i=1:m1
for j=1:n1;
b(i,j)=a(round(i/mul),round(j/mul));
end
end
elseif type==2
for i=1:m1-1
for j=1:n1-1;
u0=i/mul;v0=j/mul;
u=round(u0);v=round(v0);
s=u0-u;t=v0-v;
b(i,j)=(a(u+1,v)-a(u,v))*s+(a(u,v+1)-a(u,v))*t+(a(u+1,v+1)+a(u,v)-a(u,v+1)-a(u+1,v))*s*t+a(u,v);
end
end
end
%*****************************************************
b=uint8(b);
imshow(b);
title('处理后图像');

%%
%平移**************************************************************************
imshow(im);
init=rgb2gray(im);
% 读取图像
[R, C] = size(init); % 获取图像大小
res = zeros(R, C); % 构造结果矩阵。每个像素点默认初始化为0（黑色）
delX = -200; % 平移量X
delY = 0; % 平移量Y
tras = [1 0 delX; 0 1 delY; 0 0 1]; % 平移的变换矩阵 

for i = 1 : R
    for j = 1 : C
        temp = [i; j; 1];
        temp = tras * temp; % 矩阵乘法
        x = temp(1, 1);
        y = temp(2, 1);
        % 变换后的位置判断是否越界
        if (x <= R) & (y <= C) & (x >= 1) & (y >= 1)
            res(x, y) = init(i, j);
        end
    end
end;
imshow(uint8(res)); % 显示图像
%%
%旋转************************************************************************
% 读入图片
im = imread('F:\大三\数字图像处理\课程设计\图片\bmp\airplane.bmp');
 imshow(im);
% 求出旋转矩阵
a = 190 / 180 * pi;
R = [cos(a), -sin(a); sin(a), cos(a)];
R = R'; % 求出旋转矩阵的逆矩阵进行逆向查找
 
% 计算原图大小
sz = size(im1);
h = sz(1);
w = sz(2);
ch = sz(3);
c1 = [h; w] / 2;
 
% 计算显示完整图像需要的画布大小
hh = floor(w*sin(a)+h*cos(a))+1;
ww = floor(w*cos(a)+h*sin(a))+1;
c2 = [hh; ww] / 2;
 
% 初始化目标画布
im2 = uint8(ones(hh, ww, 3)*128);
for k = 1:ch
    for i = 1:hh
       for j = 1:ww
          p = [i; j];
          pp = (R*(p-c2)+c1);
          mn = floor(pp);
          ab = pp - mn;
          a = ab(1);
          b = ab(2);
          m = mn(1);
          n = mn(2);
          % 线性插值方法
          if (pp(1) >= 2 && pp(1) <= h-1 && pp(2) >= 2 && pp(2) <= w-1)
             im2(i, j, k) = (1-a)*(1-b)*im(m, n, k) + a*(1-b)*im(m+1, n, k)...
                          + (1-a)*b*im(m, n, k)     + a*b*im(m, n, k);
          end
       end
    end
end
% 显示图像
figure;
imshow(im2);
%%
%翻转************************************************************************
clc;
I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\airplane.bmp'); 
%I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
[ROW,COL,DIM] = size(I);
Ih = uint8(zeros(ROW, COL,DIM));%Horizontal mirroring
Iv = uint8(zeros(ROW, COL,DIM));%Vertical mirroring
Ihv = uint8(zeros(ROW, COL,DIM));

%水平镜像
for i =1:ROW
    for j=1:COL
        for k=1:DIM
        x = i;
        y = COL-j+1;
        z = k;
        Ih(x,y,z) =I(i,j,k);
        end
    end
end

%垂直镜像
for i =1:ROW
    for j=1:COL
        for k=1:DIM
        x = ROW-i+1;
        y = j;
        z = k;
        Iv(x,y,z) =I(i,j,k);
        end
    end
end

%水平垂直镜像
for i =1:ROW
    for j=1:COL
        for k=1:DIM
        x = ROW-i+1;
        y = COL-j+1;
        z = k;
        Ihv(x,y,z) =I(i,j,k);
        end
    end
end
figure
subplot(221),imshow(I);title('原图');
subplot(222),imshow(Ih);title('水平镜像');
subplot(223),imshow(Iv);title('垂直镜像');
subplot(224),imshow(Ihv);title('水平垂直镜像');
%%****************************************************************************
%灰度变换
I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
%I1=imadjust(I,[0.2,0.5],[0,1]);
I1=log(double(I))*20;
figure,
imshow(I);
figure,
imshow(uint8(I1));


I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
I1=I*2;
I2=I/2;
figure;
subplot(3,3,1);
imshow(I);
title('原图');
subplot(3,3,2);
imshow(I1);
title('线性2倍');
subplot(3,3,3);
imshow(I2);
title('线性1/2倍');
A=double(I);
B=40*(log(A+1));
I3=uint8(B);
subplot(3,3,4);
imshow(I3);
title('对数');
I_D=double(I);
C=I_D/255;
I4=uint8(255*(C.^0.7));
subplot(3,3,5);
imshow(I4);
title('幂律变换—γ=0.5');
I5=uint8(255*(C.^1.3));
subplot(3,3,6);
imshow(I5);
title('幂律变换—γ=1.5');


%将图像每一个像素增亮或变暗
I=imread('lena.bmp');
subplot(337);       %在2行2列的第1个位置
imshow(I);         %显示rich.tif灰度图像
title('原灰度图像')
J=I;              %重新定义一副图像，该图像暂时与I相同
K=I;
k=0;
add=100;           %图像亮度的改变量
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
title('增亮图像')
subplot(339);
imshow(K);
title('变暗图像')

[m,n]=size(I);
for i=1:m;
    for j=1:n;
    end
end

%%
%旋转
I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\airplane.bmp'); 
%I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
degree=30;
B=I;
[r,c,d]=size(B);                                                      %获取输入图像B的行r、列c和通道数d,为了旋转彩色图像所以有必要得到通道数d
nH=round(r*abs(cosd(degree))+c*abs(sind(degree)));                    %旋转图像后得到的新高度，“round()函数四舍五入“
nW=round(c*abs(cosd(degree))+r*abs(sind(degree)));                    %旋转图像后得到的新宽度
A=zeros(nH,nW,d);                                                     %定义生成目标图像的行列以及通道数
M1=[1 0 0;0 -1 0;-0.5*nW 0.5*nH 1 ];                                  %坐标系变换矩阵M1
M2=[cosd(degree) -sind(degree) 0;sind(degree) cosd(degree) 0;0 0 1];  %角度旋转变换矩阵M2，我用的是顺时针方向
M3=[1 0 0;0 -1 0;0.5*c 0.5*r 1];                                      %坐标系变换矩阵M3
    for i=1:nW
        for j=1:nH
            temp=[i j 1]*M1*M2*M3;                                    %得到旋转后的矩阵temp
            y=temp(1,2);                                              %y取矩阵temp的第一行第二列,y对应j，为高度
            x=temp(1,1);                                              %x取矩阵temp的第一行第一列,x对应i，为宽度
            y=round(y);                                               %y四舍五入取整
            x=round(x);                                               %x四舍五入取整
           if(x>=1&&x<=c)&&(y>=1&&y<=r)                               %判断的得到的(x,y)点是否在原图像上
               A(j,i,:)=B(y,x,:);                                     %将原图像的像素点赋值给对应的旋转后图像上的点
           end                                                        %（”有人疑惑为啥不是A(i,j,:)=B(x,y,:);因为i,x对应的是列，即宽，而j,y对应的是行，即高“），我这里以x为横坐标，y为竖向纵坐标
        end
    end
imshow(B);
figure,
imshow(A);



%%
                              %定义旋转函数，degree为要旋转的角度
I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\airplane.bmp'); 
%I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
imshow(I);
II=imrotate(I,30);
imshow(II);
degree=240;
B=I;
[r,c,d]=size(B);                                                      %获取输入图像B的行r、列c和通道数d,为了旋转彩色图像所以有必要得到通道数d
nH=round(r*abs(cosd(degree))+c*abs(sind(degree)));                    %旋转图像后得到的新高度，“round()函数四舍五入“
nW=round(c*abs(cosd(degree))+r*abs(sind(degree)));                    %旋转图像后得到的新宽度
A=zeros(nH,nW,d);                                                     %定义生成目标图像的行列以及通道数
M1=[1 0 0;0 -1 0;-0.5*nW 0.5*nH 1 ];                                  %坐标系变换矩阵M1
M2=[cosd(degree) -sind(degree) 0;sind(degree) cosd(degree) 0;0 0 1];  %角度旋转变换矩阵M2，我用的是顺时针方向
M3=[1 0 0;0 -1 0;0.5*c 0.5*r 1];                                      %坐标系变换矩阵M3
    for i=1:nW
        for j=1:nH
            temp=[i j 1]*M1*M2*M3;                                    %得到旋转后的矩阵temp
            y=temp(1,2);                                              %y取矩阵temp的第一行第二列,y对应j，为高度
            x=temp(1,1);                                    
            y=round(y);                                               %y四舍五入取整
            x=round(x);                                               %x四舍五入取整
           if(x>=1&&x<=c)&&(y>=1&&y<=r)                               %判断的得到的(x,y)点是否在原图像上
               A(j,i,:)=B(y,x,:);                                     %将原图像的像素点赋值给对应的旋转后图像上的点
           end                                                        %（”有人疑惑为啥不是A(i,j,:)=B(x,y,:);因为i,x对应的是列，即宽，而j,y对应的是行，即高“），我这里以x为横坐标，y为竖向纵坐标
        end
    end
figure,imshow(uint8(A));
 %%
                              %定义旋转函数，degree为要旋转的角度
I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\airplane.bmp'); 
%I=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
% imshow(I);
% II=imrotate(I,30);
% imshow(II);
degree=30;
B=rgb2gray(I);
[r,c]=size(B);                                                      %获取输入图像B的行r、列c和通道数d,为了旋转彩色图像所以有必要得到通道数d
nH=round(r*abs(cosd(degree))+c*abs(sind(degree)));                    %旋转图像后得到的新高度，“round()函数四舍五入“
nW=round(c*abs(cosd(degree))+r*abs(sind(degree)));                    %旋转图像后得到的新宽度
A=zeros(nH,nW);                                                     %定义生成目标图像的行列以及通道数
M1=[1 0 0;0 -1 0;-0.5*nW 0.5*nH 1 ];                                  %坐标系变换矩阵M1
M2=[cosd(degree) -sind(degree) 0;sind(degree) cosd(degree) 0;0 0 1];  %角度旋转变换矩阵M2，我用的是顺时针方向
M3=[1 0 0;0 -1 0;0.5*c 0.5*r 1];                                      %坐标系变换矩阵M3
    for i=1:nW
        for j=1:nH
            temp=[i j 1]*M1*M2*M3;                                    %得到旋转后的矩阵temp
            y=temp(1,2);                                              %y取矩阵temp的第一行第二列,y对应j，为高度
            x=temp(1,1);                                    
            y=round(y);                                               %y四舍五入取整
            x=round(x);                                               %x四舍五入取整
           if(x>=1&&x<=c)&&(y>=1&&y<=r)                               %判断的得到的(x,y)点是否在原图像上
               A(j,i)=B(y,x);                                     %将原图像的像素点赋值给对应的旋转后图像上的点
           end                                                        %（”有人疑惑为啥不是A(i,j,:)=B(x,y,:);因为i,x对应的是列，即宽，而j,y对应的是行，即高“），我这里以x为横坐标，y为竖向纵坐标
        end
    end
figure,imshow(uint8(A));

I1=imread('F:\大三\数字图像处理\课程设计\图片\bmp\plane.bmp'); 
I2=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
imshow(I1);
if(size(I1)==size(I2))
    a=2
end
[m n]=size(I1);
I3=zeros(m,n);
    for i=1:m
        for j=1:n
            I3(i,j)=I1(i,j);
        end
    end
    imshow(I3);
    
    %%
I1=imread('F:\大三\数字图像处理\课程设计\图片\bmp\plane.bmp'); 
I2=imread('F:\大三\数字图像处理\课程设计\图片\bmp\camera.bmp'); 
imshow(I1);

image=I1;
[width,height,z]=size(image);
subplot(1,2,1);
imshow(image);
title('原图');
av=0;
std=0.1;
u1=rand(width,height);
u2=rand(width,height);
x=std*sqrt(-2*log(u1)).*cos(2*pi*u2)+av;
result1=double(image)/255+x;
result1=uint8(255*result1);
subplot(1,2,2);
imshow(result1);
title('加入均值为0，标准差为0.1的高斯噪声后');



J = imnoise(I2,'salt & pepper',0.02);
%添加椒盐噪声
imshow(J);

img=I1;
degree=90;
%获取图片信息 注意三通道获取完 即定义三个变量
[m,n,dep]=size(img);

%计算出旋转之后，形成一个大矩形的长宽 可以看效果图
rm=round(m*abs(cosd(degree))+n*abs(sind(degree)));
rn=round(m*abs(sind(degree))+n*abs(cosd(degree)));

%定义一个新矩阵，三通道的，存储新图片的信息
newimage=zeros(rm,rn,dep);

%坐标变换 分三步 
m1=[1,0,0;0,1,0;-0.5*rm,-0.5*rn,1];
m2=[cosd(degree),sind(degree),0;-sind(degree),cosd(degree),0;0,0,1];
m3=[1,0,0;0,1,0;0.5*m,0.5*n,1];

%利用循环，对每一个像素点进行变换
for i=1:rm
    for j=1:rn
        tem=[i j 1];
        tem=tem*m1*m2*m3;
        x=tem(1,1);
        y=tem(1,2);
        x=round(x);
        y=round(y);
        if(x>0&&x<=m)&&(y>0&&y<=n)
        newimage(i,j,:)=img(x,y,:);
        end
        end
end
imshow(uint8(newimage));