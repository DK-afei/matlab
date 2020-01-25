global im;
I1 = imresize(im,0.5);
imshow(I1);	

%%
%�Ŵ�*****************************************************************************
I=I1;
Igray=rgb2gray(I); %rgb2gray�����ǽ���ɫͼ��תΪ�Ҷ�ͼ��
figure,imshow(Igray);
[a,b]=size(Igray);
m=1.5;
J=zeros(a*m,b*m); %����һ���Ŵ��ߴ�ľ���J
for row1=1:a*m
    for col1=1:b*m  %����J�����еĵ�
        row=round(row1/m);  %round����Ϊ�����������˼�����ڽ��ڲ�ֵ��˵������ԭ�����J�еĵ��Ӧԭʼͼ��Igray�еĵ㣬�ĸ�����ȡ�ĸ���
        col=round(col1/m);
        if row<1
            row=1;  
        end   %ԭ����row1С��mʱ��row=0��������Igray������˵����������Ϊ0�������Ϊ�ض���row1С��mʱ��row=1��
        if row>a
            row=a;
        end   %��ʵ�������г���ʱʵ����һ�£����������row����a��col����b����if��ȥ�����ǲ�Ӱ���������ģ������Ҳ��̫����Ƿ�����Ǻ����ô��ġ���С��1����������Ǳ�Ҫ�ģ�������ʾ���±���������Ϊ���������ͻ��߼����͡�����
        if col<1
            col=1;
        end
        if col>b
            col=b
        end
        J(row1,col1)=Igray(row,col);  %�ڹ涨��row��col����������󣬾Ϳ��Խ�Igray��ֵ����J��
    end
end
%J=uint8(J);  %ע����uint����֮ǰ����Ϊ��unit����������Ҫ��һ��J�����ͣ�����J��double���ͣ����޷���ʾ��ԭ����Jͼ������Ǿ���������ģ������������ͻ��uint��Ϊdouble�͡����ֱ������imshow(I)�����ǻᷢ����ʾ����һ����ɫ��ͼ��������Ϊimshow()��ʾͼ��ʱ��double������Ϊ��0~1��Χ�ڣ�������1ʱ������ʾΪ��ɫ����imshow��ʾuint8��ʱ��0~255��Χ������������ķ�Χ��0-255֮���double�����ݾͱ�����������ʾΪ��ɫͼ���ˡ�
figure,imshow(J,[]);
%�����J=uint8(J)���򣬾Ϳ���ֱ��imshow(J)��imshow(I,[])��Ϊ���Զ��������ݵķ�Χ�Ա�����ʾ 
%%

%%��С**********************************************************************

I2 = rgb2gray(I1);
a=I2;
mul=0.5;
type=1;
%****************************************************
%a:����ͼ��Ҷ�ֵ
%mul:���ű���
%type:1��ʾ���ڽ�����2��ʾ˫���Բ�ֵ��
%�������ź�ͼ�񲢷�����Ҷ�ֵ
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
title('�����ͼ��');

%%
%ƽ��**************************************************************************
imshow(im);
init=rgb2gray(im);
% ��ȡͼ��
[R, C] = size(init); % ��ȡͼ���С
res = zeros(R, C); % ����������ÿ�����ص�Ĭ�ϳ�ʼ��Ϊ0����ɫ��
delX = -200; % ƽ����X
delY = 0; % ƽ����Y
tras = [1 0 delX; 0 1 delY; 0 0 1]; % ƽ�Ƶı任���� 

for i = 1 : R
    for j = 1 : C
        temp = [i; j; 1];
        temp = tras * temp; % ����˷�
        x = temp(1, 1);
        y = temp(2, 1);
        % �任���λ���ж��Ƿ�Խ��
        if (x <= R) & (y <= C) & (x >= 1) & (y >= 1)
            res(x, y) = init(i, j);
        end
    end
end;
imshow(uint8(res)); % ��ʾͼ��
%%
%��ת************************************************************************
% ����ͼƬ
im = imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\airplane.bmp');
 imshow(im);
% �����ת����
a = 190 / 180 * pi;
R = [cos(a), -sin(a); sin(a), cos(a)];
R = R'; % �����ת��������������������
 
% ����ԭͼ��С
sz = size(im1);
h = sz(1);
w = sz(2);
ch = sz(3);
c1 = [h; w] / 2;
 
% ������ʾ����ͼ����Ҫ�Ļ�����С
hh = floor(w*sin(a)+h*cos(a))+1;
ww = floor(w*cos(a)+h*sin(a))+1;
c2 = [hh; ww] / 2;
 
% ��ʼ��Ŀ�껭��
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
          % ���Բ�ֵ����
          if (pp(1) >= 2 && pp(1) <= h-1 && pp(2) >= 2 && pp(2) <= w-1)
             im2(i, j, k) = (1-a)*(1-b)*im(m, n, k) + a*(1-b)*im(m+1, n, k)...
                          + (1-a)*b*im(m, n, k)     + a*b*im(m, n, k);
          end
       end
    end
end
% ��ʾͼ��
figure;
imshow(im2);
%%
%��ת************************************************************************
clc;
I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\airplane.bmp'); 
%I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
[ROW,COL,DIM] = size(I);
Ih = uint8(zeros(ROW, COL,DIM));%Horizontal mirroring
Iv = uint8(zeros(ROW, COL,DIM));%Vertical mirroring
Ihv = uint8(zeros(ROW, COL,DIM));

%ˮƽ����
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

%��ֱ����
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

%ˮƽ��ֱ����
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
subplot(221),imshow(I);title('ԭͼ');
subplot(222),imshow(Ih);title('ˮƽ����');
subplot(223),imshow(Iv);title('��ֱ����');
subplot(224),imshow(Ihv);title('ˮƽ��ֱ����');
%%****************************************************************************
%�Ҷȱ任
I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
%I1=imadjust(I,[0.2,0.5],[0,1]);
I1=log(double(I))*20;
figure,
imshow(I);
figure,
imshow(uint8(I1));


I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
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

[m,n]=size(I);
for i=1:m;
    for j=1:n;
    end
end

%%
%��ת
I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\airplane.bmp'); 
%I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
degree=30;
B=I;
[r,c,d]=size(B);                                                      %��ȡ����ͼ��B����r����c��ͨ����d,Ϊ����ת��ɫͼ�������б�Ҫ�õ�ͨ����d
nH=round(r*abs(cosd(degree))+c*abs(sind(degree)));                    %��תͼ���õ����¸߶ȣ���round()�����������롰
nW=round(c*abs(cosd(degree))+r*abs(sind(degree)));                    %��תͼ���õ����¿��
A=zeros(nH,nW,d);                                                     %��������Ŀ��ͼ��������Լ�ͨ����
M1=[1 0 0;0 -1 0;-0.5*nW 0.5*nH 1 ];                                  %����ϵ�任����M1
M2=[cosd(degree) -sind(degree) 0;sind(degree) cosd(degree) 0;0 0 1];  %�Ƕ���ת�任����M2�����õ���˳ʱ�뷽��
M3=[1 0 0;0 -1 0;0.5*c 0.5*r 1];                                      %����ϵ�任����M3
    for i=1:nW
        for j=1:nH
            temp=[i j 1]*M1*M2*M3;                                    %�õ���ת��ľ���temp
            y=temp(1,2);                                              %yȡ����temp�ĵ�һ�еڶ���,y��Ӧj��Ϊ�߶�
            x=temp(1,1);                                              %xȡ����temp�ĵ�һ�е�һ��,x��Ӧi��Ϊ���
            y=round(y);                                               %y��������ȡ��
            x=round(x);                                               %x��������ȡ��
           if(x>=1&&x<=c)&&(y>=1&&y<=r)                               %�жϵĵõ���(x,y)���Ƿ���ԭͼ����
               A(j,i,:)=B(y,x,:);                                     %��ԭͼ������ص㸳ֵ����Ӧ����ת��ͼ���ϵĵ�
           end                                                        %���������ɻ�Ϊɶ����A(i,j,:)=B(x,y,:);��Ϊi,x��Ӧ�����У�������j,y��Ӧ�����У����ߡ�������������xΪ�����꣬yΪ����������
        end
    end
imshow(B);
figure,
imshow(A);



%%
                              %������ת������degreeΪҪ��ת�ĽǶ�
I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\airplane.bmp'); 
%I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
imshow(I);
II=imrotate(I,30);
imshow(II);
degree=240;
B=I;
[r,c,d]=size(B);                                                      %��ȡ����ͼ��B����r����c��ͨ����d,Ϊ����ת��ɫͼ�������б�Ҫ�õ�ͨ����d
nH=round(r*abs(cosd(degree))+c*abs(sind(degree)));                    %��תͼ���õ����¸߶ȣ���round()�����������롰
nW=round(c*abs(cosd(degree))+r*abs(sind(degree)));                    %��תͼ���õ����¿��
A=zeros(nH,nW,d);                                                     %��������Ŀ��ͼ��������Լ�ͨ����
M1=[1 0 0;0 -1 0;-0.5*nW 0.5*nH 1 ];                                  %����ϵ�任����M1
M2=[cosd(degree) -sind(degree) 0;sind(degree) cosd(degree) 0;0 0 1];  %�Ƕ���ת�任����M2�����õ���˳ʱ�뷽��
M3=[1 0 0;0 -1 0;0.5*c 0.5*r 1];                                      %����ϵ�任����M3
    for i=1:nW
        for j=1:nH
            temp=[i j 1]*M1*M2*M3;                                    %�õ���ת��ľ���temp
            y=temp(1,2);                                              %yȡ����temp�ĵ�һ�еڶ���,y��Ӧj��Ϊ�߶�
            x=temp(1,1);                                    
            y=round(y);                                               %y��������ȡ��
            x=round(x);                                               %x��������ȡ��
           if(x>=1&&x<=c)&&(y>=1&&y<=r)                               %�жϵĵõ���(x,y)���Ƿ���ԭͼ����
               A(j,i,:)=B(y,x,:);                                     %��ԭͼ������ص㸳ֵ����Ӧ����ת��ͼ���ϵĵ�
           end                                                        %���������ɻ�Ϊɶ����A(i,j,:)=B(x,y,:);��Ϊi,x��Ӧ�����У�������j,y��Ӧ�����У����ߡ�������������xΪ�����꣬yΪ����������
        end
    end
figure,imshow(uint8(A));
 %%
                              %������ת������degreeΪҪ��ת�ĽǶ�
I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\airplane.bmp'); 
%I=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
% imshow(I);
% II=imrotate(I,30);
% imshow(II);
degree=30;
B=rgb2gray(I);
[r,c]=size(B);                                                      %��ȡ����ͼ��B����r����c��ͨ����d,Ϊ����ת��ɫͼ�������б�Ҫ�õ�ͨ����d
nH=round(r*abs(cosd(degree))+c*abs(sind(degree)));                    %��תͼ���õ����¸߶ȣ���round()�����������롰
nW=round(c*abs(cosd(degree))+r*abs(sind(degree)));                    %��תͼ���õ����¿��
A=zeros(nH,nW);                                                     %��������Ŀ��ͼ��������Լ�ͨ����
M1=[1 0 0;0 -1 0;-0.5*nW 0.5*nH 1 ];                                  %����ϵ�任����M1
M2=[cosd(degree) -sind(degree) 0;sind(degree) cosd(degree) 0;0 0 1];  %�Ƕ���ת�任����M2�����õ���˳ʱ�뷽��
M3=[1 0 0;0 -1 0;0.5*c 0.5*r 1];                                      %����ϵ�任����M3
    for i=1:nW
        for j=1:nH
            temp=[i j 1]*M1*M2*M3;                                    %�õ���ת��ľ���temp
            y=temp(1,2);                                              %yȡ����temp�ĵ�һ�еڶ���,y��Ӧj��Ϊ�߶�
            x=temp(1,1);                                    
            y=round(y);                                               %y��������ȡ��
            x=round(x);                                               %x��������ȡ��
           if(x>=1&&x<=c)&&(y>=1&&y<=r)                               %�жϵĵõ���(x,y)���Ƿ���ԭͼ����
               A(j,i)=B(y,x);                                     %��ԭͼ������ص㸳ֵ����Ӧ����ת��ͼ���ϵĵ�
           end                                                        %���������ɻ�Ϊɶ����A(i,j,:)=B(x,y,:);��Ϊi,x��Ӧ�����У�������j,y��Ӧ�����У����ߡ�������������xΪ�����꣬yΪ����������
        end
    end
figure,imshow(uint8(A));

I1=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\plane.bmp'); 
I2=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
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
I1=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\plane.bmp'); 
I2=imread('F:\����\����ͼ����\�γ����\ͼƬ\bmp\camera.bmp'); 
imshow(I1);

image=I1;
[width,height,z]=size(image);
subplot(1,2,1);
imshow(image);
title('ԭͼ');
av=0;
std=0.1;
u1=rand(width,height);
u2=rand(width,height);
x=std*sqrt(-2*log(u1)).*cos(2*pi*u2)+av;
result1=double(image)/255+x;
result1=uint8(255*result1);
subplot(1,2,2);
imshow(result1);
title('�����ֵΪ0����׼��Ϊ0.1�ĸ�˹������');



J = imnoise(I2,'salt & pepper',0.02);
%��ӽ�������
imshow(J);

img=I1;
degree=90;
%��ȡͼƬ��Ϣ ע����ͨ����ȡ�� ��������������
[m,n,dep]=size(img);

%�������ת֮���γ�һ������εĳ��� ���Կ�Ч��ͼ
rm=round(m*abs(cosd(degree))+n*abs(sind(degree)));
rn=round(m*abs(sind(degree))+n*abs(cosd(degree)));

%����һ���¾�����ͨ���ģ��洢��ͼƬ����Ϣ
newimage=zeros(rm,rn,dep);

%����任 ������ 
m1=[1,0,0;0,1,0;-0.5*rm,-0.5*rn,1];
m2=[cosd(degree),sind(degree),0;-sind(degree),cosd(degree),0;0,0,1];
m3=[1,0,0;0,1,0;0.5*m,0.5*n,1];

%����ѭ������ÿһ�����ص���б任
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