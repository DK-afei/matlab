im=imread('II/鄂K9M555.jpg');
figure%画一个窗口
subplot(4,1,1),imshow(im);%4 1 1 意思是这个窗口有四行，一列  最后的1的意思是这是第一个
title('原图')
im=imadjust(im,[0,0.3],[0.9,1])
gray=rgb2gray(im);%灰度化
level=graythresh(gray);%求阈值
bw=im2bw(gray,level);%转化为二值图像 
set(gcf,'Name','车牌分割')%窗口标题
%bw=bwareaopen(bw,30);
subplot(4,1,2),imshow(bw);%两行一列
title('二值化');


[m,n]=size(bw);%求二值化后车牌的长和宽，m是高，n是宽
% 求垂直投影


%对车牌进行再处理，把边缘噪声背景变为黑色，便于下一步的分割以及识别。
for x=1:m%对图片从上往下进行扫描
    count=0;
    for z=1:n-1
        if bw(x,z)*bw(x,z+1)==0
            if bw(x,z)==1 || bw(x,z+1)==1
                count=count+1;
            end
        
        end
    end
    if count<14          %跳变次数小于十四次（跳变就是从0到1或者1到0）                就认为是噪声， 把他变为背景色：黑色
        bw(x,1:n)=0;
    end;
end
%下面进行左右两边的噪声消除


for z=1:n*0.04%车牌左边
    temp=sum(bw(1:m,z));
    if temp<0.5*m
        bw(1:m,z)=0;
    end
end

for z=int32(0.96*n):n%车牌右边
    temp=sum(bw(1:m,z));
    if temp<0.5*m
        bw(1:m,z)=0;
    end
end


    
    

for y=1:n
     S(y)=sum(bw(1:m,y));
end

y=1:n;
bw=bwareaopen(bw,95);
subplot(413),imshow(bw);
title('背景处理');

for y=1:n
     S(y)=sum(bw(1:m,y));
end
y=1:n;


%=========================   字符分割   ============================
X=[];                               %用来存放水平分割线的横坐标
flag=0;
for j=1:size(bw,2)    
    sum_y=sum(bw(:,j));
    if logical(sum_y)~=flag         %列和有变化时，记录下此列
        X=[X j];
        flag=logical(sum_y);
    end
end
figure
for n=1:7                          
    char=bw(:,X(2*n-1):X(2*n)-1); %进行粗分割
    for i=1:size(char,1)            %这两个for循环对分割字符的上下进行裁剪
        if sum(char(i,:))~=0
            top=i;
            break
        end
    end
    for i=1:size(char,1)
        if sum(char(size(char,1)-i,:))~=0
            bottom=size(char,1)-i;
            break
        end
    end
    char=char(top:bottom,:);
    subplot(2,4,n);imshow(char);
    char=imresize(char,[32,16],'nearest'); %归一化为32*16的大小，以便模板匹配
    eval(strcat('Char_',num2str(n),'=char;'));  %将分割的字符放入Char_i中
end


        
    



