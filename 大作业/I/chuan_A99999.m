clear;clc;
%% 1.载入车牌图片并预处理
I=imread('I/川A99999.jpg');
subplot(231); imshow(I); title('原图');
I2 = rgb2gray(I);
subplot(232); imshow(I2); title('灰度');
%% 2.边缘算子、直线检测和旋转
I3 = edge(I2,'Sobel','horizontal');
subplot(233);imshow(I3);title('边缘检测');
se = [1 1 1;1 1 1;1 1 1];
I4 = imdilate(I3,se);
subplot(234);imshow(I4);title('膨胀');
[H,T,R] = hough(I4,'Theta',-89:89);
ccc = max(H);
[value, rot_theta] = max(ccc);
img_correction = imrotate(I , -1.2,'bilinear', 'loose');
subplot(235), imshow(img_correction);title('矫正后');
%%
%%%%%% 1.预处理图片
img=imadjust(img_correction,[0.1,0.5],[]);
subplot(331);imshow(img);title('原始');
I = rgb2gray(img);
[row,col]=size(I);
row=row*8;
col=col*8;
I=imresize(I,[row col],'bicubic');
subplot(332);imshow(I);title('灰度');
I1 = imbinarize(I);
subplot(333);imshow(I1);title('二值');
I2 = bwareaopen(I1,20);
subplot(337);imshow(I2);title('去除孤立噪声');
%%
%%%%%% 2.水平和垂直投影(去掉车牌以外的区域)
I3=remove_extra_region(I2);
%subplot(335);imshow(I3);title('去除车牌以外的区域');
 
%%
%%%%%% 3.去掉上下边框和铆钉（统计跳变次数）
%%% 定位行的起始位置(从1/3处先上扫描行)
%%% 定位行的结束位置(从2/3处先下扫描行)
diff_row = diff(I3,1,2);  % 前一列减后一列
diff_row_sum = sum(abs(diff_row), 2);  
[rows, columns] = size(I3);
trows = ceil(rows*(1/3));
j = trows;
for i=1:trows
    if diff_row_sum(j,1)<10
        plate.rowa = j;
        break;
    end
    j = trows-i;
end
 
for i=2*trows:size(diff_row_sum,1)
    if diff_row_sum(i,1)<10
        plate.rowb = i;
        break;
    end
end
I4 = I3(plate.rowa:plate.rowb, :);
I4=remove_extra_region(I4);
I4 = bwareaopen(I4,20);
%subplot(336);imshow(I4);title('去除上下边框和铆钉');

%%
%%%%%% 4.去除左右边框（投影法）
plate_projection_v = sum(I4,1);
for i=1:size(plate_projection_v, 2)
    if plate_projection_v(1,i) == 0
        plate.cola = i;
        break;
    end
end 
 
for i=1:size(plate_projection_v, 2)
    j = size(plate_projection_v, 2) - i + 1;
    if plate_projection_v(1,j) == 0
        plate.colb = j;
        break;
    end
end
I5 = I4(:,plate.cola:plate.colb);
% I5=imcrop(I4,[0 0 111 111]);
subplot(338);imshow(I5);title('去除车牌以区域');
%%
%%%%%% 5.去除字符左右背景（投影法）
ppv1 = sum(I5,1);
for i=1:size(ppv1, 2)
    if ppv1(1,i) ~= 0
        pl.cola = i;
        break;
    end
end
 
for i=1:size(ppv1, 2)
    j = size(ppv1, 2) - i + 1;
    if ppv1(1,j) ~= 0
        pl.colb = j;
        break;
    end
end
I6 = I5(:,pl.cola:pl.colb);
subplot(339);imshow(I5);title('字符车牌');
%%
%%%%%% 6.分割字符（垂直投影法）
[height, Twidth] = size(I6);
Cwidth = Twidth*47/409;  % 单一字符宽度
Cspace = Twidth*12/409;  % 字符间距
SecThspace = Twidth*34/409;  % 第二个和第三个字符间距
projection = sum(I6, 1);
%subplot(339);stem(projection,'.',...
%     'MarkerFaceColor','w',...
%     'MarkerEdgeColor','w');
%title('字符垂直方向投影');
 
figure;
for i=1:7
    if i == 1
        k = (floor(Twidth - Cwidth )); % 切换到最后一个字符起始列
        k=k-1;
    else  % 自右向左逐列扫描
        k = (floor(k - Cwidth - Cspace)); % 切换字符的起始列
    end
    
    % 对特殊情况置一处理
    if k <= 0
        k=1;
    end
    % 取当前字符
    fprintf('第%d字符起始列的大概位置:%d \n', i,k);
    fprintf('列投影值:%d \n', projection(1, k));
    character = I6(:, k:ceil(k+Cwidth)+1);
    subplot(178-i);imshow(character);
    % 第二个和第三个字符之间的空格特殊处理
    if i == 5
        k = k - SecThspace + Cspace;
    end
end