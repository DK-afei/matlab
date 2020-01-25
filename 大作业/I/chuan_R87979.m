%% 1.预处理图片
fn = '川R87979.jpg';
img=imread('I/川R87979.jpg');%输入原始图像
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

%% 2.水平和垂直投影(去掉车牌以外的区域)

I3=remove_extra_region(I2);
subplot(338);imshow(I3);title('去除车牌以外的区域');
 

%% 3.去掉上下边框和铆钉（统计跳变次数）
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
%subplot(336);imshow(I4);title('去除上下边框和铆钉');

%% 4.去除左右边框（投影法）
% plate_projection_v = sum(I4,1);
% for i=1:size(plate_projection_v, 2)
%     if plate_projection_v(1,i) == 0
%         plate.cola = i-4;
%         break;
%     end
% end
% for i=1:size(plate_projection_v, 2)
%     j = size(plate_projection_v, 2) - i + 1;
%     if plate_projection_v(1,j) == 0
%         plate.colb = j+45;
%         break;
%     end
% end
% I5 = I4(:,plate.cola:plate.colb);
% subplot(337);imshow(I5);title('去除左右边框');
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
%subplot(337);imshow(I5);title('去除左右边框');

%% 5.去除字符左右背景（投影法）
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
subplot(339);imshow(I6);title('字符车牌');


%% 6.分割字符（垂直投影法）
[height, Twidth] = size(I6);
Cwidth = Twidth* 15/116;  % 单一字符宽度
Cspace = Twidth* 1/116;  % 字符间距
SecThspace = Twidth* 7/116;  % 第二个和第三个字符间距
projection = sum(I6, 1);
 
figure;
for i=1:7
    if i == 1
        k =(floor(Twidth - Cwidth )); % 切换到最后一个字符起始列
        k=k-1;
    else  % 自右向左逐列扫描
        k =(floor(k - Cwidth - Cspace)); % 切换字符的起始列
    end
    
    % 对特殊情况置一处理
    if k <= 0
        k=1;
    end
    
    % 取当前字符
    fprintf('第%d字符起始列的大概位置:%d \n', i,k);
    fprintf('列投影值:%d \n', projection(1, k));
    character = I6(:, k:ceil(k+Cwidth));
    subplot(178-i);imshow(character);
    
    % 保存当前字符
    j = 8-i;
    cn = strcat(fn(j), '.jpg');
    fprintf('%s \n', cn);
    char = imresize(character,[32, 16],'bilinear');
    %imwrite(char,['字符/川R87979/',cn]);
    
    % 第二个和第三个字符之间的空格特殊处理
    if i == 5
        k = k - SecThspace + Cspace;
    end
end
 
 

