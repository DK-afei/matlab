clear;clc;
%% 1.���복��ͼƬ��Ԥ����
I=imread('I/��A99999.jpg');
subplot(231); imshow(I); title('ԭͼ');
I2 = rgb2gray(I);
subplot(232); imshow(I2); title('�Ҷ�');
%% 2.��Ե���ӡ�ֱ�߼�����ת
I3 = edge(I2,'Sobel','horizontal');
subplot(233);imshow(I3);title('��Ե���');
se = [1 1 1;1 1 1;1 1 1];
I4 = imdilate(I3,se);
subplot(234);imshow(I4);title('����');
[H,T,R] = hough(I4,'Theta',-89:89);
ccc = max(H);
[value, rot_theta] = max(ccc);
img_correction = imrotate(I , -1.2,'bilinear', 'loose');
subplot(235), imshow(img_correction);title('������');
%%
%%%%%% 1.Ԥ����ͼƬ
img=imadjust(img_correction,[0.1,0.5],[]);
subplot(331);imshow(img);title('ԭʼ');
I = rgb2gray(img);
[row,col]=size(I);
row=row*8;
col=col*8;
I=imresize(I,[row col],'bicubic');
subplot(332);imshow(I);title('�Ҷ�');
I1 = imbinarize(I);
subplot(333);imshow(I1);title('��ֵ');
I2 = bwareaopen(I1,20);
subplot(337);imshow(I2);title('ȥ����������');
%%
%%%%%% 2.ˮƽ�ʹ�ֱͶӰ(ȥ���������������)
I3=remove_extra_region(I2);
%subplot(335);imshow(I3);title('ȥ���������������');
 
%%
%%%%%% 3.ȥ�����±߿��í����ͳ�����������
%%% ��λ�е���ʼλ��(��1/3������ɨ����)
%%% ��λ�еĽ���λ��(��2/3������ɨ����)
diff_row = diff(I3,1,2);  % ǰһ�м���һ��
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
%subplot(336);imshow(I4);title('ȥ�����±߿��í��');

%%
%%%%%% 4.ȥ�����ұ߿�ͶӰ����
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
subplot(338);imshow(I5);title('ȥ������������');
%%
%%%%%% 5.ȥ���ַ����ұ�����ͶӰ����
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
subplot(339);imshow(I5);title('�ַ�����');
%%
%%%%%% 6.�ָ��ַ�����ֱͶӰ����
[height, Twidth] = size(I6);
Cwidth = Twidth*47/409;  % ��һ�ַ����
Cspace = Twidth*12/409;  % �ַ����
SecThspace = Twidth*34/409;  % �ڶ����͵������ַ����
projection = sum(I6, 1);
%subplot(339);stem(projection,'.',...
%     'MarkerFaceColor','w',...
%     'MarkerEdgeColor','w');
%title('�ַ���ֱ����ͶӰ');
 
figure;
for i=1:7
    if i == 1
        k = (floor(Twidth - Cwidth )); % �л������һ���ַ���ʼ��
        k=k-1;
    else  % ������������ɨ��
        k = (floor(k - Cwidth - Cspace)); % �л��ַ�����ʼ��
    end
    
    % �����������һ����
    if k <= 0
        k=1;
    end
    % ȡ��ǰ�ַ�
    fprintf('��%d�ַ���ʼ�еĴ��λ��:%d \n', i,k);
    fprintf('��ͶӰֵ:%d \n', projection(1, k));
    character = I6(:, k:ceil(k+Cwidth)+1);
    subplot(178-i);imshow(character);
    % �ڶ����͵������ַ�֮��Ŀո����⴦��
    if i == 5
        k = k - SecThspace + Cspace;
    end
end