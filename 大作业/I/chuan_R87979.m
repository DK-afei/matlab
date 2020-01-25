%% 1.Ԥ����ͼƬ
fn = '��R87979.jpg';
img=imread('I/��R87979.jpg');%����ԭʼͼ��
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

%% 2.ˮƽ�ʹ�ֱͶӰ(ȥ���������������)

I3=remove_extra_region(I2);
subplot(338);imshow(I3);title('ȥ���������������');
 

%% 3.ȥ�����±߿��í����ͳ�����������
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
%subplot(336);imshow(I4);title('ȥ�����±߿��í��');

%% 4.ȥ�����ұ߿�ͶӰ����
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
% subplot(337);imshow(I5);title('ȥ�����ұ߿�');
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
%subplot(337);imshow(I5);title('ȥ�����ұ߿�');

%% 5.ȥ���ַ����ұ�����ͶӰ����
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
subplot(339);imshow(I6);title('�ַ�����');


%% 6.�ָ��ַ�����ֱͶӰ����
[height, Twidth] = size(I6);
Cwidth = Twidth* 15/116;  % ��һ�ַ����
Cspace = Twidth* 1/116;  % �ַ����
SecThspace = Twidth* 7/116;  % �ڶ����͵������ַ����
projection = sum(I6, 1);
 
figure;
for i=1:7
    if i == 1
        k =(floor(Twidth - Cwidth )); % �л������һ���ַ���ʼ��
        k=k-1;
    else  % ������������ɨ��
        k =(floor(k - Cwidth - Cspace)); % �л��ַ�����ʼ��
    end
    
    % �����������һ����
    if k <= 0
        k=1;
    end
    
    % ȡ��ǰ�ַ�
    fprintf('��%d�ַ���ʼ�еĴ��λ��:%d \n', i,k);
    fprintf('��ͶӰֵ:%d \n', projection(1, k));
    character = I6(:, k:ceil(k+Cwidth));
    subplot(178-i);imshow(character);
    
    % ���浱ǰ�ַ�
    j = 8-i;
    cn = strcat(fn(j), '.jpg');
    fprintf('%s \n', cn);
    char = imresize(character,[32, 16],'bilinear');
    %imwrite(char,['�ַ�/��R87979/',cn]);
    
    % �ڶ����͵������ַ�֮��Ŀո����⴦��
    if i == 5
        k = k - SecThspace + Cspace;
    end
end
 
 

