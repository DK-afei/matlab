im=imread('II/��K9M555.jpg');
figure%��һ������
subplot(4,1,1),imshow(im);%4 1 1 ��˼��������������У�һ��  ����1����˼�����ǵ�һ��
title('ԭͼ')
im=imadjust(im,[0,0.3],[0.9,1])
gray=rgb2gray(im);%�ҶȻ�
level=graythresh(gray);%����ֵ
bw=im2bw(gray,level);%ת��Ϊ��ֵͼ�� 
set(gcf,'Name','���Ʒָ�')%���ڱ���
%bw=bwareaopen(bw,30);
subplot(4,1,2),imshow(bw);%����һ��
title('��ֵ��');


[m,n]=size(bw);%���ֵ�����Ƶĳ��Ϳ�m�Ǹߣ�n�ǿ�
% ��ֱͶӰ


%�Գ��ƽ����ٴ����ѱ�Ե����������Ϊ��ɫ��������һ���ķָ��Լ�ʶ��
for x=1:m%��ͼƬ�������½���ɨ��
    count=0;
    for z=1:n-1
        if bw(x,z)*bw(x,z+1)==0
            if bw(x,z)==1 || bw(x,z+1)==1
                count=count+1;
            end
        
        end
    end
    if count<14          %�������С��ʮ�ĴΣ�������Ǵ�0��1����1��0��                ����Ϊ�������� ������Ϊ����ɫ����ɫ
        bw(x,1:n)=0;
    end;
end
%��������������ߵ���������


for z=1:n*0.04%�������
    temp=sum(bw(1:m,z));
    if temp<0.5*m
        bw(1:m,z)=0;
    end
end

for z=int32(0.96*n):n%�����ұ�
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
title('��������');

for y=1:n
     S(y)=sum(bw(1:m,y));
end
y=1:n;


%=========================   �ַ��ָ�   ============================
X=[];                               %�������ˮƽ�ָ��ߵĺ�����
flag=0;
for j=1:size(bw,2)    
    sum_y=sum(bw(:,j));
    if logical(sum_y)~=flag         %�к��б仯ʱ����¼�´���
        X=[X j];
        flag=logical(sum_y);
    end
end
figure
for n=1:7                          
    char=bw(:,X(2*n-1):X(2*n)-1); %���дַָ�
    for i=1:size(char,1)            %������forѭ���Էָ��ַ������½��вü�
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
    char=imresize(char,[32,16],'nearest'); %��һ��Ϊ32*16�Ĵ�С���Ա�ģ��ƥ��
    eval(strcat('Char_',num2str(n),'=char;'));  %���ָ���ַ�����Char_i��
end


        
    



