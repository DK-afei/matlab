%�г̱���
I =imread('lena.bmp');
%����ͼ��
[m n]=size(I);
J=[];
for i=1:m
    value=I(i,1);
    num=1;
    for j=2:n
        if I(i,j)==value
            num=num+1;
        else
            J=[J num value];
            num=1;
            value=I(i,j);
        end
    end    
    J=[J num value 0 0];
    %��ӵ����ж�λ 0 0
end
disp('ԭͼ���С��')
whos('I');
disp('ѹ��ͼ���С��')
whos('J');
disp('ͼ���ѹ���ȣ�')
disp(m*n/length(J))
