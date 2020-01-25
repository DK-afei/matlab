%行程编码
I =imread('lena.bmp');
%棋盘图像
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
    %添加的行判断位 0 0
end
disp('原图像大小：')
whos('I');
disp('压缩图像大小：')
whos('J');
disp('图像的压缩比：')
disp(m*n/length(J))
