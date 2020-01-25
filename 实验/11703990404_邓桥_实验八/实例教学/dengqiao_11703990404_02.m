
%ssbm()函数文件ssbm.m
function image1702(str)
l=0;r=1;d=1;
%初始间隔
about={...
        '本实例说明：=====>>'
    '字符串不能太长，程序不加判断，请注意溢出；'
    '本实例只限定少数字符串；'
    '实例只是说明一下算术编码过程。'};
disp(about);
%程序限定字符为：a、b、c、d、e
p=[0.2 0.3 0.1 0.15 0.25];
%字符的概率分布，sum(p)=1
n=length(str);
disp('a b c d e')
disp(num2str(p))
for i=1:n
    switch str(i)
        case 'a' 
            m=1;
        case 'b'
            m=2;
        case 'c'
            m=3;
        case 'd'
            m=4;
        case 'e'
            m=5;
        otherwise
            error('请不要输入其它字符！');
    end
    %判断字符
    pl=0;pr=0;
    for j=1:m-1 
        pl=pl+p(j);
    end
    for j=1:m 
        pr=pr+p(j);
    end
    %概率统计
    l=l+d*pl;
    r=l+d*(pr-pl);
    strl=strcat('输入第',int2str(i),'符号的间隔左边界：');    
    disp(strl);
    format long
    disp(l)
    strr=strcat('输入第',int2str(i),'符号的间隔右边界：');
    disp(strr);
    disp(r)
    d=r-l;
end

