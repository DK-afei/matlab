
%ssbm()�����ļ�ssbm.m
function image1702(str)
l=0;r=1;d=1;
%��ʼ���
about={...
        '��ʵ��˵����=====>>'
    '�ַ�������̫�������򲻼��жϣ���ע�������'
    '��ʵ��ֻ�޶������ַ�����'
    'ʵ��ֻ��˵��һ������������̡�'};
disp(about);
%�����޶��ַ�Ϊ��a��b��c��d��e
p=[0.2 0.3 0.1 0.15 0.25];
%�ַ��ĸ��ʷֲ���sum(p)=1
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
            error('�벻Ҫ���������ַ���');
    end
    %�ж��ַ�
    pl=0;pr=0;
    for j=1:m-1 
        pl=pl+p(j);
    end
    for j=1:m 
        pr=pr+p(j);
    end
    %����ͳ��
    l=l+d*pl;
    r=l+d*(pr-pl);
    strl=strcat('�����',int2str(i),'���ŵļ����߽磺');    
    disp(strl);
    format long
    disp(l)
    strr=strcat('�����',int2str(i),'���ŵļ���ұ߽磺');
    disp(strr);
    disp(r)
    d=r-l;
end

