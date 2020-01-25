
clear
close all;
%����HufData/LenΪȫ�ֱ����Ľṹ��
global HufData;
global Len
disp('���������׼�����������������,�����ĵȴ�����');
%ԭʼ���ֵĻҶ�
a=imread('lena.bmp');
%��������ԭʼͼ��ͻҶ�ֱ��ͼ
figure;
subplot(1,2,1)
imshow(a);
%ȡ��������ͱ߿�
axis off
box off
title('MATLAB�Դ�ͼ��','fontsize',13);
subplot(1,2,2);
axis off
box off
imhist(a);
title('ͼ��Ҷ�ֱ��ͼ','fontsize',13);
%ͼ��ĻҶ�ͳ��
GrayStatistics=imhist(a);
GrayStatistics=GrayStatistics';
GrayRatioo=GrayStatistics/sum(GrayStatistics);
GrayRatioNO=find(GrayRatioo~=0);
Len=length(GrayRatioNO);
%��ʼ���Ҷȼ�����ֹϵͳ�漴����������ֵ
GrayRatio=ones(1,Len);
for i=1:Len
GrayRatio(i)=GrayRatioo(i);
end
GrayRatio=abs(sort(-GrayRatio));
%��ͼ��Ҷȸ��ʸ���ṹ��
for i=1:Len
HufData(i).value=GrayRatio(i);
end
% ����������/����������
HuffmanCode(Len);
%�������
zippedHuffman=1;
for i=1:Len
tmpData=HufData(i).code;
str='';
for j=1:length(tmpData)
str=strcat(str,num2str(tmpData(j)));
zippedHuffman=zippedHuffman+1;
end
disp(strcat('a',num2str(i),'= ',str))
end
i;
%��������һ��������ٸ�����������/����������
zippedHuffman;
%������ɾȥ0�Ҷȼ�ѹ��֮ǰ��ԭʼͼ���ֽ�����
unzipped_delete=i*8;
%����ѹ������
ratio_delete=zippedHuffman/unzipped_delete;
%����ͼ���ѹ������
ad=num2str(ratio_delete*100);
str2=strcat(ad,'%');
disp(strcat('����������ѹ������','= ',str2))
%�ӳ��򣺹���������/���������뺯��HuffmanCode.m
function HuffmanCode(OriginSize)
global HufData;
global Len
for i=1:Len
%%��������������߼�¼Ϊ1
HufData(i).left=1;
%%�������������ұ߼�¼Ϊ0
HufData(i).right=0;
%%������ʼ��Ϊ0
HufData(i).code=[];
%%�����б��ʼ��
SortList(i).symbol=i;
SortList(i).value=HufData(i).value;
end
%��ʼ��ԭʼ��Ϣ��Ŀ
newsymbol=OriginSize;
for n=OriginSize:-1:2
%��N����Ϣ��������
SortList=sortdata(SortList,n);
%������������ָ�����С����Ϣ�ϳ�һ����Ϣ
newsymbol=newsymbol+1;
HufData(newsymbol).value=SortList(n-1).value+SortList(n).value;
HufData(newsymbol).left=SortList(n-1).symbol;
HufData(newsymbol).right=SortList(n).symbol;
%����Ϣ��ӵ��жӵ����ΪN-1����Ϣ������������׼��
SortList(n-1).symbol=newsymbol;
SortList(n-1).value=HufData(newsymbol).value;
end
%����������������û���������/����������
visit(newsymbol,Len,[]);
end
%�ӳ���ð�����򷨺���sortdata.m
function reData=sortdata(SortList,n)
%������Ϣ���ʽ�������
for k=n:-1:2
for j=1:k-1
min=SortList(j).value;
sbl=SortList(j).symbol;
if(min<SortList(j+1).value)
SortList(j).value=SortList(j+1).value;
SortList(j+1).value=min;
SortList(j).symbol=SortList(j+1).symbol;
SortList(j+1).symbol=sbl;
end
end
end
reData=SortList;
end
%�ӳ��򣺱�������������/��������������������visit.m
function visit(node,n,ocode)
    global HufData
    if node<=n
    %���û�й���������/���������������ӽӵ�ֱ�����ԭʼ�룬����Ϊ���루[]��
    HufData(node).code=ocode;
    else
    if(HufData(node).left>0)
    %�������֧�ӵ����1����������Ӻ���Ƕ�׵���
    ocode1=[ocode 1];
    visit(HufData(node).left,n,ocode1);
    end
    if(HufData(node).right>0)
    %�����ҷ�֧�ӵ����0����������Ӻ���Ƕ�׵���
    ocode2=[ocode 0];
    visit(HufData(node).right,n,ocode2);
    end
    end
    end