
clear
close all;
%定义HufData/Len为全局变量的结构体
global HufData;
global Len
disp('计算机正在准备输出哈夫曼编码结果,请耐心等待……');
%原始码字的灰度
a=imread('lena.bmp');
%分区画出原始图像和灰度直方图
figure;
subplot(1,2,1)
imshow(a);
%取消坐标轴和边框
axis off
box off
title('MATLAB自带图像','fontsize',13);
subplot(1,2,2);
axis off
box off
imhist(a);
title('图像灰度直方图','fontsize',13);
%图像的灰度统计
GrayStatistics=imhist(a);
GrayStatistics=GrayStatistics';
GrayRatioo=GrayStatistics/sum(GrayStatistics);
GrayRatioNO=find(GrayRatioo~=0);
Len=length(GrayRatioNO);
%初始化灰度集，防止系统随即赋予其垃圾值
GrayRatio=ones(1,Len);
for i=1:Len
GrayRatio(i)=GrayRatioo(i);
end
GrayRatio=abs(sort(-GrayRatio));
%将图像灰度概率赋予结构体
for i=1:Len
HufData(i).value=GrayRatio(i);
end
% 哈夫曼编码/霍夫曼编码
HuffmanCode(Len);
%输出码字
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
%计算计算机一共输出多少个哈夫曼编码/霍夫曼编码
zippedHuffman;
%计算在删去0灰度级压缩之前的原始图像字节容量
unzipped_delete=i*8;
%计算压缩比率
ratio_delete=zippedHuffman/unzipped_delete;
%计算图像的压缩比率
ad=num2str(ratio_delete*100);
str2=strcat(ad,'%');
disp(strcat('哈夫曼编码压缩比率','= ',str2))
%子程序：哈夫曼编码/霍夫曼编码函数HuffmanCode.m
function HuffmanCode(OriginSize)
global HufData;
global Len
for i=1:Len
%%霍夫曼编码树左边纪录为1
HufData(i).left=1;
%%霍夫曼编码树右边纪录为0
HufData(i).right=0;
%%输出码初始化为0
HufData(i).code=[];
%%排序列表初始化
SortList(i).symbol=i;
SortList(i).value=HufData(i).value;
end
%初始化原始消息数目
newsymbol=OriginSize;
for n=OriginSize:-1:2
%将N个消息进行排序
SortList=sortdata(SortList,n);
%将最后两个出现概率最小的消息合成一个消息
newsymbol=newsymbol+1;
HufData(newsymbol).value=SortList(n-1).value+SortList(n).value;
HufData(newsymbol).left=SortList(n-1).symbol;
HufData(newsymbol).right=SortList(n).symbol;
%将消息添加到列队的最后，为N-1个消息重新排序作好准备
SortList(n-1).symbol=newsymbol;
SortList(n-1).value=HufData(newsymbol).value;
end
%遍历霍夫曼树，获得霍夫曼编码/哈夫曼编码
visit(newsymbol,Len,[]);
end
%子程序：冒泡排序法函数sortdata.m
function reData=sortdata(SortList,n)
%根据消息概率进行排序
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
%子程序：遍历哈夫曼编码/霍夫曼编码树搜索函数visit.m
function visit(node,n,ocode)
    global HufData
    if node<=n
    %如果没有哈夫曼编码/霍夫曼编码树的子接点直接输出原始码，这里为空码（[]）
    HufData(node).code=ocode;
    else
    if(HufData(node).left>0)
    %遍历左分支接点输出1，这里采用子函数嵌套调用
    ocode1=[ocode 1];
    visit(HufData(node).left,n,ocode1);
    end
    if(HufData(node).right>0)
    %遍历右分支接点输出0，这里采用子函数嵌套调用
    ocode2=[ocode 0];
    visit(HufData(node).right,n,ocode2);
    end
    end
    end