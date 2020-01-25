%霍夫曼编码压缩-----------------------------------------------------------------
function Huffen_Compression
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
X = I1;
[mm,nn]=size(X);

tic;
vec=reshape(X,1,mm*nn);
[zipped, info] = huffencode(vec);
t1=toc;
cr = imageratio(X, zipped);
xl1=sprintf('图像编码用时%4.2f秒\n图像压缩比为：%4.1f：1\n',t1,cr);

tic;
af_vec = huffdecode(zipped, info);
XX=reshape(af_vec,mm,nn);
t2=toc;

PSNR = psnr(X, XX);
xl2=sprintf('图像解码用时%4.2f秒\n两幅图像的峰值信噪比为：%4.2f dB\n\n',t2,PSNR);

figure;subplot(121);imshow(X);
title('原始图像');xlabel(xl1);
disp('原图像大小：')
whos('X');
subplot(122);imshow(uint8(XX));
title('解码图像');xlabel(xl2);
disp('解码图像大小：')
whos('XX');
end

%%编码过程%%%%
function [zipped, info] = huffencode(vector)
% 输入和输出都是 uint8 格式
% info 返回解码需要的结构信息
% info.pad 是添加的比特数
% info.huffcodes 是 Huffman 码字
% info.rows 是原始图像行数
% info.cols 是原始图像列数
% info.length 是原始图像数据长度
% info.maxcodelen 是最大码长

if ~isa(vector, 'uint8')
    error('input argument must be a uint8 vector');
end

[m, n] = size(vector);
vector = vector(:)';
f = frequency(vector);      %计算各符号出现的概率
symbols = find(f~=0);
f = f(symbols);
[f, sortindex] = sort(f);    %将符号按照出现的概率大小排列
symbols = symbols(sortindex);
len = length(symbols);
symbols_index = num2cell(1:len);
codeword_tmp = cell(len, 1);

% 生成 Huffman 树，得到码字编码表
while length(f)>1
    index1 = symbols_index{1};
    index2 = symbols_index{2};
    codeword_tmp(index1) = addnode(codeword_tmp(index1), uint8(0));
    codeword_tmp(index2) = addnode(codeword_tmp(index2), uint8(1));
    f = [sum(f(1:2)),f(3:end)];
    symbols_index = [{[index1, index2]},symbols_index(3:end)];
    [f, sortindex] = sort(f);
    symbols_index = symbols_index(sortindex);
end
codeword = cell(256, 1);
codeword(symbols) = codeword_tmp;
len = 0;
for index = 1:length(vector)       %得到整个图像所有比特数
    len = len + length(codeword{double(vector(index))+1});
end
string = repmat(uint8(0), 1, len);
pointer = 1;
for index = 1:length(vector)       %对输入图像进行编码
    code = codeword{double(vector(index))+1};
    len = length(code);
    string(pointer + (0:len-1))=code;
    pointer = pointer + len;
end
len = length(string);
pad = 8-mod(len, 8);
if pad > 0
    string = [string uint8(zeros(1, pad))];
end
codeword = codeword(symbols);
codelen = zeros(size(codeword));
weights = 2.^(0:23);
maxcodelen = 0;
for index = 1:length(codeword)
    len = length(codeword{index});
    if len > maxcodelen;
        maxcodelen = len;
    end
    if len > 0
        code = sum(weights(codeword{index} == 1));
        code = bitset(code, len + 1);
        codeword{index} = code;
        codelen(index) = len;
    end
end
codeword = [codeword{:}];
    
%计算压缩的向量
cols = length(string)/8;
string = reshape(string, 8, cols);
weights = 2.^(0: 7);
zipped = uint8(weights * double(string));
    
%码表存储到一个希疏矩阵
huffcodes = sparse(1, 1);
for index = 1:nnz(codeword)   % length(codeword)  %numel(codeword)
    huffcodes(codeword(index), 1) = symbols(index);
end
    
%填写解码时所需的结构信息
info.pad = pad;
info.huffcodes = huffcodes;
info.ratio = cols./length(vector);
info.length = length(vector);
info.maxcodelen = maxcodelen;
info.rows = m;
info.cols = n;
end
   
%函数addnode添加节点
function codeword_new = addnode(codeword_old, item)
codeword_new = cell(size(codeword_old));
for index = 1:length(codeword_old)
    codeword_new{index} = [item codeword_old{index}];
end
end

%函数frequency计算各符号出现的概率
function f = frequency(vector)
if ~isa(vector, 'uint8')
    error('input argument must be a uint8 vector');
end
f = repmat(0, 1, 256);
len = length(vector);
for index = 0:255
    f(index+1) = sum(vector == uint8(index));
end
f = f./len;   %归一化
end

%%解码过程%%%%%%%%%%%%
function vector = huffdecode(zipped, info)
% 函数对输入矩阵vector进行Huffman解码，返回解压后的图像数据

if ~isa(zipped, 'uint8')
    error('input argument must be be a uint8 vector');
end

%产生0，1序列，每位占一个字节
len = length(zipped);
string = repmat(uint8(0), 1, len.*8);
bitindex = 1:8;
for index = 1:len
    string(bitindex + 8.*(index-1)) = uint8(bitget(zipped(index), bitindex));
end
string = logical(string(:)');
len = length(string);
string ((len-info.pad+1):end)=[];
len = length(string);

%开始解码
weights = 2.^(0:51);
vector = repmat(uint8(0), 1, info.length);
vectorindex = 1;
codeindex = 1;
code = 0;
for index = 1:len
    code = bitset(code, codeindex, string(index));
    codeindex = codeindex+1;
    byte = decode(bitset(code, codeindex), info);
    if byte > 0
        vector(vectorindex) = byte-1;
        codeindex = 1;
        code = 0;
        vectorindex = vectorindex + 1;
    end
end
vector = reshape(vector, info.rows, info.cols);
end

%函数decode返回码字对应的符号
function byte = decode(code, info)
byte = info.huffcodes(code);
end
