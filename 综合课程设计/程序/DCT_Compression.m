%DCTѹ��--------------------------------------------------------------------
function DCT_Compression(handles)
% ����ѹ���ȣ�cr=0.5Ϊ2:1ѹ����cr=0.1250Ϊ8:1ѹ��
% fprintf('������ѹ����Ϊ��1��0.5��0.25��0.125\n')
global im;
mysize=size(im);
if numel(mysize)>2
    I1=rgb2gray(im);
else
    I1=im;
end
data1=get(handles.edit9,'string');
data1_num=str2num(data1);
cr=data1_num;
initialimage = I1;
initialimage = double(initialimage)/255;
[mm,nn]=size(initialimage);
figure;
subplot(121);
imshow(initialimage);
title('ԭʼͼ��');
disp('ԭͼ���С��')
whos('initialimage');
tic;
%��ͼ�����DCT�任
t = dctmtx(8);
dctcoe = blkproc(initialimage, [8 8], 'P1*x*P2', t, t');
%��DCT�任��ľ���ת�����У�������������
coevar = im2col(dctcoe, [8 8], 'distinct');
coe = coevar;
[y, ind] = sort(coevar);
[m, n] = size(coevar);
%��ȥ����Ҫ��ϵ��
snum = 64-64 * cr;
for i = 1:n
    coe(ind(1:snum), i) = 0;
end
%���б任Ϊ��ά����
b2 = col2im(coe, [8 8], [mm nn], 'distinct');
t1=toc;
xl1=sprintf('ͼ�������ʱ%4.2f��\nͼ��ѹ����Ϊ��%4.1f��1\n',t1,1/cr);
xlabel(xl1);

%��DCT�任

tic;
i2 = blkproc(b2, [8 8], 'P1*x*P2', t', t);
t2=toc;

PSNR = psnr(initialimage, i2);
xl2=sprintf('ͼ�������ʱ%4.2f��\n����ͼ��ķ�ֵ�����Ϊ��%4.2f dB\n\n',t2,PSNR);
subplot(122);
imshow(i2);
title('����ͼ��');xlabel(xl2);
disp('����ͼ���С��')
whos('i2');
end

