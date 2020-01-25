function cr = imageratio(f1, f2)
%计算两幅图像压缩比
% error(nargchk(2, 2, nargin));
cr = bytes(f1)/bytes(f2);
disp('压缩前图像数据量（比特数）:');
bytes(f1)
disp('压缩前图像数据量（比特数）:');
bytes(f2)
