%图像文件的读取、写盘
[x,map] = imread('canoe.tif');
imwrite(x,map, 'canoe.tif', 'Compression', 'none', 'WriteMode', 'append');
