%显示一副灰度图像
I = imread('moon.tif');                        
%读入图像数据
imagesc(I,[0 255]);                          
%预处理
colormap(gray);                            
%灰度处理 显示灰度图像
pause                                    
%暂停
imshow(I);                                
%显示灰度图像
