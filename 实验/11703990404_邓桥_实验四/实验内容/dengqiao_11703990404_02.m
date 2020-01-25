%% ��ͨ�˲�����ͼ��
moon=imread('moon.tif');
[m,n]=size(moon);
I=im2double(moon);
figure;
subplot(1,2,1),imshow(moon),title('moon');
I_spectrum=fft2(I);
I_spectrum=fftshift(I_spectrum);
subplot(1,2,2),imshow(log(1+abs(I_spectrum)),[]),title('fourier spectrum of moon');

%% highpass filter
H=zeros(m,n);
centerx=m/2;
centery=n/2;
D0=60;  % �ɵ���ͨ���뾶������ͨ���ĸ�Ƶ����
for x=1:m
    for y=1:n
        H(x,y)=exp(-((x-centerx)^2+(y-centery)^2)/(2*D0^2));  %�����˹�˲�ģ��
    end
end
H=1-H;
figure,imshow(H,[]),title('gaussian filter');

g1=H.*I_spectrum;% ��Ƶ����ͼ�񣬱�Եͼ��
g2=g1+I_spectrum;
% �ȷ����Ļ�����ת������ͼ��
g3=ifftshift(g2);
I2=real(ifft2(g3));

figure,imshow(real(ifft2(ifftshift(g1))),[]),title('edge of moon');
figure;
subplot(1,2,1),imshow(log(1+abs(g2)),[]),title('fouier image sharpen by gaussian filter');
subplot(1,2,2),imshow(I2,[]),title('sharpen image');