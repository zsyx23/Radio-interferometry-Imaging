%% 一维高斯函数
x = -10:0.1:10;
% x=-350:1:350;
mu=0;
sigma=1;
y = 1/(sqrt(2*pi)*sigma)*exp(-(x-mu).^2/(2*sigma^2));
plot(x,y,'r');

%% 二维高斯函数    fspecial函数    h = fspecial('gaussian',hsize,sigma);
sigma1=2;  % 方差大小
k=700;    % kernel大小
kernel=zeros(k);
m=(k+1)/2;
sigma=2*sigma1*sigma1;
for i=-1*(k-1)/2:(k-1)/2
   for j=-1*(k-1)/2:(k-1)/2
      kernel(i+m,j+m)=(-1/(pi*sigma))*exp(-1*(i^2+j^2)/(sigma));
   end
end
kernel=kernel./sum(kernel,'all');  % 归一化
figure;pcolor(kernel);shading interp;colorbar;  %绘图

h = fspecial('gaussian',700,100);
figure;pcolor(h);shading interp;colorbar;  %绘图
figure; mesh(h);
max(max(h))

