clc
clear
close all

% 多尺度 CLEAN 函数的卷积核函数
%% 球面波函数
lamda=6328e-10;                      %波长，单位：米
k=2*pi/lamda;                        %波数
x0=0.0001;                             %点光源的x坐标，单位：米
y0=0.0001;                            %点光源的y坐标，单位：米
z=0.3;                                %观察面到点光源的垂直距离，单位：米
L=0.005;                              %观察面的尺寸，单位：米
x=linspace(-L/2,L/2,700);y=x;        %构建x坐标和y坐标
[x,y]=meshgrid(x,y);                 %构建二维坐标网格
U1=exp(j*k*z).*exp(j*k.*((x-x0).^2+(y-y0).^2)/2/z);  %发散球面光波
figure;pcolor(real(U1));shading interp;colorbar;
Fai=real(U1);%取实数部分    
%% 抛物线函数
% r=sqrt(x^2+y^2);
% m=Fai.*(1-(r/a)^2);%抛物线函数

%x,y为矩阵位置坐标
a = 24 ;%a为当前尺度大小
m=zeros(700,700);
for i=1:700
    for j=1:700
        r=sqrt(i^2+j^2);
        m(i,j)=Fai(i,j).*(1-(r/a)^2)/700;%抛物线函数
%         if m(i,j)<0
%             m(i,j)=0;
%         end
    end
end
figure;pcolor(real(m));shading interp;colorbar;
m_max=max(max(real(m)))

%x,y为像素点位置坐标
% a = 2 ;%a为当前尺度大小
% % m=real(exp(j*k*z).*exp(j*k.*((x-x0).^2+(y-y0).^2)/2/z)).*(1-(sqrt(x^2+y^2)./a)^2);
% Fai=real(U1);%取实数部分   
% m=zeros(700,700);
% for i=1:700
%     for j=1:700
%         r=sqrt(x(i,j)^2+y(i,j)^2);
%         m(i,j)=Fai(i,j)*(1-(r/a)^2);%抛物线函数
%     end
% end
% figure;pcolor(real(m));shading interp;colorbar;
% max(max(real(m)))

