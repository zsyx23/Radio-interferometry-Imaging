addpath('C:\Users\DELL\Desktop\L1\function')

load exampleimages

% [u,mu,K,beta]=L1Adpt_tiduWWSplitBregmanIteration(f,A,mu,beta,Niter,a)

[m n]=find(PSF==max(max(PSF)));
center=[m n];
[m,n]=size(Dirtymap);
%Computng the UV mask with the psf
UV=fft2(circshift(PSF,1-center));
maxUV=max(max(abs(UV)));
indx=find(abs(UV/maxUV)>0.05);
invuv=zeros(m,n);
invuv(indx)=1./(UV(indx).^2);

% W = @(x) x;   % 傅里叶变换
% WT = @(x) x;
level=4;
W = @(x) IUWT(x,level);% W means IUWT  小波变换
WT = @(x) IUWT(x,-level);% WT  means inverse IUWT

%Initialization
im_temp=Dirtymap;% 初始脏图
X_temp=W(zeros(m,n));% 脏图的小波变换域下的系数
weight=ones(size(X_temp));

Fdirty=fft2(Dirtymap);%脏束
Residual=zeros(m,n);%残差
X=X_temp;
t_new=1;

X_old=X_temp;
t_old=t_new;
Residual_old=Residual;
im_temp=WT(X_temp); %小波逆变换   im_temp为时域下的脏图    最终输出的结果
positiveflg=1;
if positiveflg
    im_temp=im_temp.*(im_temp>0);
end
fim_temp=fft2(im_temp);%脏图的傅里叶变换
residual0=-Fdirty+UV.*fim_temp;  % residual0为脏束的残差   UV为M矩阵   M.*F
D=conj(UV).*residual0.*invuv;  % D为V可见度函数  conj计算共轭值
X=X-W(real(ifft2(D)));    %  更新小波变换域下的系数   W(real(ifft2(D)))为可见度函数傅里叶变换（脏图）的小波变换系数




%L1 sparse recovery
% u = arg min |u|+0.5*mu||Au-f||_2^2
A=W(real(ifft2(D)));




N=size(f,1)
d=zeros(N,1);
b=zeros(N,1);
u=zeros(N,1);
Z=zeros(N,1);
Ft=mu*A'*f;
IV=inv(mu*A'*A+lambda*eye(N));
up=ones(N,1);
while ((u-up)'*(u-up))>eps,
up=u; %store the previous iteration
u=IV*(Ft+lambda*(d-b)); %update u
tmp=u+b;
d=sign(tmp).*max(Z,abs(tmp)-1/lambda); %update d
b=tmp-d; %update b
end






