%========================================================
%
% APPLICATION EXAMPLES OF THE BREGMAN COOKBOOK FOR IMAGES
%
% -v 1.2: 06/20/2013
% -v 1.0: 03/01/2012
% 
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
% 
%========================================================
clear;

%Choose which experiment you want to perform
expe=1; % 0 = denoising ; 1 = deblurring

%load the image
load lena;

if expe==0
%======================
%
% DENOISING EXPERIMENTS
%
%======================

%build a noisy image
Noise=rand(size(lena));
ND=lena+0.1*Noise;

%process the denoising
mu=100;
lambda = 50;
Niter = 10;

%Choose which method you want to apply
%DND=ATV_ROF(ND,mu,lambda,Niter);
DND=ITV_ROF(ND,mu,lambda,Niter);

subplot(1,2,1);imshow(ND,[]);xlabel('Noisy');subplot(1,2,2);imshow(DND,[]);xlabel('Denoised')

elseif expe==1
%=======================
%
% DEBLURRING EXPERIMENTS
%
%=======================

%we extend the image by mirroring to get rid of boundaries artifacts
W0=size(lena,2);
H0=size(lena,1);
W1=20;
H1=20;
Dp=[lena(H1:-1:1,W1:-1:1) lena(H1:-1:1,:) lena(H1:-1:1,end:-1:end-W1+1) ; lena(:,W1:-1:1) lena lena(:,end:-1:end-W1+1) ; lena(end:-1:end-H1+1,W1:-1:1) lena(end:-1:end-H1+1,:) lena(end:-1:end-H1+1,end:-1:end-W1+1)];

%build a gaussian kernel
gauss=fspecial('gaussian',9,1.5);
BD=conv2(Dp,gauss,'same');

%process the deblurring
Debmu=100000;
Deblambda=50;
DebNiter=20;

%Choose which method you want to apply
%DBD=ATV_NB_Deconvolution(BD,gauss,Debmu,Deblambda,DebNiter,0);
%DBD=ITV_NB_Deconvolution(BD,gauss,Debmu,Deblambda,DebNiter,0);
%DBD=Framelet_NB_Deconvolution(BD,gauss,Debmu,Deblambda,DebNiter,0,3,0);
%DBD=Framelet_NB_Deconvolution2(BD,gauss,Debmu,Deblambda,1,DebNiter,0,3,0);
DBD=Curvelet_NB_Deconvolution(BD,gauss,Debmu,Deblambda,DebNiter,3,0);

DBD=DBD(H1+1:H1+H0,W1+1:W1+W0);
BD=BD(H1+1:H1+H0,W1+1:W1+W0);
subplot(1,3,1);imshow(lena,[]);xlabel('Original');subplot(1,3,2);imshow(BD,[]);xlabel('Blurry');subplot(1,3,3);imshow(abs(DBD),[]);xlabel('Deblurred');

end