%=========================================================
%
% APPLICATION EXAMPLES OF THE BREGMAN COOKBOOK FOR 3D DATA
%
% -v 1.0: 03/01/2012
% 
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
% 
%=========================================================
%clear;

%Choose the experiment you want to perform
expe=1; % 0=denoising ; 1=deblurring

%load the MRI data
load mri
D = double(squeeze(D))/256;

if expe==0
%======================
%
% DENOISING EXPERIMENTS
%
%======================

%build a noisy 3D datacube
Noise=rand(size(D));
ND=D+0.1*Noise;

%process the denoising
mu=500;
lambda=100;
Niter=10;

%Choose which method you want to apply
DND=ATV_ROF_3D(ND,mu,lambda,Niter);
%DND=ITV_ROF_3D(ND,mu,lambda,Niter);

subplot(1,2,1);imshow(ND(:,:,10));subplot(1,2,2);imshow(DND(:,:,10));

elseif expe==1
%=======================
%
% DEBLURRING EXPERIMENTS
%
%=======================
Dp=D;

%build a gaussian blurred 3D datacube
gauss3d=fspecial3('gaussian');
BD=convn(Dp,gauss3d,'same');

%process the deblurring
Debmu=1000000;
Deblambda=30;
DebNiter=10;

%Choose which method you want to apply
%DBD=ATV_NB_Deconvolution_3D(BD,gauss3d,Debmu,Deblambda,DebNiter,0);
%DBD=ITV_NB_Deconvolution_3D(BD,gauss3d,Debmu,Deblambda,DebNiter,0);
DBD=Curvelet_NB_Deconvolution_3D(BD,gauss3d,Debmu,Deblambda,DebNiter,0);

slice=11; %Choose the slice of the cube you want show
subplot(1,3,1);imshow(Dp(:,:,slice));subplot(1,3,2);imshow(BD(:,:,slice));subplot(1,3,3);imshow(DBD(:,:,slice));

end

