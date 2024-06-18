close all;
clear;

airi_model_image=fitsread('airi_model_image.fits','raw');
airi_residual_dirty_image=fitsread('airi_residual_dirty_image.fits','raw');
airi_residual_dirty_image_normalised=fitsread('airi_residual_dirty_image_normalised.fits','raw');
dirty=fitsread('dirty.fits','raw');
PSF=fitsread('PSF.fits','raw');


%% 结果
figure; imagesc(real(log10((airi_model_image))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-3 0]);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);


figure; imagesc(real(((airi_model_image))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);


%% 残差
figure; imagesc(real(((airi_residual_dirty_image))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);

figure; imagesc(real((log10(airi_residual_dirty_image))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);


%% 脏图
figure; imagesc(real((log10(dirty))));set(gca,'YDir','normal');axis off;
colorbar, axis image; % caxis([-5.6 -3.6]);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);

figure; imagesc(real((log10(dirty))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);

figure; imagesc(real(((dirty))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%% 点扩散函数
figure; imagesc(real((log10(PSF))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);

figure; imagesc(real(((PSF))));set(gca,'YDir','normal');axis off;
colorbar, axis image; 
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
colormap(cubehelix);






