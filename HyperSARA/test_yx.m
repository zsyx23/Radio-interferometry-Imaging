
close all;
clear;


% % x0 原始图像
% [~,~,len]=size(x0);
% for iii=1:len
%     figure;
%     pcolor(real(x0(:,:,iii)));
%     shading interp;axis xy;colorbar;colormap(jet);  
% end

% 可见度函数
% [len2,~]=size(X0);
% for iii=1:len
%     if iii<1
%     figure; plot(1:len2,X0(:,iii));hold on;
%     else
%         plot(1:len2,X0(:,iii));hold on;
%     end
% end


% % xsol 重构图像
% [~,~,len]=size(xsol);
% for iii=1:len
%     figure;
%     pcolor(real(xsol(:,:,iii)));
%     shading interp;axis xy;colorbar;colormap(jet);  
% end


% x0和xsol重构图像对比
[~,~,len]=size(xsol);
SNR_sum=[];
for iii=1:len
%     figure;
%     pcolor(real(x0(:,:,iii)));
%     shading interp;axis xy;colorbar;colormap(jet);  
         figure;
    pcolor(real(log10(x0(:,:,iii))));title('x0');
    shading interp;axis xy;colorbar;colormap(jet); caxis([-3.5,1]);
%         figure;
%     pcolor(real(xsol(:,:,iii)));
%     shading interp;axis xy;colorbar;colormap(jet); 
     figure;
    pcolor(real(log10(xsol(:,:,iii))));title('xsol');
    shading interp;axis xy;colorbar;colormap(jet); caxis([-3.5,1]);
    
    
    %信噪比
    SNR=10*log10((sum(sum(x0(:,:,iii).^2)))/(sum(sum((xsol(:,:,iii)-x0(:,:,iii)).^2))));
    SNR_sum=[SNR_sum,SNR];

end

SNR_ave=sum(SNR_sum)./length(SNR_sum);

% figure;
% plot(u1,v1,'.r');
% figure;
% plot(u{1,1},v{1,1},'.r');







