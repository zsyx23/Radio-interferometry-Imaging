function [thr_coef,TTT ]= thr_bishrink222(coef)

%set the window size and the corresponding filter
% windowsize = 7;
windowsize = 7 ;
windowfilt = ones(1,windowsize)/windowsize;

n_Lvl = length(coef)-1;

level_length=length(coef{1});
sum=0;
TTT=[];

for i=1:level_length
sum=sum+median(abs(coef{1}{i}(:)));
end
sum=sum/level_length;
% Nsig=sum/0.6745;
Nsig=sum/0.6745;
thr_coef = coef;
for scale = 1:n_Lvl-1
   L = length(coef{scale});
   for l  = 1:L
            Y_coef = coef{scale}{l};
            Y_parent = coef{scale+1}{l};
            Y_parent = expandx(Y_parent);
            Wsig = conv2(windowfilt,windowfilt,(abs(Y_coef)).^2,'same');     
            
            Ssig = sqrt(max(Wsig - Nsig.^2, eps));
%             Ssig = sqrt(max(Wsig - Nsig.^2, 1e-50));
%             Ssig = sqrt((Wsig - Nsig.^2));
%             Ssig = sqrt(max(Wsig - Nsig.^2, (eps/1e-5)));
%             Ssig = sqrt(max(Wsig - Nsig.^2, (eps/1e-3)));
            
%             T = sqrt(3)*Nsig^2./Ssig; 
%             T = ones(size(real((sqrt(3)*Nsig^2./Ssig))))*1e-6 ;
%             T = T *(1e2) ; 
%             T = sqrt(3)*Nsig^(1.5)./Ssig; % 1.5

            T = sqrt(3)*Nsig^(2)./Ssig;
% T=10;
% TTT=[TTT,T];

% figure; imagesc(real(((T))));%set(gca,'YDir','normal');
% axis off;
% colorbar, axis image; 
% % caxis([-5.6 -3.6]);
% % caxis([-5.7 -3.5]);
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

            Y_coef = bishrink(Y_coef, Y_parent, T);
            thr_coef{scale}{l} = Y_coef;
   end
end
    

end