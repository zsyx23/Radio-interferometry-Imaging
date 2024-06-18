function extended_test()
%This test is conducted for the extended source in the datase
%"exampleimages", in which the Dirtymap, PSF, Original images are included
load exampleimages



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz = get(0,'ScreenSize');
figure(1)
set(1,'Position',[20 scrsz(4)*0.05 scrsz(3)/4 1*scrsz(3)])
subplot(3,1,1)
imshow(Original,[0,2])
title('Original Model','FontName','Times','FontSize',14)
subplot(3,1,2)
imshow(PSF,[0,1])
title('Point Spread Function','FontName','Times','FontSize',14)
subplot(3,1,3)
imshow(Dirtymap,[-2,2])
title('Dirtymap','FontName','Times','FontSize',14)
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


%Test for the PF algorithm
fprintf('Starting PF based FISTA deconvolver');
[PF_Model PF_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.01,100,1);


%Or Test for the IUWT algorithm
fprintf('Starting IUWT based FISTA deconvolver');
[IUWT_Model IUWT_Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,0.0001,100,1,1,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot figures
scrsz = get(0,'ScreenSize');
figure(2)
set(2,'Position',[20 scrsz(4)*0.05 scrsz(3)/4 1*scrsz(4)])

subplot(4,1,1)
imshow(PF_Model,[0,2]);
title('PF Model','FontName','Times','FontSize',14)

subplot(4,1,2)
imshow(PF_Residual,[-0.02,0.02]);
title('PF residual','FontName','Times','FontSize',14)

subplot(4,1,3)
imshow(IUWT_Model,[0,2]);
title('IUWT Model','FontName','Times','FontSize',14)

subplot(4,1,4)
imshow(IUWT_Residual,[-0.02,0.02]);
title('IUWT residual','FontName','Times','FontSize',14)
axis off