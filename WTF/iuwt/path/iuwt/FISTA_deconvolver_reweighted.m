
function [Model Residual]=FISTA_deconvolver_reweighted(Dirtymap,PSF,lambda,niter,positiveflg,waveletflg,level)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA -- a L1 norm based algorithom for solving the deconvolution problem in
%radio astronomy
% any kind of UV coverage can be accepted.
% Details can be found in paper 
% F. Li, T. Cornwell and F. de Hoog, ?The Application of Compressive Sampling to Radio Astronomy I: Deconvolution?, Astronomy & Astrophysics, Volume 528, 2011"
% The only difference with the paper is that a reweighted version is used
% to improve the results.

%
% Dirtymap.......The blurred image
%
% PSF .......... Point Spread Function of the Dirymap,here the psf is surposed to
%                be the same size of Dirtymap
%
%
% lambda........ Regularization parameter,  in general, it should be smaller than
%                 the followed parameter 'threshold'
% niter......... The number of iterations
%
% positiveflg....The flag of the positive prior. If positiveflg=1, the source are
%                 all positive in the model, by default positiveflg=0;

% waveletflg ....The flag of Istropic Undecimated Wavelet Transform (IUWT). If
%                waveflg=1, IUWT is adopted, otherwise the partial Fourier
% level  ....... The level of the wavelet transform, for this case, the level is
%                no larger than 6. By default, it is set to  4


% Model .........The output cleaned image
% Residual ......The residual image which is: Dirtymap-Model*PSF;
%
wh=size(Dirtymap) == size(PSF);

if wh(1) & wh(2)
else
    error(' The dirtymap and the dirty beam have to be the same size');
end

if nargin < 7
    level=4;
end

if nargin < 6
    waveletflg=0;
    level=4;
end

if nargin < 5
    positiveflg=0;
    waveletflg=0;
    level=4;
end



[m,n]=find(PSF==max(max(PSF)));
center=[m n];


if waveletflg
    W = @(x) IUWT(x,level);% W means IUWT  С���任
    WT = @(x) IUWT(x,-level);% WT  means inverse IUWT
else
    W = @(x) x;   % ����Ҷ�任
    WT = @(x) x;
end
[m,n]=size(Dirtymap);


%Computng the UV mask with the psf
UV=fft2(circshift(PSF,1-center));
maxUV=max(max(abs(UV)));
indx=find(abs(UV/maxUV)>0.05);
invuv=zeros(m,n);
invuv(indx)=1./(UV(indx).^2);



%Initialization
im_temp=Dirtymap;% ��ʼͼ��Ϊ��ͼ
X_temp=W(zeros(m,n));% С���任���µ�ϵ��
weight=ones(size(X_temp));  %Ȩ��
for count=1:3
    Fdirty=fft2(Dirtymap);%����
    Residual=zeros(m,n);%�в�
    X=X_temp;
    t_new=1;
    for i=1:niter  % niterΪ��������
        X_old=X_temp;
        t_old=t_new;
        Residual_old=Residual;      
        im_temp=WT(X_temp); %С����任   X_tempΪϵ����im_tempΪʱ���µ�ͼ��    ��������Ľ��
        if positiveflg
                im_temp=im_temp.*(im_temp>0);
        end
        fim_temp=fft2(im_temp);% ��ǰͼ��ĸ���Ҷ�任  
        residual0=-Fdirty+UV.*fim_temp;  % residual0Ϊ�����Ĳв�   UVΪM����   M.*F
        D=conj(UV).*residual0.*invuv;  % DΪV�ɼ��Ⱥ���  conj���㹲��ֵ   D��500*700
        X=X-W(real(ifft2(D)));    %  ����С���任���µ�ϵ��   W(real(ifft2(D)))Ϊ�ɼ��Ⱥ�������Ҷ�任����ͼ����С���任ϵ��

        % Soft thresholding
        Shrink=abs(X)-lambda*weight;
        X_temp=sign(X).*((Shrink>0).*Shrink); % X_tempΪ����������   �����еĹ�ʽ��23��

        %Updating t and X
        t_new=(1+sqrt(1+4*t_old^2))/2; % t_newΪt��k+1��   �����еĹ�ʽ��24��
        X=X_temp+(t_old-1)/t_new*(X_temp-X_old);  % XΪ�£�k+1��������IUWT������ͼ��С��ϵ�� ά��Ϊ500*3500     �����еĹ�ʽ��25��
               
      
        %Evaluating
        Residual=Dirtymap-real(ifft2(UV.*fft2(im_temp)));%��ͼ�Ĳв�
        likelyhood=norm(Residual,'fro')^2;
        total=likelyhood+sum(sum(abs(W(im_temp))));
        maxresidual=max(max(abs(Residual)));
        
        %Print the progress on the screen
%         fprintf('%3d    %15.5f   %15.5f   %15.5f \n',i,likelyhood,total,maxresidual);

        changes=norm(Residual-Residual_old,'fro')^2/likelyhood;
        if i~=1
            if changes<1e-7 
                break;
            end
        end
    end
    weight=1./(abs(W(im_temp))+.001);  % Ȩ��
    %Here we assume 0.001 is the minimum value in the original signal
    %Bring the current value as the weight for the next external iterations

end
Model=im_temp;