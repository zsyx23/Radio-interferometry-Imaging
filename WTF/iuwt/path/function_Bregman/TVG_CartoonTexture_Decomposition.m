function [u,v]=TVG_CartoonTexture_Decomposition(f,mu,lambda,Breglambda)
%=====================================================================
% function u=TVG_CartoonTexture_Decomposition(f,mu,lambda,Breglambda)
%
% TV-G Cartoon + Texture decomposition
% Version:
% -v1.0 - 10/30/2011
%
% This function performs Meyer's Cartoon + Texture decomposition of an
% image based on Aujol's formulation by Split Bregman Iterations
%
% Reference: J.Gilles, S.Osher, "Bregman implementation of Meyer's 
%            G-norm for cartoon + textures decomposition" available as 
%            UCLA-CAM Report CAM11-73
%
% f = input image (must be normalized between 0 and 1)
% mu = texture regularization parameter
% lambda = cartoon regularization 
% Breglambda =  regularization parameter for the Bregman variable
%               (a value of 50 works for most of images)
%
% The return quantities are:
% u = the cartoon part
% v = the texture part
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%=====================================================================
u=zeros(size(f));
v=zeros(size(f));

err=norm(f(:),2);
tol=1e-3*norm(f(:),2);
Niter=10;
while err>tol,
    up=u;
    vp=v;
    
    %update u
    tmp=f-v;
    u=ITV_ROF(tmp,lambda,Breglambda,Niter);

    %update v
    tmp=f-u;
    v=ITV_ROF(tmp,1/mu,Breglambda,Niter);
    v=tmp-v;
    
    err=max(sum(sum((u-up).^2)),sum(sum((v-vp).^2)));
end
