function u=Curvelet_NB_Deconvolution(f,A,mu,lambda,Niter,NLevel,typeKernel)

%============================================================================
%
% function u=Curvelet_NB_Deconvolution(f,A,mu,lambda,Niter,NLevel,typeKernel)
%
% This function performs the Nonblind Curvelet deconvolution by the analysis 
% approach.
% Version:
% -v1.4 - 06/20/2013
% -v1.2 - 07/18/2011
% -v1.0 - 05/04/2011
%
% Parameters:
% f: input blurred image
% A: input blur kernel
% mu,lambda: regularization parameters
% Niter: number of iterations
% NLevel: number of scale used in the Curvelet decomposition
% typeKernel: 0 <=> A is considered in the spatial domain
%             1 <=> A is considered in the Fourier domain and must be of
%             the same size of the image with frequency (0,0) at location
%             (1,1)
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%============================================================================


%structures initialization
[M,N]=size(f);
U=fdct_wrapping(zeros(M,N),1,1,NLevel);
d=U;
b=U;

%Fourier mask + initialization of Fourier constant quantities
if typeKernel==0
    Mask=zeros(M,N);
    [H,L]=size(A);
    Mask([end+1-floor(H/2):end,1:ceil(H/2)],[end+1-floor(L/2):end,1:ceil(L/2)]) = A;
    FMask=fft2(Mask);
else
    FMask = A;
end

FW=(mu*abs(FMask).^2+lambda).^-1;
FF=mu*conj(FMask).*fft2(f);

K=0;
err=norm(f(:),2);
tol=1e-3*err;
u=zeros(size(f));
while ((err>tol) && (K<Niter)),
    K=K+1;
    
    up=u;
    %update u
    tx=ifdct_wrapping(SubCurveletArray(d,b),1,M,N);
    u=real(ifft2(FW.*(FF+lambda*fft2(tx))));
    
    %update d
    U=AddCurveletArray(fdct_wrapping(u,1,1,NLevel),b);
    d=ShrinkCurvelet(U,1/lambda);
    
    %update b
    b=SubCurveletArray(U,d);
    
    err=sum(sum((up-u).^2));
end
