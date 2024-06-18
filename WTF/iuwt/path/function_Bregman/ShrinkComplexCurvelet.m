function R=ShrinkComplexCurvelet(A,tau)

%========================================================
%
%  Execute the shrinkage of complex Curvelet coefficients
%  stored in array A with threshhold tau. The low
%  "frequencies" coefficients are keep without any
%  changes.
%  
%  Parameters:
%       A: input array of Framelet coeeficients
%       tau: threshhold
%       R: output array of threshholded coefficients
%
%   Author: Jerome Gilles
%   Institution: UCLA - Math Department
%   email: jegilles@math.ucla.edu
%=========================================================

L=length(A);

for l=1:L
   N=length(A{l});
   for nh=1:N
       if (l==1)
           R{l}{nh}=A{l}{nh};
       else
           R{l}{nh}=max(zeros(size(A{l}{nh})),abs(A{l}{nh})-tau).*exp(1i*angle(A{l}{nh}));
       end
   end
end