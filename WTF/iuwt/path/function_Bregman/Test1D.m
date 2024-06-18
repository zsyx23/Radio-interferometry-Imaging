%=========================================================
%
% APPLICATION EXAMPLES OF THE BREGMAN COOKBOOK FOR SIGNALS
%
% -v 1.0: 06/20/2013
% 
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
% 
%=========================================================
clear;

%Build a sparse signal
f=zeros(500,1);
f(55)=0.3;
f(132)=0.7;
f(178)=0.6;
f(221)=0.4;
f(269)=0.85;
f(341)=0.76;
f(429)=0.51;

subplot(311);plot(f);

%Apply a random linear operator
A=rand(500,500);
f=A*f;
subplot(312);plot(f);

%L1 minimization
mu=10;
lambda=1;
Niter=10;

u=L1_SplitBregmanIteration(f,A,mu,lambda,Niter);

subplot(313);plot(u);