
clear;

load('3c353_meas_dt_1_seed_0.mat')

uvw(:,1)=u;
uvw(:,2)=v;
uvw(:,3)=w;
weights=ones(size(uvw(:,1)));
y_I=y;
y_V=double(zeros(size(y_I)));

save('3c353_data.mat','uvw','weights','y_I','y_V');

