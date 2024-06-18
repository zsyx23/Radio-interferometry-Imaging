%% 来源 https://blog.csdn.net/weixin_28475533/article/details/115996055?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522165095365616782395363619%2522%252C%2522scm%2522%253A%252220140713.130102334.pc%255Fall.%2522%257D&request_id=165095365616782395363619&biz_id=0&utm_medium=distribute.pc_search_result.none-task-blog-2~all~first_rank_ecpm_v1~rank_v31_ecpm-1-115996055.142^v9^control,157^v4^control&utm_term=matlab+%E8%AF%BB%E5%8F%96fit%E6%96%87%E4%BB%B6&spm=1018.2226.3001.4187
clc,clear,close all;

addpath('C:\Users\DELL\Desktop\dec_class_220426_141057-myopic');

filename = 'Rpsf_dec_class_141057_titanhe_153_IF_scaled.fits';

fr_id = fopen(filename,'r','s'); %FITfile

% info=FITinfo(filename);
info=fitsinfo(filename);

%FIT file head info, which includeprimarydata and image

oset=info.PrimaryData.Offset;

%thesize of primarydata head and it is also the start of primarydata

cols =info.PrimaryData.Size(1);%the columns of primary data

rows =info.PrimaryData.Size(2);%the rows of primary data

SIZE_UNIT = 2880; %The size of a FIT logical record

SIZE_TYPE_NAME= 8; % size of keyword, Keyword= Value/Comment

SIZE_TYPE = 80; % size of (Keyword= Value/Comment)

num_type_prim = oset / SIZE_TYPE;%the numbers of primarydata head type

%-------get the start position oflamda and its step(log10)

%%%%% 这一段出错，注释掉
% for i=1:num_type_prim
%     fseek(fr_id, (i-1)*SIZE_TYPE,-1);
%     type_name = fread(fr_id, SIZE_TYPE_NAME,'*char');
%     switch type_name'
%         case 'BITPIX ' %data dispersion
%             fseek(fr_id,2,0); %skip two bytes,cause there are '= ' after type_name
%             str_tmp = fread(fr_id, SIZE_TYPE -SIZE_TYPE_NAME - 2,'*char');
%             data_percision =str2num(str_tmp(1:20)');
%             data_percision = abs(data_percision);
%         case 'CRVAL1 ' %thestart position of lamda
%             fseek(fr_id,2,0);
%             str_tmp = fread(fr_id, SIZE_TYPE -SIZE_TYPE_NAME - 2,'*char');
%             lamda_start =str2num(str_tmp(1:20)');
%         case 'CD1_1 ' %the step of lamda
%             fseek(fr_id,2,0);
%             str_tmp = fread(fr_id, SIZE_TYPE -SIZE_TYPE_NAME - 2,'*char');
%             lamda_deta =str2num(str_tmp(1:20)');
%     end
% end

% calc the array of lamda(not logformat)

% lamda = lamda_start:lamda_deta:lamda_start+lamda_deta*(cols-1);%%%%出错

% for i=1:length(lamda)
%     
%     lamda(i) = 10^lamda(i);
%     
% end

%--------------------write primaryhead----------------------

fph_id =fopen('primary_head_info.txt','w','b');

fseek(fr_id,0,-1);

for i=1:num_type_prim
    
    head_info = fread(fr_id,SIZE_TYPE);
    
    fwrite(fph_id,head_info);
    
    fprintf(fph_id,'\r\n');
    
end

fclose(fph_id);

%--------------------write primarydata----------------------

fpd_id =fopen('prim_data.txt','w','s');

fseek(fr_id,oset,-1);

lamda=linspace(1,512,512);%%%%自己加的
for i=1:5 %because of the openspeed in txt, this only write five rows data into txt file.If you want all, youcan replace 5 with rows
    
    data_prim = fread(fr_id, cols, 'double');
    
    figure
    
    plot(lamda, data_prim);
    
    fprintf(fpd_id,'%f ',data_prim);
    
end

fclose(fpd_id);

